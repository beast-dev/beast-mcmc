/*
 * GMRFSkygridLikelihood.java
 *
 * Copyright (c) 2002-2020 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.coalescent;

import dr.evolution.coalescent.TreeIntervalList;
import dr.evolution.coalescent.TreeIntervals;
import dr.evomodelxml.coalescent.GMRFSkyrideLikelihoodParser;
import dr.inference.model.*;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Arrays;
import java.util.List;

/**
 * @author Xiang Ji
 * @author Filippo Monti
 * @author Marc A. Suchard
 */

public class MultilocusNonparametricCoalescentLikelihood extends AbstractModelLikelihood implements Citable {

    private final int numGridPoints;

    private double[] numCoalEvents;
    private double[] storedNumCoalEvents;
    private double[] sufficientStatistics;
    private double[] storedSufficientStatistics;

    private final Parameter ploidyFactors;

    private double[] ploidySums;
    private double[] storedPloidySums;

    private final List<TreeIntervals> intervalsList;

    private final Parameter logPopSizes;
    private final Parameter gridPoints;

    private double logLikelihood;
    private double storedLogLikelihood;

    private boolean intervalsKnown;
    private boolean likelihoodKnown;

    public MultilocusNonparametricCoalescentLikelihood(List<TreeIntervals> intervalLists,
                                                       Parameter logPopSizes,
                                                       Parameter gridPoints,
                                                       Parameter ploidyFactors) {

        super(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD);

        // adding the key word to the the model means the keyword will be logged in the
        // header of the logfile.
        this.addKeyword("skygrid");
        if (intervalLists.size() > 1) {
            this.addKeyword("multilocus");
        }

        this.intervalsList = intervalLists;
        this.logPopSizes = logPopSizes;
        this.gridPoints = gridPoints;
        this.ploidyFactors = ploidyFactors;
        this.numGridPoints = gridPoints.getDimension();

        addVariable(logPopSizes);
        addVariable(gridPoints);
        addVariable(ploidyFactors);

        if (ploidyFactors.getDimension() != intervalLists.size()) {
            throw new IllegalArgumentException("Ploidy factors parameter should have length " + intervalLists.size());
        }

        int fieldLength = logPopSizes.getDimension();

        this.sufficientStatistics = new double[fieldLength];
        this.numCoalEvents = new double[fieldLength];
        this.storedNumCoalEvents = new double[fieldLength];
        this.ploidySums = new double[fieldLength];
        this.storedPloidySums = new double[fieldLength];
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

        if (model instanceof TreeIntervalList) {
            TreeIntervalList treeModel = (TreeIntervalList) model;
            int tn = intervalsList.indexOf(treeModel);
            if (tn >= 0) {
                intervalsKnown = false;
                likelihoodKnown = false;
            } else {
                throw new RuntimeException("Unknown tree");
            }
        } else {
            throw new RuntimeException("Unknown object");
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    private int moveToNextTimeIndex(int treeIndex, int lastTimeIndex, double[] times) {
        int currentTimeIndex = lastTimeIndex;
        double currentTime = intervalsList.get(treeIndex).getIntervalTime(currentTimeIndex);
        double nextTime = intervalsList.get(treeIndex).getIntervalTime(currentTimeIndex + 1);
        while (nextTime <= currentTime && currentTimeIndex + 2 < intervalsList.get(treeIndex).getIntervalCount()) {
            currentTimeIndex++;
            currentTime = intervalsList.get(treeIndex).getIntervalTime(currentTimeIndex);
            nextTime = intervalsList.get(treeIndex).getIntervalTime(currentTimeIndex + 1);
        }
        times[0] = currentTime;
        times[1] = nextTime;
        return currentTimeIndex;
    }

    protected void setupSufficientStatistics() {

        Arrays.fill(numCoalEvents, 0);
        Arrays.fill(sufficientStatistics, 0.0);
        Arrays.fill(ploidySums, 0);

//        //index of the smallest grid point greater than at least one sampling/coalescent time in current tree
//        int minGridIndex;
//        //index of the greatest grid point less than at least one sampling/coalescent time in current tree
//        int maxGridIndex;

        double[] gridPoints = this.gridPoints.getParameterValues(); // TODO delegate to interface

        for (int i = 0; i < intervalsList.size(); i++) {

            double ploidyFactor = 1 / getPopulationFactor(i);

            double[] currentAndNextTime = new double[2];

            int currentTimeIndex = moveToNextTimeIndex(i, 0, currentAndNextTime);

            int numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
            int minGridIndex = 0;
            while (minGridIndex < numGridPoints - 1 && gridPoints[minGridIndex] <= currentAndNextTime[0]) { // MAS: Unclear about need for -1
                minGridIndex++;
            }
            int currentGridIndex = minGridIndex;

            double lastCoalescentTime = currentAndNextTime[0] + intervalsList.get(i).getTotalDuration();

            int maxGridIndex = numGridPoints - 1;
            while ((maxGridIndex >= 0) && (gridPoints[maxGridIndex] >= lastCoalescentTime)) {
                maxGridIndex = maxGridIndex - 1;
            }

            if (maxGridIndex >= 0 && minGridIndex < numGridPoints) {


                //from likelihood of interval between first sampling time and gridPoints[minGridIndex]

                while (currentAndNextTime[1] < gridPoints[currentGridIndex]) {

                    //check to see if interval ends with coalescent event
                    if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {

                        numCoalEvents[currentGridIndex]++;
                    }
                    sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (currentAndNextTime[1] - currentAndNextTime[0]) * numLineages * (numLineages - 1) * 0.5 * ploidyFactor;
                    currentTimeIndex++;
                    currentTimeIndex = moveToNextTimeIndex(i, currentTimeIndex, currentAndNextTime);

                    numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);

                }

                sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (gridPoints[currentGridIndex] - currentAndNextTime[0]) * numLineages * (numLineages - 1) * 0.5 * ploidyFactor;
                ploidySums[currentGridIndex] = ploidySums[currentGridIndex] + Math.log(ploidyFactor) * numCoalEvents[currentGridIndex];

                currentGridIndex++;


                //from likelihood of intervals between gridPoints[minGridIndex] and gridPoints[maxGridIndex]

                while (currentGridIndex <= maxGridIndex) {
                    if (currentAndNextTime[1] >= gridPoints[currentGridIndex]) {
                        sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (gridPoints[currentGridIndex] - gridPoints[currentGridIndex - 1]) * numLineages * (numLineages - 1) * 0.5 * ploidyFactor;
                        ploidySums[currentGridIndex] = ploidySums[currentGridIndex] + Math.log(ploidyFactor) * numCoalEvents[currentGridIndex];

                        currentGridIndex++;
                    } else {

                        sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (currentAndNextTime[1] - gridPoints[currentGridIndex - 1]) * numLineages * (numLineages - 1) * 0.5 * ploidyFactor;

                        //check to see if interval ends with coalescent event
                        if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
                            numCoalEvents[currentGridIndex]++;
                        }
                        currentTimeIndex++;
                        currentTimeIndex = moveToNextTimeIndex(i, currentTimeIndex, currentAndNextTime);

                        numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);

                        while (currentAndNextTime[1] < gridPoints[currentGridIndex]) {
                            //check to see if interval is coalescent interval or sampling interval
                            if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
                                numCoalEvents[currentGridIndex]++;
                            }
                            sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (currentAndNextTime[1] - currentAndNextTime[0]) * numLineages * (numLineages - 1) * 0.5 * ploidyFactor;

                            currentTimeIndex++;
                            currentTimeIndex = moveToNextTimeIndex(i, currentTimeIndex, currentAndNextTime);

                            numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);

                        }
                        sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (gridPoints[currentGridIndex] - currentAndNextTime[0]) * numLineages * (numLineages - 1) * 0.5 * ploidyFactor;
                        ploidySums[currentGridIndex] = ploidySums[currentGridIndex] + Math.log(ploidyFactor) * numCoalEvents[currentGridIndex];

                        currentGridIndex++;
                    }
                }

                //from likelihood of interval between gridPoints[maxGridIndex] and lastCoalescentTime

                sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (currentAndNextTime[1] - gridPoints[currentGridIndex - 1]) * numLineages * (numLineages - 1) * 0.5 * ploidyFactor;

                //check to see if interval ends with coalescent event
                if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
                    numCoalEvents[currentGridIndex]++;
                }

                currentTimeIndex++;

                while ((currentTimeIndex + 1) < intervalsList.get(i).getIntervalCount()) {

                    currentTimeIndex = moveToNextTimeIndex(i, currentTimeIndex, currentAndNextTime);

                    numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);

                    //check to see if interval is coalescent interval or sampling interval

                    if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
                        numCoalEvents[currentGridIndex]++;
                    }
                    sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (currentAndNextTime[1] - currentAndNextTime[0]) * numLineages * (numLineages - 1) * 0.5 * ploidyFactor;
                    currentAndNextTime[0] = currentAndNextTime[1];
                    currentTimeIndex++;

                }

                // if tree does not overlap with any gridpoints/change-points, in which case logpopsize is constant

            } else {
                while ((currentTimeIndex + 1) < intervalsList.get(i).getIntervalCount()) {
                    //check to see if interval is coalescent interval or sampling interval
                    if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
                        numCoalEvents[currentGridIndex]++;
                    }
                    sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (currentAndNextTime[1] - currentAndNextTime[0]) * numLineages * (numLineages - 1) * 0.5 * ploidyFactor;

                    currentTimeIndex++;
                    if ((currentTimeIndex + 1) < intervalsList.get(i).getIntervalCount()) {
                        currentTimeIndex = moveToNextTimeIndex(i, currentTimeIndex, currentAndNextTime);

                        numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);

                    }

                }
                ploidySums[currentGridIndex] = ploidySums[currentGridIndex] + Math.log(ploidyFactor) * numCoalEvents[currentGridIndex];

            }
        }
    }

    private double calculateLogCoalescentLikelihood() {

        checkIntervals();

        double currentLike = 0;
        for (int i = 0; i < logPopSizes.getDimension(); i++) {
            double currentGamma = logPopSizes.getParameterValue(i);
            currentLike += -numCoalEvents[i] * currentGamma + ploidySums[i] - sufficientStatistics[i] * Math.exp(-currentGamma);
        }

        return currentLike;
    }


    @Override
    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogCoalescentLikelihood();
            likelihoodKnown = true;
        }

        return logLikelihood;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
    }

    private double getPopulationFactor(int nt) {
        return ploidyFactors.getParameterValue(nt);
    }

    protected void storeState() {
        System.arraycopy(numCoalEvents, 0, storedNumCoalEvents, 0, numCoalEvents.length);
        System.arraycopy(ploidySums, 0, storedPloidySums, 0, ploidySums.length);

        storedLogLikelihood = logLikelihood;
    }


    protected void restoreState() {
        // Swap pointers
        double[] tmp = numCoalEvents;
        numCoalEvents = storedNumCoalEvents;
        storedNumCoalEvents = tmp;

        double[] tmp2 = ploidySums;
        ploidySums = storedPloidySums;
        storedPloidySums = tmp2;

        logLikelihood = storedLogLikelihood;
    }

    @Override
    protected void acceptState() {

    }

//    private double[] getGradientLogDensity() {
//
//        checkIntervals();
//
//        final int dim = popSizeParameter.getSize();
//        double[] gradLogDens = new double[dim];
//        double[] gamma = getMeanAdjustedGamma();
//
//        double currentPrec = precisionParameter.getParameterValue(0);
//
//        gradLogDens[0] = -currentPrec * (gamma[0] - gamma[1])
//                - numCoalEvents[0] + sufficientStatistics[0]
//                * Math.exp(-popSizeParameter.getParameterValue(0));
//
//        for (int i = 1; i < (dim - 1); i++) {
//            gradLogDens[i] = -currentPrec * (-gamma[i - 1] + 2 * gamma[i] - gamma[i + 1])
//                    - numCoalEvents[i] + sufficientStatistics[i]
//                    * Math.exp(-popSizeParameter.getParameterValue(i));
//        }
//
//        gradLogDens[dim - 1] = -currentPrec * (gamma[dim - 1] - gamma[dim - 2])
//                - numCoalEvents[dim - 1] + sufficientStatistics[dim - 1]
//                * Math.exp(-popSizeParameter.getParameterValue(dim - 1));
//
//        return gradLogDens;
//
//    }

    private void checkIntervals() {
        if (!intervalsKnown) {
            setupSufficientStatistics();
            intervalsKnown = true;
        }
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "non-parametric coalescent";
    }

    @Override
    public List<Citation> getCitations() {
        return Arrays.asList(new Citation(
                        new Author[]{
                                new Author("MS", "Gill"),
                                new Author("P", "Lemey"),
                                new Author("NR", "Faria"),
                                new Author("A", "Rambaut"),
                                new Author("B", "Shapiro"),
                                new Author("MA", "Suchard")
                        },
                        "Improving Bayesian population dynamics inference: a coalescent-based model for multiple loci",
                        2013,
                        "Mol Biol Evol",
                        30, 713, 724
                ),
                new Citation(
                        new Author[]{
                                new Author("VN", "Minin"),
                                new Author("EW", "Bloomquist"),
                                new Author("MA", "Suchard")
                        },
                        "Smooth skyride through a rough skyline: Bayesian coalescent-based inference of population dynamics",
                        2008,
                        "Mol Biol Evol",
                        25, 1459, 1471,
                        "10.1093/molbev/msn090"
                )
        );
    }
}
