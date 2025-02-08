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
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.inference.model.*;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.Reportable;

import java.util.Arrays;
import java.util.List;

/**
 * @author Xiang Ji
 * @author Filippo Monti
 * @author Marc A. Suchard
 */

public class MultilocusNonparametricCoalescentLikelihood extends AbstractModelLikelihood implements Citable, Reportable {

    private final int numGridPoints;

    private int[] numCoalEvents;
    private int[] storedNumCoalEvents;
    private double[] sufficientStatistics;
    private double[] storedSufficientStatistics;

    private final Parameter ploidyFactors;

    private double[] ploidySums;
    private double[] storedPloidySums;

    private final List<BigFastTreeIntervals> intervalsList;

    private final Parameter logPopSizes;
    private final Parameter gridPoints;

    private double logLikelihood;
    private double storedLogLikelihood;

    private boolean intervalsKnown;
    private boolean storedIntervalsKnown;
    private boolean likelihoodKnown;

    double[] fullTimeLine;
    int[] gridIndices;
    int[] numLineages;

    public MultilocusNonparametricCoalescentLikelihood(List<BigFastTreeIntervals> intervalLists,
                                                       Parameter logPopSizes,
                                                       Parameter gridPoints,
                                                       Parameter ploidyFactors) {

        super("Multilocus Nonparametric Coalescent Likelihood");

        // adding the key word to the model means the keyword will be logged in the
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

        for (BigFastTreeIntervals intervals : intervalLists) {
            addModel(intervals);
        }

        if (ploidyFactors.getDimension() != intervalLists.size()) {
            throw new IllegalArgumentException("Ploidy factors parameter should have length " + intervalLists.size());
        }

        int fieldLength = logPopSizes.getDimension();

        this.sufficientStatistics = new double[fieldLength];
        this.storedSufficientStatistics = new double[fieldLength];
        this.numCoalEvents = new int[fieldLength];
        this.storedNumCoalEvents = new int[fieldLength];
        this.ploidySums = new double[fieldLength];
        this.storedPloidySums = new double[fieldLength];
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model instanceof BigFastTreeIntervals) {
            BigFastTreeIntervals treeModel = (BigFastTreeIntervals) model;
            int tn = intervalsList.indexOf(treeModel);
            if (tn >= 0) {
                intervalsKnown = false; // TODO This should only fire the change for one tree model
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
        likelihoodKnown = false;
    }

    private void setupSufficientStatistics() {
        Arrays.fill(sufficientStatistics, 0);
        Arrays.fill(ploidySums, 0.0);
        Arrays.fill(numCoalEvents, 0);

        for (int treeIndex = 0; treeIndex < intervalsList.size(); treeIndex++) {
            double ploidyFactor = 1 / getPopulationFactor(treeIndex);

            SingleTreeNodesTimeline singleTreeNodesTimeLine = new
                    SingleTreeNodesTimeline(intervalsList.get(treeIndex), gridPoints);
            fullTimeLine = singleTreeNodesTimeLine.getMergedTimeLine();
            numLineages = singleTreeNodesTimeLine.getMergedNumLineages();
            gridIndices = singleTreeNodesTimeLine.gridIndices;

            int[] tempNumCoalEvents = singleTreeNodesTimeLine.getNumCoalEvents();
            for (int j = 0; j < numCoalEvents.length; j++)  numCoalEvents[j] += tempNumCoalEvents[j];

            int i = 0; // index for the fullTimeLine
            while (gridIndices[i] == i) i++; // choose the first grid point bigger than the most recent sampling time

            boolean skipFirstSamplingTime = true; // necessary since "i" cannot be incremented here

            double interval = 0.0;
            for (int gridIndex = i; gridIndex < numGridPoints; gridIndex++) {
                while (i <= gridIndices[gridIndex]) {
                    if (!skipFirstSamplingTime) {
                            interval = (fullTimeLine[i] - fullTimeLine[i - 1]);
                            sufficientStatistics[gridIndex] += 0.5 * numLineages[i - 1] * (numLineages[i - 1] - 1) *
                                    interval * ploidyFactor;
                            ploidySums[gridIndex] += Math.log(ploidyFactor) * tempNumCoalEvents[gridIndex];
                    } else {
                        skipFirstSamplingTime = false;
                    }
                    i++;
                }
            }
            // manage events after last grid point
            if (i < fullTimeLine.length) {
                while (i < fullTimeLine.length) {
                    interval = (fullTimeLine[i] - fullTimeLine[i - 1]);
                    sufficientStatistics[numGridPoints] += 0.5 * numLineages[i - 1] * (numLineages[i - 1] - 1) *
                            interval * ploidyFactor;
                    ploidySums[numGridPoints] += Math.log(ploidyFactor) * tempNumCoalEvents[numGridPoints];
                    i++;
                }
            }
        }
    }



    // this builds a timeline with coalescent and sampling events for a single tree
    private class SingleTreeNodesTimeline {
        final int nNodes;
        final BigFastTreeIntervals treeIntervals;

        final double[] timeLine;
        double[] mergedTimeLine;

        final boolean[] flagCoalescentEvent;
//        boolean[] mergedFlagCoalescentEvent;

        int[] mergedNumLineages;
        int[] numCoalEvents;
        final Parameter gridPoints;
        int[] gridIndices;

        private SingleTreeNodesTimeline(BigFastTreeIntervals treeIntervals) {
            this(treeIntervals, null);
        }

        private SingleTreeNodesTimeline(BigFastTreeIntervals treeIntervals, Parameter gridPoints) {
            this.treeIntervals = treeIntervals;
            this.nNodes = treeIntervals.getIntervalCount() + 1;
            this.timeLine = new double[nNodes];
            this.flagCoalescentEvent = new boolean[nNodes];
            this.gridPoints = gridPoints;

            makeLine();

            if (gridPoints != null) { // "merged" =  nodes' times and grid points
                this.mergedTimeLine = new double[nNodes + gridPoints.getDimension()];
//                this.mergedFlagCoalescentEvent = new boolean[nNodes + gridPoints.getDimension()];
                this.gridIndices = new int[gridPoints.getDimension()];
                numCoalEvents = new int[gridPoints.getDimension() + 1]; // "+1" to account for the events after the last grid point
                mergedNumLineages = new int[nNodes + gridPoints.getDimension() + 1];
                integrateGridPoints(gridPoints.getParameterValues());
            }
        }

        private void makeLine() {
            timeLine[0] = treeIntervals.getStartTime();
            flagCoalescentEvent[0] = false;
            for (int nodeIndex = 1; nodeIndex < nNodes; nodeIndex++) {
                timeLine[nodeIndex] = treeIntervals.getIntervalTime(nodeIndex);
//                       timeLine[nodeIndex - 1];
                flagCoalescentEvent[nodeIndex] =
                        String.valueOf(treeIntervals.getIntervalType(nodeIndex - 1)).equals("coalescent"); //TODO this is hard coded ...
            }
        }

        private void computeNumLineages(int index, boolean isCoalescentEvent) {
            if (isCoalescentEvent) {
                mergedNumLineages[index] = mergedNumLineages[index - 1] - 1;
            } else {
                if (index == 0) { // the first node is a tip node
                    mergedNumLineages[index] = 1;
                } else {
                    mergedNumLineages[index] = mergedNumLineages[index - 1] + 1;
                }
            }
        }

        private void integrateGridPoints(double[] gridPointsVector) {
            int nGridPoints = gridPointsVector.length;
            Arrays.fill(mergedNumLineages, 0);
            Arrays.fill(numCoalEvents, 0);

            int i = 0, j = 0;
            while (i < nNodes && j < nGridPoints) {
                if (timeLine[i] <= gridPointsVector[j]) { // grid points are set after node times if equal
                    mergedTimeLine[i + j] = timeLine[i];
                    computeNumLineages(i + j, flagCoalescentEvent[i]);
                    if (flagCoalescentEvent[i]) numCoalEvents[j] += 1;
                    i++;
                } else {
                    mergedTimeLine[i + j] = gridPointsVector[j];
                    gridIndices[j] = i + j;
                    mergedNumLineages[i + j] = mergedNumLineages[i + j - 1];
                    j++;
                }
            }

            // Add remaining elements
            while (i < nNodes) {
                mergedTimeLine[i + j] = timeLine[i];
                computeNumLineages(i + j, flagCoalescentEvent[i]);
                if (flagCoalescentEvent[i]) numCoalEvents[j] += 1;
                i++;
            }
            while (j < nGridPoints) {
                mergedTimeLine[i + j] = gridPointsVector[j];
                gridIndices[j] = i + j;
                mergedNumLineages[i + j] = mergedNumLineages[i + j - 1];
                j++;
            }
        }
//        private boolean[] getFlagCoalescentEvent() { return flagCoalescentEvent; }
//        private double[] getTimeLine() { return timeLine; }
//        private boolean[] getmergedFlagCoalescentEvent() { return mergedFlagCoalescentEvent; }
        private double[] getMergedTimeLine() { return mergedTimeLine; }
        private int[] getMergedNumLineages() {return mergedNumLineages;}
        private int[] getNumCoalEvents() {return numCoalEvents;}
    }

    protected double calculateLogCoalescentLikelihood() {

        computeSufficientStatistics();

        double currentLike = 0;
        for (int i = 0; i < logPopSizes.getDimension(); i++) {
            double currentGamma = logPopSizes.getParameterValue(i);
            currentLike += - numCoalEvents[i] * currentGamma
                    + ploidySums[i]
                    - sufficientStatistics[i] * Math.exp(-currentGamma);
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
        intervalsKnown = false;
        likelihoodKnown = false;
    }

    private double[] getFullTimeLine() {return fullTimeLine;}

    private double getPopulationFactor(int nt) {return ploidyFactors.getParameterValue(nt);}

    protected int getNGridPoints() { return gridPoints.getDimension();}

    protected double getPopulationSize(int i) { return Math.exp(logPopSizes.getParameterValue(i)); }

    protected double getLogPopulationSize(int i) { return logPopSizes.getParameterValue(i); }

    protected Parameter getLogPopSizes() { return logPopSizes; }

    protected int getPopSizeDimension() { return logPopSizes.getDimension(); }

    protected void storeState() {
        System.arraycopy(numCoalEvents, 0, storedNumCoalEvents, 0, numCoalEvents.length);
        System.arraycopy(ploidySums, 0, storedPloidySums, 0, ploidySums.length);
        System.arraycopy(sufficientStatistics, 0, storedSufficientStatistics, 0,
                sufficientStatistics.length);

        storedIntervalsKnown = intervalsKnown;
        storedLogLikelihood = logLikelihood;
    }

    protected void restoreState() {
        // Swap pointers
        int[] tmp = numCoalEvents;
        numCoalEvents = storedNumCoalEvents;
        storedNumCoalEvents = tmp;

        double[] tmp2 = ploidySums;
        ploidySums = storedPloidySums;
        storedPloidySums = tmp2;

        double[] tmp3 = sufficientStatistics;
        sufficientStatistics = storedSufficientStatistics;
        storedSufficientStatistics = tmp3;

        intervalsKnown = storedIntervalsKnown;

        logLikelihood = storedLogLikelihood;
    }

    @Override
    protected void acceptState() {

    }

    @SuppressWarnings("unused")
    private double[] getGradientLogDensity() {
        return getGradientLogDensity(logPopSizes.getParameterValues());
    }

    public double[] getGradientLogDensity(Object x) {
        double[] field = (double[]) x;
        computeSufficientStatistics();

        final int dim = field.length;
        double[] gradLogDens = new double[dim];

        for (int i = 0; i < dim; ++i) {
            gradLogDens[i] = -numCoalEvents[i] + sufficientStatistics[i]
                    * Math.exp(-field[i]);
        }

        return gradLogDens;
    }

    private void computeSufficientStatistics() {
        if (!intervalsKnown) {
            setupSufficientStatistics();
            intervalsKnown = true;
        }
    }

    @Override
    public String getReport() {
        return "Non-parametric Coalescent LogLikelihood: " + getLogLikelihood() + "\n" +
                "Sufficient statistics: " + Arrays.toString(sufficientStatistics) + "\n" +
                "Ploidy factors: " + Arrays.toString(ploidyFactors.getParameterValues()) + "\n" +
                "Grid points: " + Arrays.toString(gridPoints.getParameterValues()) + "\n" +
                "Log population sizes: " + Arrays.toString(logPopSizes.getParameterValues()) + "\n" +
                "Full time line: " + Arrays.toString(fullTimeLine) + "\n";
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "Non-parametric coalescent";
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