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

import dr.evolution.coalescent.IntervalType;
import dr.evolution.coalescent.TreeIntervals;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.GMRFSkyrideLikelihoodParser;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.SymmTridiagMatrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Mandev Gill
 * @author Marc A. Suchard
 */

public class GMRFMultilocusSkyrideLikelihood extends OldGMRFSkyrideLikelihood
        implements MultiLociTreeSet, CoalescentIntervalProvider, Citable {

    public static final boolean DEBUG = false;

    private final double cutOff;
    private final int numGridPoints;
    private final int oldFieldLength;
    private final int numTrees;

    // number of coalescent events which occur in an interval with constant population size
    private double[] numCoalEvents;
    private double[] storedNumCoalEvents;
    private double[] gridPoints;

    // sortedPoints[i][0] is the time of the i-th grid point or sampling or coalescent event
    // sortedPoints[i][1] is 0 if the i-th point is a grid point, 1 if it's a sampling point, and 2 if it's a coalescent point
    // sortedPoints[i][2] is the number of lineages present in the interval starting at time sortedPoints[i][0]

    private final Parameter phiParameter;
    private final Parameter ploidyFactors;
    private double[] ploidySums;
    private double[] storedPloidySums;

//    protected SymmTridiagMatrix precMatrix;
//    protected SymmTridiagMatrix storedPrecMatrix;

    private final SkygridHelper skygridHelper;
    //    protected List<Parameter> missingCov;
    private final List<MatrixParameter> covariates;
    private final List<Parameter> beta;
    private final List<Parameter> delta;
    private final List<Parameter> covPrecParametersRecent;
    private final List<Parameter> covPrecParametersDistant;

    private List<SymmTridiagMatrix> weightMatricesForMissingCovRecent;
    private List<SymmTridiagMatrix> weightMatricesForMissingCovDistant;

    private int[] firstObservedIndex;
    private int[] lastObservedIndex;
    private int[] recIndices;
    private int[] distIndices;

    private double[] coalescentEventStatisticValues;

    private List<Tree> treeList;
    private List<TreeIntervals> intervalsList;

    public GMRFMultilocusSkyrideLikelihood(List<Tree> treeList,
                                           Parameter popParameter,
                                           Parameter groupParameter,
                                           Parameter precParameter,
                                           Parameter lambda,
                                           Parameter beta,
                                           MatrixParameter dMatrix,
                                           boolean timeAwareSmoothing,
                                           double cutOff,
                                           int numGridPoints,
                                           Parameter phi,
                                           Parameter ploidyFactorsParameter) {

        super(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD);

        // adding the key word to the the model means the keyword will be logged in the
        // header of the logfile.
        this.addKeyword("skygrid");
        if (treeList.size() > 1) {
            this.addKeyword("multilocus");
        }

        this.popSizeParameter = popParameter;
        this.groupSizeParameter = groupParameter;
        this.precisionParameter = precParameter;
        this.lambdaParameter = lambda;
        this.betaParameter = beta;
        this.dMatrix = dMatrix;
        if (dMatrix != null) {
            addVariable(dMatrix);
        }
        this.timeAwareSmoothing = timeAwareSmoothing;

        this.cutOff = cutOff;
        this.numGridPoints = numGridPoints;
        this.phiParameter = phi;
        this.ploidyFactors = ploidyFactorsParameter;

        setupGridPoints();

        addVariable(popSizeParameter);
        addVariable(precisionParameter);
        addVariable(lambdaParameter);
        if (betaParameter != null) {
            addVariable(betaParameter);
            skygridHelper = new SkygridCovariateHelper();
        } else {
            skygridHelper = new SkygridHelper();
        }
        if (phiParameter != null) {
            addVariable(phiParameter);
        }
        addVariable(ploidyFactors);


        this.numTrees = setTree(treeList);

        int correctFieldLength = getCorrectFieldLength();

        if (popSizeParameter.getDimension() <= 1) {
            // popSize dimension hasn't been set yet, set it here:
            popSizeParameter.setDimension(correctFieldLength);
        }

        fieldLength = popSizeParameter.getDimension();
        if (correctFieldLength != fieldLength) {
            throw new IllegalArgumentException("Population size parameter should have length " + correctFieldLength);
        }

        oldFieldLength = getCorrectOldFieldLength();


        if (ploidyFactors.getDimension() != treeList.size()) {
            throw new IllegalArgumentException("Ploidy factors parameter should have length " + treeList.size());
        }


        // Field length must be set by this point
        wrapSetupIntervals();

        coalescentIntervals = new double[oldFieldLength];
        storedCoalescentIntervals = new double[oldFieldLength];
        sufficientStatistics = new double[fieldLength];
        storedSufficientStatistics = new double[fieldLength];

        numCoalEvents = new double[fieldLength];
        storedNumCoalEvents = new double[fieldLength];
        ploidySums = new double[fieldLength];
        storedPloidySums = new double[fieldLength];

        setupGMRFWeights();
        setupSufficientStatistics();

        addStatistic(new DeltaStatistic());

        initializationReport();

        /* Force all entries in groupSizeParameter = 1 for compatibility with Tracer */
        if (groupSizeParameter != null) {
            for (int i = 0; i < groupSizeParameter.getDimension(); i++)
                groupSizeParameter.setParameterValue(i, 1.0);
        }

        this.coalescentEventStatisticValues = new double[getNumberOfCoalescentEvents()];

        this.covariates = null;
        this.beta = null;
        this.delta = null;
        this.covPrecParametersRecent = null;
        this.covPrecParametersDistant = null;
    }


    //rewrite this constructor without duplicating so much code
    public GMRFMultilocusSkyrideLikelihood(List<Tree> treeList,
                                           Parameter popParameter,
                                           Parameter groupParameter,
                                           Parameter precParameter,
                                           Parameter lambda,
                                           Parameter betaParameter,
                                           MatrixParameter dMatrix,
                                           boolean timeAwareSmoothing,
                                           Parameter specGridPoints,
                                           List<MatrixParameter> covariates,
                                           Parameter ploidyFactorsParameter,
                                           List<Parameter> firstObservedIndexParameter,
                                           List<Parameter> lastObservedIndexParameter,
                                           List<Parameter> covPrecParametersRecent,
                                           List<Parameter> covPrecParametersDistant,
                                           Parameter recentIndices,
                                           Parameter distantIndices,
                                           List<Parameter> betaList,
                                           List<Parameter> deltaList) {

        super(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD);

        // adding the key word to the the model means the keyword will be logged in the
        // header of the logfile.
        this.addKeyword("skygrid");
        if (treeList.size() > 1) {
            this.addKeyword("multilocus");
        }

        this.gridPoints = specGridPoints.getParameterValues();
        this.numGridPoints = gridPoints.length;
        this.cutOff = gridPoints[numGridPoints - 1];

        if (firstObservedIndexParameter != null) {
            this.firstObservedIndex = new int[firstObservedIndexParameter.size()];
            for (int i = 0; i < firstObservedIndexParameter.size(); i++) {
                this.firstObservedIndex[i] = (int) firstObservedIndexParameter.get(i).getParameterValue(0);
            }

            if (recentIndices != null) {
                // indices specify which covariates require default unobserved covariate data prior
                this.recIndices = new int[firstObservedIndexParameter.size()];
                for (int i = 0; i < firstObservedIndexParameter.size(); i++) {
                    this.recIndices[i] = (int) recentIndices.getParameterValue(i);
                }
            } else {
                // If specific covariates not specified by indices, need default unobserved covariate data prior for all covariates
                this.recIndices = new int[firstObservedIndexParameter.size()];
                for (int i = 0; i < firstObservedIndexParameter.size(); i++) {
                    this.recIndices[i] = i + 1;
                }
            }

        }

        if (lastObservedIndexParameter != null) {
            this.lastObservedIndex = new int[lastObservedIndexParameter.size()];
            for (int i = 0; i < lastObservedIndexParameter.size(); i++) {
                this.lastObservedIndex[i] = (int) lastObservedIndexParameter.get(i).getParameterValue(0);
            }

            if (distantIndices != null) {
                // indices specify which covariates require default unobserved covariate data prior
                this.distIndices = new int[lastObservedIndexParameter.size()];
                for (int i = 0; i < lastObservedIndexParameter.size(); i++) {
                    this.distIndices[i] = (int) distantIndices.getParameterValue(i);
                }
            } else {
                // If specific covariates not specified by indices, need default unobserved covariate data prior for all covariates
                this.distIndices = new int[lastObservedIndexParameter.size()];
                for (int i = 0; i < lastObservedIndexParameter.size(); i++) {
                    this.distIndices[i] = i + 1;
                }
            }

        }

        /*else{
            for(int i=0; i < beta.getDimension(); i++) {
                this.lastObservedIndex[i] = popParameter.getDimension();
            }
        }*/

        this.betaParameter = betaParameter;
        if (betaParameter != null) {
            addVariable(betaParameter);
        }

        this.popSizeParameter = popParameter;
        this.groupSizeParameter = groupParameter;
        this.precisionParameter = precParameter;
        this.lambdaParameter = lambda;
        this.beta = betaList;

        if (deltaList != null) {
            this.delta = deltaList;
        } else {
            this.delta = new ArrayList<Parameter>();
            if (betaList != null) {
                for (int i = 0; i < betaList.size(); i++) {
                    Parameter deltaParam = new Parameter.Default(1.0);
                    deltaParam.setParameterValue(0, 1.0);
                    this.delta.add(deltaParam);
                }
            }
        }

        this.dMatrix = dMatrix;
        if (dMatrix != null) {
            addVariable(dMatrix);
        }
        this.timeAwareSmoothing = timeAwareSmoothing;
        this.ploidyFactors = ploidyFactorsParameter;
        this.covariates = covariates;
        if (covariates != null) {
            for (MatrixParameter cov : covariates) {
                addVariable(cov);
            }
        }
        this.covPrecParametersRecent = covPrecParametersRecent;
        if (covPrecParametersRecent != null) {
            for (Parameter covPrecRecent : covPrecParametersRecent) {
                addVariable(covPrecRecent);
            }
        }
        this.covPrecParametersDistant = covPrecParametersDistant;
        if (covPrecParametersDistant != null) {
            for (Parameter covPrecDistant : covPrecParametersDistant) {
                addVariable(covPrecDistant);
            }
        }

        addVariable(popSizeParameter);
        addVariable(precisionParameter);
        addVariable(lambdaParameter);

        addVariable(ploidyFactors);

        this.numTrees = setTree(treeList);

        int correctFieldLength = getCorrectFieldLength();

        if (popSizeParameter.getDimension() <= 1) {
            // popSize dimension hasn't been set yet, set it here:
            popSizeParameter.setDimension(correctFieldLength);
        }

        fieldLength = popSizeParameter.getDimension();
        if (correctFieldLength != fieldLength) {
            throw new IllegalArgumentException("Population size parameter should have length " + correctFieldLength);
        }

        oldFieldLength = getCorrectOldFieldLength();

        if (ploidyFactors.getDimension() != treeList.size()) {
            throw new IllegalArgumentException("Ploidy factor parameter should have length " + treeList.size());
        }

        // Field length must be set by this point

        if (betaList != null || betaParameter != null) {
            if (betaList != null) {
                for (Parameter betaParam : betaList) {
                    addVariable(betaParam);
                }
            }
            if (deltaList != null) {
                for (Parameter dParam : deltaList) {
                    addVariable(dParam);
                }
            }
            if (lastObservedIndexParameter != null || firstObservedIndexParameter != null) {
                setupGMRFWeightsForMissingCov();
                skygridHelper = new SkygridMissingCovariateHelper();
            } else {
                skygridHelper = new SkygridCovariateHelper();
            }
        } else {
            skygridHelper = new SkygridHelper();
        }

        wrapSetupIntervals();
        coalescentIntervals = new double[oldFieldLength];
        storedCoalescentIntervals = new double[oldFieldLength];
        sufficientStatistics = new double[fieldLength];
        storedSufficientStatistics = new double[fieldLength];

        numCoalEvents = new double[fieldLength];
        storedNumCoalEvents = new double[fieldLength];
        ploidySums = new double[fieldLength];
        storedPloidySums = new double[fieldLength];

        setupGMRFWeights();

        addStatistic(new DeltaStatistic());

        initializationReport();

        phiParameter = null;

        this.coalescentEventStatisticValues = new double[getNumberOfCoalescentEvents()];
    }

    protected int setTree(List<Tree> treeList) {
        treesSet = this;
        this.treeList = treeList;
        makeTreeIntervalList(treeList, true);
        return treeList.size();
    }

    private void makeTreeIntervalList(List<Tree> treeList, boolean add) {
        if (intervalsList == null) {
            intervalsList = new ArrayList<>();
        } else {
            intervalsList.clear();
        }
        for (Tree tree : treeList) {
            intervalsList.add(new TreeIntervals(tree));
            if (add && tree instanceof TreeModel) {
                addModel((TreeModel) tree);
            }
        }
    }

    protected int getCorrectFieldLength() {

        return numGridPoints + 1;
    }

    private int getCorrectOldFieldLength() {
        int tips = 0;
        for (Tree tree : treeList) {
            tips += tree.getExternalNodeCount();
        }
        return tips - treeList.size();
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

        if (model instanceof TreeModel) {
            TreeModel treeModel = (TreeModel) model;
            int tn = treeList.indexOf(treeModel);
            if (tn >= 0) {
                //   intervalsList.get(tn).setIntervalsUnknown();  // TODO Why is this slower (?) than remaking whole list?
                makeTreeIntervalList(treeList, false);
                intervalsKnown = false;
                likelihoodKnown = false;
            } else {
                throw new RuntimeException("Unknown tree modified in GMRFSkygridLikelihood");
            }
        } else {
            throw new RuntimeException("Unknown object modified in GMRFSkygridLikelihood");
        }
    }

    public void initializationReport() {
        System.out.println("Creating a GMRF smoothed skyride model for multiple loci (SkyGrid)");
        System.out.println("\tPopulation sizes: " + popSizeParameter.getDimension());
    }

    public void wrapSetupIntervals() {
        // Do nothing
    }

    private void setupGridPoints() {
        if (gridPoints == null) {
            gridPoints = new double[numGridPoints];
        } else {
            Arrays.fill(gridPoints, 0);
        }

        for (int pt = 0; pt < numGridPoints; pt++) {
            gridPoints[pt] = (pt + 1) * (cutOff / numGridPoints);
        }
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

//    private class CurrentAndNextTime {
//        private double currentTime;
//        private double nextTime;
//
//        public double getCurrentTime() {
//            return currentTime;
//        }
//
//        public double getNextTime() {
//            return nextTime;
//        }
//
//        public void setCurrentTime(int treeIndex, int timeIndex) {
//            currentTime = getTime(treeIndex, timeIndex);
//        }
//
//        private double getTime(int treeIndex, int timeIndex) {
//            return intervalsList.get(treeIndex).getIntervalTime(timeIndex);
//        }
//
//        public void setNextTime(int treeIndex, int timeIndex) {
//            nextTime = getTime(treeIndex, timeIndex);
//        }
//    }

    protected void setupSufficientStatistics() {

        Arrays.fill(numCoalEvents, 0);
        Arrays.fill(sufficientStatistics, 0.0);
        Arrays.fill(ploidySums, 0);
        //index of smallest grid point greater than at least one sampling/coalescent time in current tree
        int minGridIndex;
        //index of greatest grid point less than at least one sampling/coalescent time in current tree
        int maxGridIndex;

        int numLineages;

        int currentGridIndex;
        int currentTimeIndex;

//        double currentTime;
//        double nextTime;
        double[] currentAndNextTime = new double[2];
        double ploidyFactor;

        //time of last coalescent event in tree
        double lastCoalescentTime;

        for (int i = 0; i < numTrees; i++) {
            ploidyFactor = 1 / getPopulationFactor(i);
            currentTimeIndex = moveToNextTimeIndex(i, 0, currentAndNextTime);

            numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
            minGridIndex = 0;
            while (minGridIndex < numGridPoints - 1 && gridPoints[minGridIndex] <= currentAndNextTime[0]) { // MAS: Unclear about need for -1
                minGridIndex++;
            }
            currentGridIndex = minGridIndex;

            lastCoalescentTime = currentAndNextTime[0] + intervalsList.get(i).getTotalDuration();

//            theLastTime = lastCoalescentTime;

            maxGridIndex = numGridPoints - 1;
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

    public double[] getNumCoalEvents() {
        return numCoalEvents;
    }

    public int getNumberOfCoalescentEvents() {
        return getCorrectOldFieldLength();
    }

    public double getCoalescentEventsStatisticValue(int i) {
        if (i == 0) {

            if (DEBUG) {
                System.err.println("numTrees: " + numTrees);
                System.err.println("getCoalescentIntervalDimension(): " + super.getCoalescentIntervalDimension());
                System.err.println("getNumberOfCoalescentEvents(): " + getNumberOfCoalescentEvents());
                System.err.println("getIntervalCount(): " + getIntervalCount());
                System.err.println("intervalsList.size(): " + intervalsList.size());
                System.err.println("intervalsList.get(0).getIntervalCount(): " + intervalsList.get(0).getIntervalCount());
            }

            if (numTrees > 1) {
                throw new RuntimeException("Generalized stepping-stone sampling for the Skygrid not implemented for #trees > 1");
            }
            for (int j = 0; j < coalescentEventStatisticValues.length; j++) {
                coalescentEventStatisticValues[j] = 0.0;
            }
            int counter = 0;

            for (int j = 0; j < intervalsList.get(0).getIntervalCount(); j++) {
                if (intervalsList.get(0).getIntervalType(j) == IntervalType.COALESCENT) {
                    //this.coalescentEventStatisticValues[counter] += getCoalescentInterval(j) * (getLineageCount(j) * (getLineageCount(j) - 1.0)) / 2.0;
                    this.coalescentEventStatisticValues[counter] += intervalsList.get(0).getInterval(j) * (intervalsList.get(0).getLineageCount(j) * (intervalsList.get(0).getLineageCount(j) - 1.0)) / 2.0;
                    counter++;
                } else {
                    //this.coalescentEventStatisticValues[counter] += getCoalescentInterval(j) * (getLineageCount(j) * (getLineageCount(j) - 1.0)) / 2.0;
                    this.coalescentEventStatisticValues[counter] += intervalsList.get(0).getInterval(j) * (intervalsList.get(0).getLineageCount(j) * (intervalsList.get(0).getLineageCount(j) - 1.0)) / 2.0;
                }
            }
        }
        return coalescentEventStatisticValues[i];
        //throw new RuntimeException("getCoalescentEventsStatisticValue(int i) not implemented for Bayesian Skygrid");
        //return sufficientStatistics[i];
    }

    protected double calculateLogCoalescentLikelihood() {

        checkIntervals();

        // Matrix operations taken from block update sampler to calculate data likelihood and field prior

        double currentLike = 0;
        double[] currentGamma = popSizeParameter.getParameterValues();

        for (int i = 0; i < fieldLength; i++) {
            currentLike += -numCoalEvents[i] * currentGamma[i] + ploidySums[i] - sufficientStatistics[i] * Math.exp(-currentGamma[i]);
        }

        return currentLike;
    }


    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogCoalescentLikelihood();
            logFieldLikelihood = skygridHelper.getLogFieldLikelihood();
            likelihoodKnown = true;
        }

        return logLikelihood + logFieldLikelihood;
    }

    protected void setupGMRFWeights() {

        //setupSufficientStatistics();

        //Set up the weight Matrix
        double[] offdiag = new double[fieldLength - 1];
        double[] diag = new double[fieldLength];

        //    private double theLastTime;
        double diagonalValue = 2;
        //First set up the offdiagonal entries;

        for (int i = 0; i < fieldLength - 1; i++) {
            offdiag[i] = -1;
        }

        //Then set up the diagonal entries;
        for (int i = 1; i < fieldLength - 1; i++) {
            //	diag[i] = -(offdiag[i] + offdiag[i - 1]);
            diag[i] = diagonalValue;
        }
        //Take care of the endpoints
        //diag[0] = -offdiag[0];
        //diag[fieldLength - 1] = -offdiag[fieldLength - 2];
        diag[0] = diagonalValue - 1.0;
        diag[fieldLength - 1] = diagonalValue - 1.0;


        weightMatrix = new SymmTridiagMatrix(diag, offdiag);

    }

    protected double getFieldScalar() {
        return 1.0;
    }

    private void setupGMRFWeightsForMissingCov() {

        if (firstObservedIndex != null) {
            weightMatricesForMissingCovRecent = new ArrayList<>();

            for (int i = 0; i < covPrecParametersRecent.size(); i++) {
                double[] offdiagRec = new double[firstObservedIndex[i] - 2];
                double[] diagRec = new double[firstObservedIndex[i] - 1];

                for (int k = 0; k < firstObservedIndex[i] - 2; k++) {
                    offdiagRec[k] = -1;
                }

                for (int k = 1; k < firstObservedIndex[i] - 1; k++) {
                    diagRec[k] = 2.0;
                }
                diagRec[0] = 1.0;

                weightMatricesForMissingCovRecent.add(i, new SymmTridiagMatrix(diagRec, offdiagRec));
            }

        }

        if (lastObservedIndex != null) {
            weightMatricesForMissingCovDistant = new ArrayList<>();

            for (int i = 0; i < covPrecParametersDistant.size(); i++) {
                double[] offdiag = new double[fieldLength - lastObservedIndex[i] - 1];
                double[] diag = new double[fieldLength - lastObservedIndex[i]];

                //First set up the offdiagonal entries;

                for (int k = 0; k < fieldLength - lastObservedIndex[i] - 1; k++) {
                    offdiag[k] = -1;
                }

                //Then set up the diagonal entries;
                for (int k = 0; k < fieldLength - lastObservedIndex[i] - 1; k++) {
                    //	diag[i] = -(offdiag[i] + offdiag[i - 1]);
                    diag[k] = 2.0;
                }
                //Take care of the endpoint
                diag[fieldLength - lastObservedIndex[i] - 1] = 1.0;

                weightMatricesForMissingCovDistant.add(i, new SymmTridiagMatrix(diag, offdiag));
            }
        }

    }

    private SymmTridiagMatrix getScaledWeightMatrixForMissingCovRecent(double precision, int covIndex, int firstObs) {
        SymmTridiagMatrix a = weightMatricesForMissingCovRecent.get(covIndex).copy();
        for (int i = 0; i < a.numRows() - 1; i++) {
            a.set(i, i, a.get(i, i) * precision);
            a.set(i + 1, i, a.get(i + 1, i) * precision);
        }
        a.set(firstObs - 2, firstObs - 2,
                a.get(firstObs - 2, firstObs - 2) * precision);
        return a;
    }

    private SymmTridiagMatrix getScaledWeightMatrixForMissingCovDistant(double precision, int covIndex, int lastObs) {
        SymmTridiagMatrix a = weightMatricesForMissingCovDistant.get(covIndex).copy();
        for (int i = 0; i < a.numRows() - 1; i++) {
            a.set(i, i, a.get(i, i) * precision);
            a.set(i + 1, i, a.get(i + 1, i) * precision);
        }
        a.set(fieldLength - lastObs - 1, fieldLength - lastObs - 1,
                a.get(fieldLength - lastObs - 1, fieldLength - lastObs - 1) * precision);
        return a;
    }

    public int nLoci() {
        return treeList.size();
    }

    public Tree getTree(int nt) {
        return treeList.get(nt);
    }

    public TreeIntervals getTreeIntervals(int nt) {
        return intervalsList.get(nt);
    }

    public double getPopulationFactor(int nt) {
        return ploidyFactors.getParameterValue(nt);
    }

    @Deprecated
    public List<Parameter> getBetaListParameter() {
        return beta;
    }

    public List<MatrixParameter> getCovariates() {
        return covariates;
    }

    public void storeTheState() {
        for (TreeIntervals intervals : intervalsList) {
            intervals.storeState();
        }
    }

    public void restoreTheState() {
        for (TreeIntervals intervals : intervalsList) {
            intervals.restoreState();
        }
    }

    protected void storeState() {
        // System.arraycopy(numCoalEvents, 0, storedNumCoalEvents, 0, numCoalEvents.length);
        super.storeState();
        System.arraycopy(numCoalEvents, 0, storedNumCoalEvents, 0, numCoalEvents.length);
        // storedPrecMatrix = precMatrix.copy();
        System.arraycopy(ploidySums, 0, storedPloidySums, 0, ploidySums.length);
    }


    protected void restoreState() {
        super.restoreState();

        // Swap pointers
        double[] tmp = numCoalEvents;
        numCoalEvents = storedNumCoalEvents;
        storedNumCoalEvents = tmp;
        double[] tmp2 = ploidySums;
        ploidySums = storedPloidySums;
        storedPloidySums = tmp2;
    }

    // Implementation of GradientWrtParameterProvider
    // Need to update to include beta and log-transformed precision parameter

    public Parameter getParameter() {
        return popSizeParameter;
    }

    public int getDimension() {
        return getParameter().getDimension();
    }

    public Likelihood getLikelihood() {
        return this;
    }

    public double[] getGradientWrtLogPopulationSize() {
        return getGradientLogDensity();
    }

    public double[] getDiagonalHessianWrtLogPopulationSize() {
        return getDiagonalHessianLogDensity();
    }

    private double[] getMeanAdjustedGamma() {
        return skygridHelper.getMeanAdjustedGamma().getData();
    }

    public double[] getGradientWrtPrecision() {

        double[] gamma = getMeanAdjustedGamma();
        double currentPrec = precisionParameter.getParameterValue(0);

        double grad = numGridPoints / (2 * currentPrec);
        for (int i = 0; i < numGridPoints; i++) {
            grad -= 0.5 * (gamma[i + 1] - gamma[i]) * (gamma[i + 1] - gamma[i]);
        }

        return new double[]{grad};
    }

    public double[] getDiagonalHessianWrtPrecision() {

        double currentPrec = precisionParameter.getParameterValue(0);
        double hessian = -numGridPoints / (2 * currentPrec * currentPrec);

        return new double[]{hessian};
    }

    public double[] getGradientWrtRegressionCoefficients() {

        if (this.beta == null) return null;

        if (this.beta.size() > 1 || covariates.size() > 1) {
            throw new RuntimeException("This is not the way to handle multidimensional parameters");
        }

        Parameter beta = this.beta.get(0);
        MatrixParameter covk = covariates.get(0);
        boolean transposed = isTransposed(numGridPoints + 1, beta.getDimension(), covk);

        // TODO I believe we need one Parameter beta (is the same across loci) and List covariates (can differ across loci)

        // TODO Need to delegate to SkygridCovariateHelper

        double[] gradLogDens = new double[beta.getDimension()];
        double[] gamma = getMeanAdjustedGamma();

        double currentPrec = precisionParameter.getParameterValue(0);

        for (int k = 0; k < beta.getDimension(); k++) {

            double gradient = // numGridPoints / 2
                    +currentPrec * (gamma[0] - gamma[1]) * getCovariateValue(covk, 0, k, transposed) //covk.getParameterValue(k,0)
                            + currentPrec * (gamma[numGridPoints] - gamma[numGridPoints - 1]) * getCovariateValue(covk, numGridPoints, k, transposed) //covk.getParameterValue(k, numGridPoints)
//                    - 0.5 * currentPrec * (gamma[1] - gamma[0]) * (gamma[1] - gamma[0])
                    ;

            for (int i = 1; i < numGridPoints; i++) {
                gradient +=
//                        -0.5 * currentPrec * (gamma[i + 1] - gamma[i]) * (gamma[i + 1] - gamma[i])
                        +currentPrec * (-gamma[i - 1] + 2 * gamma[i] - gamma[i + 1]) * getCovariateValue(covk, i, k, transposed); //covk.getParameterValue(k, i);
            }

            gradLogDens[k] = gradient;
        }

        return gradLogDens;
    }

    public double[] getDiagonalHessianWrtRegressionCoefficients() {

        if (this.beta == null) return null;

        if (this.beta.size() > 1 || covariates.size() > 1) {
            throw new RuntimeException("This is not the way to handle multidimensional parameters");
        }

        Parameter beta = this.beta.get(0);
        MatrixParameter covk = covariates.get(0);
        boolean transposed = isTransposed(numGridPoints + 1, beta.getDimension(), covk);

        double[] hessian = new double[beta.getDimension()];
        double currentPrec = precisionParameter.getParameterValue(0);

        for (int k = 0; k < beta.getDimension(); k++) {

            double h = // numGridPoints / 2
                    -currentPrec * (
                            getCovariateValue(covk, 0, k, transposed) //covk.getParameterValue(k, 0)
                                    - getCovariateValue(covk, 1, k, transposed)) //covk.getParameterValue(k,1))
                            * getCovariateValue(covk, 0, k, transposed) //covk.getParameterValue(k,0)
                            - currentPrec * (
                            getCovariateValue(covk, numGridPoints, k, transposed) //covk.getParameterValue(k, numGridPoints)
                                    - getCovariateValue(covk, numGridPoints - 1, k, transposed)) //covk.getParameterValue(k, numGridPoints - 1))
                            * getCovariateValue(covk, numGridPoints, k, transposed) //covk.getParameterValue(k, numGridPoints)
                    ;

            for (int i = 1; i < numGridPoints; i++) {
                h += -currentPrec * (
                        -getCovariateValue(covk, i - 1, k, transposed) //covk.getParameterValue(k, i - 1)
                                + 2 * getCovariateValue(covk, i, k, transposed) //covk.getParameterValue(k, i)
                                - getCovariateValue(covk, i + 1, k, transposed)) //covk.getParameterValue(k, i + 1))
                        * getCovariateValue(covk, i, k, transposed) //covk.getParameterValue(k, i)
                ;

            }

            hessian[k] = h;
        }

        return hessian;
    }

    private double[] getGradientLogDensity() {

        checkIntervals();

        final int dim = popSizeParameter.getSize();
        double[] gradLogDens = new double[dim];
        double[] gamma = getMeanAdjustedGamma();

        double currentPrec = precisionParameter.getParameterValue(0);

        gradLogDens[0] = -currentPrec * (gamma[0] - gamma[1])
                - numCoalEvents[0] + sufficientStatistics[0]
                * Math.exp(-popSizeParameter.getParameterValue(0));

        for (int i = 1; i < (dim - 1); i++) {
            gradLogDens[i] = -currentPrec * (-gamma[i - 1] + 2 * gamma[i] - gamma[i + 1])
                    - numCoalEvents[i] + sufficientStatistics[i]
                    * Math.exp(-popSizeParameter.getParameterValue(i));
        }

        gradLogDens[dim - 1] = -currentPrec * (gamma[dim - 1] - gamma[dim - 2])
                - numCoalEvents[dim - 1] + sufficientStatistics[dim - 1]
                * Math.exp(-popSizeParameter.getParameterValue(dim - 1));

        return gradLogDens;

//        // TODO Am unclear why the code below does not work for Hessian (???)
//        DenseVector currentGamma = new DenseVector(popSizeParameter.getParameterValues());
//
//        // Use same code as skygridHelper -- TODO remove code duplication
//        skygridHelper.updateGammaWithCovariates(currentGamma);
//        SymmTridiagMatrix currentQ = getScaledWeightMatrix(
//                precisionParameter.getParameterValue(0),
//                lambdaParameter.getParameterValue(0));
//
//        DenseVector g = new DenseVector(fieldLength);
//        currentQ.mult(currentGamma, g);
//
//        double[] gradient = g.getData();
//
//        for (int i = 0; i < fieldLength; ++i) {
//            gradient[i] += -numCoalEvents[i] + sufficientStatistics[i] * Math.exp(-popSizeParameter.getParameterValue(i));
//        }
//
//        return gradient;
    }

    private double[] getDiagonalHessianLogDensity() {

        checkIntervals();

        double[] hessianLogDens = new double[popSizeParameter.getSize()];
        double[] currentGamma = popSizeParameter.getParameterValues();
        double currentPrec = precisionParameter.getParameterValue(0);
        int popSizeDim = popSizeParameter.getSize();

        hessianLogDens[0] = -currentPrec
                - sufficientStatistics[0] * Math.exp(-currentGamma[0]);

        hessianLogDens[popSizeDim - 1] = -currentPrec
                - sufficientStatistics[popSizeDim - 1] * Math.exp(-currentGamma[popSizeDim - 1]);

        for (int i = 1; i < (popSizeDim - 1); i++) {
            hessianLogDens[i] = -2 * currentPrec - sufficientStatistics[i] * Math.exp(-currentGamma[i]);
        }

        return hessianLogDens;
    }

    public double[] getGridPoints() {
        return gridPoints.clone();
    }

    private void checkIntervals() {
        if (!intervalsKnown) {
            //intervalsKnown -> false when handleModelChanged event occurs in super.
            wrapSetupIntervals();
            setupSufficientStatistics();
            intervalsKnown = true;
        }
    }

    class SkygridHelper {

        SkygridHelper() {
        }

        void updateGammaWithCovariates(DenseVector currentGamma) {
            // Do nothing
        }

        private DenseVector getMeanAdjustedGamma() {
            DenseVector currentGamma = new DenseVector(popSizeParameter.getParameterValues());
            updateGammaWithCovariates(currentGamma);
            return currentGamma;
        }

        double handleMissingValues() {
            return 0.0;
        }

        double getLogFieldLikelihood() {

            checkIntervals(); // TODO Is this really necessary?  Computation below does not appear to depend on intervals.

            DenseVector diagonal1 = new DenseVector(fieldLength);
            DenseVector currentGamma = getMeanAdjustedGamma();

            double currentLike = handleMissingValues();

            SymmTridiagMatrix currentQ = getScaledWeightMatrix(precisionParameter.getParameterValue(0), lambdaParameter.getParameterValue(0));
            currentQ.mult(currentGamma, diagonal1);

            currentLike += 0.5 * (fieldLength - 1) * Math.log(precisionParameter.getParameterValue(0)) - 0.5 * currentGamma.dot(diagonal1);
            if (lambdaParameter.getParameterValue(0) == 1) {
                currentLike -= (fieldLength - 1) / 2.0 * LOG_TWO_TIMES_PI;
            } else {
                currentLike -= fieldLength / 2.0 * LOG_TWO_TIMES_PI;
            }

            return currentLike;
        }
    }

    class SkygridCovariateHelper extends SkygridHelper {

        SkygridCovariateHelper() {
        }

        @Override
        protected void updateGammaWithCovariates(DenseVector currentGamma) {

            assert (beta != null);

            // Handle betaParameter / designMatrix

            if (NEW_APPROACH) {

                final int N = currentGamma.size();
                double[] update = new double[N];

                if (dMatrix != null) {
                    final int K = dMatrix.getColumnDimension();

                    if (N != dMatrix.getRowDimension()) {
                        throw new RuntimeException("Incorrect covariate dimensions (" + N + " != "
                                + dMatrix.getRowDimension() + ")");
                    }

                    for (int i = 0; i < N; ++i) {
                        for (int j = 0; j < K; ++j) {
                            update[i] += dMatrix.getParameterValue(i, j) * betaParameter.getParameterValue(j);
                        }
                    }
                }

                if (covariates != null) {
                    if (beta.size() != covariates.size()) {
                        throw new RuntimeException("beta.size(" + beta.size() + ") != covariates.size(" + covariates.size() + ")");
                    }

                    for (int k = 0; k < beta.size(); ++k) {

                        Parameter b = beta.get(k);
                        Parameter d = delta.get(k);
                        final int J = b.getDimension();
                        MatrixParameter covariate = covariates.get(k);
                        boolean transposed = isTransposed(N, J, covariate);

                        for (int i = 0; i < N; ++i) {
                            for (int j = 0; j < J; ++j) {
                                update[i] += getCovariateValue(covariate, i, j, transposed) * b.getParameterValue(j) *
                                        d.getParameterValue(j);
                            }
                        }
                    }
                }

                for (int i = 0; i < N; ++i) {
                    currentGamma.set(i, currentGamma.get(i) - update[i]);
                }

            } else {
                DenseVector currentBeta = new DenseVector(beta.size());

                for (int i = 0; i < beta.size(); i++) {
                    currentBeta.set(i, beta.get(i).getParameterValue(0) * delta.get(i).getParameterValue(0));
                }

                //int numMissing = fieldLength - lastObservedIndex;
                //DenseVector tempVectCov = new DenseVector(numMissing);

                //System.err.println("covariates.size(): " + covariates.size());
                //System.err.println("covariates.get(0).getColumnDimension: " + covariates.get(0).getColumnDimension());
                //System.err.println("covariates.get(0).getRowDimension: " + covariates.get(0).getRowDimension());

                if (covariates != null) {

                    for (int i = 0; i < covariates.size(); i++) {
                        for (int j = 0; j < covariates.get(i).getColumnDimension(); j++) {
                            // System.err.println("j: " + j);
                            // System.err.println("covariates.get(i).getParameterValue(0,j): " + covariates.get(i).getParameterValue(0,j));
                            currentGamma.set(j, currentGamma.get(j) - covariates.get(i).getParameterValue(0, j) * currentBeta.get(i));
                        }
                    }
                }

                throw new RuntimeException("Should not get here.");
            }
        }
    }

    private double getCovariateValue(MatrixParameter matrix, int i, int j, boolean transposed) {
        if (transposed) {
            return matrix.getParameterValue(j, i);
        } else {
            return matrix.getParameterValue(i, j);
        }
    }

    private boolean isTransposed(int N, int J, MatrixParameter matrix) {
        if (J == matrix.getRowDimension() && N == matrix.getColumnDimension()) {
            return true;
        } else if (J == matrix.getColumnDimension() && N == matrix.getRowDimension()) {
            return false;
        } else {
            throw new RuntimeException("Incorrect dimensions in " + matrix.getId() + " (r=" + matrix.getRowDimension() +
                    ",c=" + matrix.getColumnDimension() + ")");
        }
    }

    private static final boolean NEW_APPROACH = true;

    class SkygridMissingCovariateHelper extends SkygridCovariateHelper {

        SkygridMissingCovariateHelper() {
        }

        @Override
        protected double handleMissingValues() {

            assert (covPrecParametersRecent != null);
            assert (covariates != null);
            assert (covPrecParametersDistant != null);

            int numMissing;
            DenseVector tempVectMissingCov;
            SymmTridiagMatrix missingCovQ;
            DenseVector tempVectMissingCov2;
            int numMissingRecent;

            double currentLike = 0.0;

            if (lastObservedIndex != null) {
                for (int i = 0; i < covPrecParametersDistant.size(); i++) {

                    numMissing = fieldLength - lastObservedIndex[i];
                    tempVectMissingCov = new DenseVector(numMissing);
                    tempVectMissingCov2 = new DenseVector(numMissing);

                    missingCovQ = getScaledWeightMatrixForMissingCovDistant(covPrecParametersDistant.get(i).getParameterValue(0), i,
                            lastObservedIndex[i]);

                    for (int j = 0; j < numMissing; j++) {
                        tempVectMissingCov.set(j, covariates.get(distIndices[i] - 1).getParameterValue(0, lastObservedIndex[i] + j) -
                                covariates.get(distIndices[i] - 1).getParameterValue(0, lastObservedIndex[i] - 1));
                    }

                    missingCovQ.mult(tempVectMissingCov, tempVectMissingCov2);
                    currentLike += 0.5 * (numMissing) * Math.log(covPrecParametersDistant.get(i).getParameterValue(0))
                            - 0.5 * tempVectMissingCov.dot(tempVectMissingCov2);
                }
            }

            if (firstObservedIndex != null) {

                for (int i = 0; i < covPrecParametersRecent.size(); i++) {

                    numMissingRecent = firstObservedIndex[i] - 1;
                    tempVectMissingCov = new DenseVector(numMissingRecent);
                    tempVectMissingCov2 = new DenseVector(numMissingRecent);

                    missingCovQ = getScaledWeightMatrixForMissingCovRecent(covPrecParametersRecent.get(i).getParameterValue(0), i,
                            firstObservedIndex[i]);

                    for (int j = 0; j < numMissingRecent; j++) {
                        tempVectMissingCov.set(j, covariates.get(recIndices[i] - 1).getParameterValue(0, j) -
                                covariates.get(recIndices[i] - 1).getParameterValue(0, firstObservedIndex[i] - 1));
                    }

                    missingCovQ.mult(tempVectMissingCov, tempVectMissingCov2);
                    currentLike += 0.5 * (numMissingRecent) * Math.log(covPrecParametersRecent.get(i).getParameterValue(0))
                            - 0.5 * tempVectMissingCov.dot(tempVectMissingCov2);
                }
            }
            return currentLike;
        }

    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "Skyride coalescent";
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
