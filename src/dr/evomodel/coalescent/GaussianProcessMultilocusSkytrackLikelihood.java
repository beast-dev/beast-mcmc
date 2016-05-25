/*
 * GaussianProcessMultilocusSkytrackLikelihood.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

//import dr.evolution.coalescent.IntervalType;
import dr.evolution.coalescent.TreeIntervals;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
//import dr.evomodelxml.coalescent.GMRFSkyrideLikelihoodParser;
import dr.evomodelxml.coalescent.GaussianProcessSkytrackLikelihoodParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
//import no.uib.cipr.matrix.DenseVector;
//import no.uib.cipr.matrix.SymmTridiagMatrix;

import java.util.ArrayList;
//import java.util.Arrays;
import java.util.List;

/**
 * @author Julia Palacios
 * @author Marc A. Suchard
 * @author Vladimir Minin
 * @author Michael Karcher
 */

public class GaussianProcessMultilocusSkytrackLikelihood extends GaussianProcessSkytrackLikelihood implements MultiLociTreeSet {


//    private double cutOff;
//    private int numGridPoints;
//    protected int oldFieldLength;
    // number of coalescent events which occur in an interval with constant population size
//    protected double[] numCoalEvents;
//    protected double[] storedNumCoalEvents;
//    protected double[] gridPoints;
//    protected double theLastTime;
//    protected double diagonalValue;
    // sortedPoints[i][0] is the time of the i-th grid point or sampling or coalescent event
    // sortedPoints[i][1] is 0 if the i-th point is a grid point, 1 if it's a sampling point, and 2 if it's a coalescent point
    // sortedPoints[i][2] is the number of lineages present in the interval starting at time sortedPoints[i][0]

//    protected Parameter phiParameter;
//    protected SymmTridiagMatrix precMatrix;
//    protected SymmTridiagMatrix storedPrecMatrix;

    public GaussianProcessMultilocusSkytrackLikelihood(List<Tree> treeList,
                                                       Parameter precParameter,
                                                       boolean rescalebyRootHeight,
                                                       Parameter lambda_bound,
                                                       Parameter lambda_parameter,
                                                       Parameter popParameter,
                                                       Parameter alpha_parameter,
                                                       Parameter beta_parameter,
                                                       Parameter change_points,
                                                       Parameter GPtype,
                                                       Parameter GPcounts,
                                                       Parameter coalfactor,
                                                       Parameter CoalCounts,
                                                       Parameter numPoints,
                                                       Parameter Tmrca) {

        super(GaussianProcessSkytrackLikelihoodParser.SKYTRACK_LIKELIHOOD);


        this.popSizeParameter = popParameter;
        this.Tmrca = Tmrca;
        this.changePoints=change_points;
        this.numPoints=numPoints;
        this.precisionParameter = precParameter;
        this.lambdaParameter = lambda_parameter;
        this.betaParameter = beta_parameter;
        this.alphaParameter=alpha_parameter;
        this.rescaleByRootHeight=rescalebyRootHeight;
        this.lambda_boundParameter=lambda_bound;
        this.GPcounts=GPcounts;
        this.GPtype=GPtype;
        this.coalfactor=coalfactor;
        this.CoalCounts=CoalCounts;

//        System.out.println("numGridPoints: " + numGridPoints);
//        setupGridPoints();

        addVariable(popSizeParameter);
        addVariable(precisionParameter);
        addVariable(changePoints);
        addVariable(numPoints);
        addVariable(GPcounts);
        addVariable(GPtype);
        addVariable(coalfactor);
        addVariable(lambda_boundParameter);
        addVariable(CoalCounts);


        setTree(treeList);

        numintervals=getIntervalCount(); //Not sure, check. I want the number of intervals defined by either
//        sampling times or coalescent times of all trees
        System.err.println("by getIntervalCount, I get:"+getIntervalCount());

        numcoalpoints = getCorrectFieldLength();
        System.err.println("by getCorrectField"+getCorrectFieldLength());

        GPcoalfactor = new double[numintervals];
        backupIntervals=new double[numintervals];
        GPCoalInterval=new double[numcoalpoints];
        storedGPCoalInterval=new double[numcoalpoints];
        CoalPosIndicator= new int[numcoalpoints];
        storedCoalPosIndicator=new int[numcoalpoints];
        CoalTime=new double[numcoalpoints];
        storedCoalTime=new double[numcoalpoints];
        storedGPcoalfactor = new double[numintervals];
        GPcounts.setDimension(numintervals);
        CoalCounts.setDimension(numcoalpoints);
        GPtype.setDimension(numcoalpoints);
        numPoints.setParameterValue(0,numcoalpoints);
        popSizeParameter.setDimension(numcoalpoints);
        changePoints.setDimension(numcoalpoints);
        coalfactor.setDimension(numcoalpoints);

        initializationReport();

//        wrapSetupIntervals();

//        coalescentIntervals = new double[oldFieldLength];
//        storedCoalescentIntervals = new double[oldFieldLength];
//        sufficientStatistics = new double[fieldLength];
//        storedSufficientStatistics = new double[fieldLength];
//        numCoalEvents = new double[fieldLength];
//        storedNumCoalEvents = new double[fieldLength];

//        setupGMRFWeights();
        setupSufficientStatistics();

//        addStatistic(new DeltaStatistic());



        /* Force all entries in groupSizeParameter = 1 for compatibility with Tracer */

    }



//                   OK
    protected void setTree(List<Tree> treeList) {
        treesSet = this;
        this.treeList = treeList;
        makeTreeIntervalList(treeList, true);
        numTrees = treeList.size();
    }
//                  OK
    private void makeTreeIntervalList(List<Tree> treeList, boolean add) {
        if (intervalsList == null) {
            intervalsList = new ArrayList<TreeIntervals>();
        } else {
            intervalsList.clear();
        }
        for (Tree tree : treeList) {
            numIntervals+=tree.getExternalNodeCount();
            intervalsList.add(new TreeIntervals(tree));
            if (add && tree instanceof TreeModel) {
                addModel((TreeModel) tree);
            }
        }
    }

    protected int getCorrectFieldLength() {

        return numIntervals-treeList.size();
//        TODO add correction when not all samples are gathered at the same time
    }
//
//    protected int getCorrectOldFieldLength() {
//        int tips = 0;
//        for (Tree tree : treeList) {
//            tips += tree.getExternalNodeCount();
//        }
//        return tips - treeList.size();
//    }

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
                throw new RuntimeException("Unknown tree modified in GMRFMultilocusSkyrideLikelihood");
            }
        } else {
            throw new RuntimeException("Unknown object modified in GMRFMultilocusSkyrideLikelihood");
        }
    }

    public void initializationReport() {
        System.out.println("Creating a GP based estimation of effective population size trajectories for multiple loci:");
        System.out.println("\tIf you publish results using this model, please reference: ");
        System.out.println("\t\tPalacios, Minin and Suchard (XXXX)");
       }
//
//    public void wrapSetupIntervals() {
//        // Do nothing
//    }

    int numTrees;
    int numIntervals;

//
//    protected void setupGridPoints() {
//        if (gridPoints == null) {
//            gridPoints = new double[numGridPoints];
//        } else {
//            Arrays.fill(gridPoints, 0);
//        }
//
//        for (int pt = 0; pt < numGridPoints; pt++) {
//            gridPoints[pt] = (pt + 1) * (cutOff / numGridPoints);
//        }
//    }

    protected void setupSufficientStatistics() {

//        //numCoalEvents = new double[fieldLength];
//        //sufficientStatistics = new double[fieldLength];
//
//        Arrays.fill(numCoalEvents, 0);
//        Arrays.fill(sufficientStatistics, 0);
//
//        //index of smallest grid point greater than at least one sampling/coalescent time in current tree
//        int minGridIndex;
//        //index of greatest grid point less than at least one sampling/coalescent time in current tree
//        int maxGridIndex;
//
//        int numLineages;
//
//        int currentGridIndex;
//        int currentTimeIndex;
//
//        double currentTime;
//        double nextTime;
//
//        //time of last coalescent event in tree
//        double lastCoalescentTime;
//
//        for (int i = 0; i < numTrees; i++) {
//            currentTimeIndex = 0;
//            currentTime = intervalsList.get(i).getIntervalTime(currentTimeIndex);
//            nextTime = intervalsList.get(i).getIntervalTime(currentTimeIndex + 1);
//            while (nextTime <= currentTime) {
//                currentTimeIndex++;
//                currentTime = intervalsList.get(i).getIntervalTime(currentTimeIndex);
//                nextTime = intervalsList.get(i).getIntervalTime(currentTimeIndex + 1);
//            }
//
//
//            numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
//            minGridIndex = 0;
//            while (gridPoints[minGridIndex] <= currentTime) {
//                minGridIndex++;
//            }
//            currentGridIndex = minGridIndex;
//
//            lastCoalescentTime = currentTime + intervalsList.get(i).getTotalDuration();
//
//            theLastTime = lastCoalescentTime;
//
//            maxGridIndex = numGridPoints - 1;
//            while ((maxGridIndex >= 0) && (gridPoints[maxGridIndex] >= lastCoalescentTime)) {
//                maxGridIndex = maxGridIndex - 1;
//            }
//
//            if (maxGridIndex >= 0) {
//
//                //from likelihood of interval between first sampling time and gridPoints[minGridIndex]
//
//                while (nextTime < gridPoints[currentGridIndex]) {
//
//                    //check to see if interval ends with coalescent event
//                    if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
//
//                        numCoalEvents[currentGridIndex]++;
//                    }
//                    sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (nextTime - currentTime) * numLineages * (numLineages - 1) * 0.5;
//                    currentTime = nextTime;
//                    currentTimeIndex++;
//                    nextTime = intervalsList.get(i).getIntervalTime(currentTimeIndex + 1);
//
//                    while (nextTime <= currentTime) {
//                        currentTimeIndex++;
//                        currentTime = intervalsList.get(i).getIntervalTime(currentTimeIndex);
//                        nextTime = intervalsList.get(i).getIntervalTime(currentTimeIndex + 1);
//                    }
//
//                    numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
//
//                }
//
//                sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (gridPoints[currentGridIndex] - currentTime) * numLineages * (numLineages - 1) * 0.5;
//
//                currentGridIndex++;
//
//
//                //from likelihood of intervals between gridPoints[minGridIndex] and gridPoints[maxGridIndex]
//
//                while (currentGridIndex <= maxGridIndex) {
//                    if (nextTime >= gridPoints[currentGridIndex]) {
//                        sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (gridPoints[currentGridIndex] - gridPoints[currentGridIndex - 1]) * numLineages * (numLineages - 1) * 0.5;
//
//                        currentGridIndex++;
//                    } else {
//
//                        sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (nextTime - gridPoints[currentGridIndex - 1]) * numLineages * (numLineages - 1) * 0.5;
//
//                        //check to see if interval ends with coalescent event
//                        if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
//                            numCoalEvents[currentGridIndex]++;
//                        }
//                        currentTime = nextTime;
//                        currentTimeIndex++;
//                        nextTime = intervalsList.get(i).getIntervalTime(currentTimeIndex + 1);
//                        while (nextTime <= currentTime) {
//                            currentTimeIndex++;
//                            currentTime = intervalsList.get(i).getIntervalTime(currentTimeIndex);
//                            nextTime = intervalsList.get(i).getIntervalTime(currentTimeIndex + 1);
//                        }
//
//                        numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
//
//                        while (nextTime < gridPoints[currentGridIndex]) {
//                            //check to see if interval is coalescent interval or sampling interval
//                            if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
//                                numCoalEvents[currentGridIndex]++;
//                            }
//                            sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (nextTime - currentTime) * numLineages * (numLineages - 1) * 0.5;
//
//                            currentTime = nextTime;
//                            currentTimeIndex++;
//                            nextTime = intervalsList.get(i).getIntervalTime(currentTimeIndex + 1);
//                            while (nextTime <= currentTime) {
//                                currentTimeIndex++;
//                                currentTime = intervalsList.get(i).getIntervalTime(currentTimeIndex);
//                                nextTime = intervalsList.get(i).getIntervalTime(currentTimeIndex + 1);
//                            }
//
//                            numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
//
//                        }
//                        sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (gridPoints[currentGridIndex] - currentTime) * numLineages * (numLineages - 1) * 0.5;
//
//                        currentGridIndex++;
//                    }
//                }
//
//                //from likelihood of interval between gridPoints[maxGridIndex] and lastCoalescentTime
//
//                sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (nextTime - gridPoints[currentGridIndex - 1]) * numLineages * (numLineages - 1) * 0.5;
//
//                //check to see if interval ends with coalescent event
//                if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
//                    numCoalEvents[currentGridIndex]++;
//                }
//
//                currentTime = nextTime;
//                currentTimeIndex++;
//
//                while ((currentTimeIndex + 1) < intervalsList.get(i).getIntervalCount()) {
//                    // currentTime = nextTime;
//                    // currentTimeIndex++;
//
//                    nextTime = intervalsList.get(i).getIntervalTime(currentTimeIndex + 1);
//                    while (nextTime <= currentTime) {
//                        currentTimeIndex++;
//                        currentTime = intervalsList.get(i).getIntervalTime(currentTimeIndex);
//                        nextTime = intervalsList.get(i).getIntervalTime(currentTimeIndex + 1);
//                    }
//
//                    numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
//
//                    //check to see if interval is coalescent interval or sampling interval
//
//
//                    if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
//                        numCoalEvents[currentGridIndex]++;
//                    }
//                    sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (nextTime - currentTime) * numLineages * (numLineages - 1) * 0.5;
//                    currentTime = nextTime;
//                    currentTimeIndex++;
//
//                }
//            } else {
//                while ((currentTimeIndex + 1) < intervalsList.get(i).getIntervalCount()) {
//                    //check to see if interval is coalescent interval or sampling interval
//                    if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
//                        numCoalEvents[currentGridIndex]++;
//                    }
//                    sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (nextTime - currentTime) * numLineages * (numLineages - 1) * 0.5;
//
//                    currentTime = nextTime;
//                    currentTimeIndex++;
//                    if ((currentTimeIndex + 1) < intervalsList.get(i).getIntervalCount()) {
//                        nextTime = intervalsList.get(i).getIntervalTime(currentTimeIndex + 1);
//
//                        while (nextTime <= currentTime) {
//                            currentTimeIndex++;
//                            currentTime = intervalsList.get(i).getIntervalTime(currentTimeIndex);
//                            nextTime = intervalsList.get(i).getIntervalTime(currentTimeIndex + 1);
//                        }
//
//                        numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
//
//                    }
//
//                }
//
//
//            }
//        }

    }
//
//    public double[] getNumCoalEvents() {
//        return numCoalEvents;
//    }

    protected double calculateLogCoalescentLikelihood() {

//        if (!intervalsKnown) {
//            // intervalsKnown -> false when handleModelChanged event occurs in super.
//            wrapSetupIntervals();
//            setupSufficientStatistics();
//            intervalsKnown = true;
//        }

        // Matrix operations taken from block update sampler to calculate data likelihood and field prior

        double currentLike = 0;
//        double[] currentGamma = popSizeParameter.getParameterValues();
//
//        for (int i = 0; i < fieldLength; i++) {
//            currentLike += -numCoalEvents[i] * currentGamma[i] - sufficientStatistics[i] * Math.exp(-currentGamma[i]);
//        }

        return currentLike;
    }


    protected double calculateLogFieldLikelihood() {

//        if (!intervalsKnown) {
//            //intervalsKnown -> false when handleModelChanged event occurs in super.
//            wrapSetupIntervals();
//            setupSufficientStatistics();
//            intervalsKnown = true;
//        }

        double currentLike = 0;
//        DenseVector diagonal1 = new DenseVector(fieldLength);
//        DenseVector currentGamma = new DenseVector(popSizeParameter.getParameterValues());
//
//        SymmTridiagMatrix currentQ = getScaledWeightMatrix(precisionParameter.getParameterValue(0), lambdaParameter.getParameterValue(0));
//        currentQ.mult(currentGamma, diagonal1);
//
//        //        currentLike += 0.5 * logGeneralizedDeterminant(currentQ) - 0.5 * currentGamma.dot(diagonal1);
//
//        currentLike += 0.5 * (fieldLength - 1) * Math.log(precisionParameter.getParameterValue(0)) - 0.5 * currentGamma.dot(diagonal1);
//        if (lambdaParameter.getParameterValue(0) == 1) {
//            currentLike -= (fieldLength - 1) / 2.0 * LOG_TWO_TIMES_PI;
//        } else {
//            currentLike -= fieldLength / 2.0 * LOG_TWO_TIMES_PI;
//        }

        return currentLike;
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogCoalescentLikelihood();
//            logFieldLikelihood = calculateLogFieldLikelihood();
            likelihoodKnown = true;
        }

        return logLikelihood;
//        return logLikelihood + logFieldLikelihood;
    }

//    protected void setupGMRFWeights() {
//
//        //setupSufficientStatistics();
//
//        //Set up the weight Matrix
//        double[] offdiag = new double[fieldLength - 1];
//        double[] diag = new double[fieldLength];
//
//        diagonalValue = 2;
//        //First set up the offdiagonal entries;
//
//        for (int i = 0; i < fieldLength - 1; i++) {
//            offdiag[i] = -1;
//        }
//
//        //Then set up the diagonal entries;
//        for (int i = 1; i < fieldLength - 1; i++) {
//            //	diag[i] = -(offdiag[i] + offdiag[i - 1]);
//            diag[i] = diagonalValue;
//        }
//        //Take care of the endpoints
//        //diag[0] = -offdiag[0];
//        //diag[fieldLength - 1] = -offdiag[fieldLength - 2];
//        diag[0] = diagonalValue - 1.0;
//        diag[fieldLength - 1] = diagonalValue - 1.0;
//
//
//        weightMatrix = new SymmTridiagMatrix(diag, offdiag);
//
//    }

    protected double getFieldScalar() {
        return 1.0;
    }


    private List<Tree> treeList;
    private List<TreeIntervals> intervalsList;

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
        return 1.0;
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

//    protected void storeState() {
//        // System.arraycopy(numCoalEvents, 0, storedNumCoalEvents, 0, numCoalEvents.length);
//        super.storeState();
//        System.arraycopy(numCoalEvents, 0, storedNumCoalEvents, 0, numCoalEvents.length);
//        // storedPrecMatrix = precMatrix.copy();
//    }


//    protected void restoreState() {
//        super.restoreState();
//
//        // Swap pointers
//        double[] tmp = numCoalEvents;
//        numCoalEvents = storedNumCoalEvents;
//        storedNumCoalEvents = tmp;
//    }

//    public int getCoalescentIntervalLineageCount(int i) {
//        return 0;  //To change body of implemented methods use File | Settings | File Templates.
//    }
//
//    public IntervalType getCoalescentIntervalType(int i) {
//        return null;  //To change body of implemented methods use File | Settings | File Templates.
//    }
}

