/*
 * ContinuousDiffusionStatistic.java
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

package dr.evomodel.continuous;

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.treelikelihood.MarkovJumpsBeagleTreeLikelihood;
import dr.app.util.Arguments;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeStatistic;
import dr.geo.KMLCoordinates;
import dr.geo.Polygon2D;
import dr.geo.contouring.ContourMaker;
import dr.geo.contouring.ContourPath;
import dr.geo.contouring.ContourWithSynder;
import dr.geo.math.SphericalPolarCoordinates;
import dr.inference.model.Statistic;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.stats.DiscreteStatistics;
import dr.stats.Regression;
import dr.xml.*;
import org.jdom.Element;
import org.apache.commons.math.stat.ranking.NaturalRanking;

import java.util.*;

/**
 * @author Marc Suchard
 * @author Philippe Lemey
 * @author Andrew Rambaut
 */
public class ContinuousDiffusionStatistic extends Statistic.Abstract {

    public static final String CONTINUOUS_DIFFUSION_STATISTIC = "continuousDiffusionStatistic";
    public static final String DIFFUSION_RATE_STATISTIC = "diffusionRateStatistic";
    public static final String TREE_DISPERSION_STATISTIC = "treeDispersionStatistic";

    public static final String USE_GREATCIRCLEDISTANCES = "greatCircleDistance";
    public static final String MODE = "mode";
    public static final String MEDIAN = "median";
    public static final String AVERAGE = "average";  // average over all branches
    public static final String WEIGHTED_AVERAGE = "weightedAverage"; // weighted average (=total distance/total time)
    public static final String COEFFICIENT_OF_VARIATION = "coefficientOfVariation"; // weighted average (=total distance/total time)
    public static final String SPEARMAN = "spearman";
    public static final String CORRELATION_COEFFICIENT = "correlationCoefficient";
    public static final String DISTANCE_TIME_CORRELATION = "distanceTimeCorrelation";
    public static final String R_SQUARED = "Rsquared";
    public static final String STATISTIC = "statistic";
    public static final String TRAIT = "trait";
    public static final String TRAIT2DAREA = "trait2Darea";
    public static final String DIMENSION = "dimension";
    public static final String DIFFUSION_TIME = "diffusionTime";
    public static final String DIFFUSION_DISTANCE = "diffusionDistance";
    public static final String DIFFUSION_RATE = "diffusionRate"; // weighted average (=total distance/total time)
    public static final String WAVEFRONT_DISTANCE = "wavefrontDistance"; // weighted average (=total distance/total time)
    public static final String WAVEFRONT_DISTANCE_PHYLO = "wavefrontDistancePhylo"; // weighted average (=total brnach distance/total time)
    public static final String WAVEFRONT_RATE = "wavefrontRate"; // weighted average (=total distance/total time)
    public static final String DIFFUSION_COEFFICIENT = "diffusionCoefficient";
    public static final String HEIGHT_UPPER = "heightUpper";
    public static final String HEIGHT_LOWER = "heightLower";
    public static final String HEIGHT_LOWER_SERIE = "heightLowerSerie";
    public static final String CUMULATIVE = "cumulative";
    public static final String DISCRETE_STATE = "discreteState";
    public static final Integer SITE = 0;
    public static final Integer NUMBER_OF_HISTORY_ENTRIES = 3;
    public static final String NOISE = "noise";
    public static final String TAXA = "taxa";
    public static final String BRANCHSET = "branchSet";
    public static final String ALL = "all";
    public static final String CLADE = "clade";
    public static final String BACKBONE = "backbone";
    public static final String BACKBONE_TIME = "backboneTime";

    public ContinuousDiffusionStatistic(String name, List<AbstractMultivariateTraitLikelihood> traitLikelihoods,
                                        boolean greatCircleDistances, Mode mode,
                                        summaryStatistic statistic, double heightUpper, double heightLower,
                                        double[] lowerHeights, boolean cumulative, boolean trueNoise, int dimension,
                                        TaxonList taxonList, BranchSet branchset, Double backboneTime,
                                        String stateString, MarkovJumpsBeagleTreeLikelihood markovJumpLikelihood) {
        super(name);
        this.traitLikelihoods = traitLikelihoods;
        this.useGreatCircleDistances = greatCircleDistances;
        summaryMode =  mode;
        summaryStat = statistic;
        this.heightUpper = heightUpper;

        if (lowerHeights == null){
            heightLowers =  new double[]{heightLower};
        } else {
            heightLowers = extractUnique(lowerHeights);
            Arrays.sort(heightLowers);
            reverse(heightLowers);
        }
        this.cumulative = cumulative;
        this.trueNoise = trueNoise;

        this.dimension = dimension;

        this.taxonList = taxonList;
        this.branchset = branchset;
        this.backboneTime = backboneTime;

        this.stateString =  stateString;
//        this.stateInt = stateInt;
        this.markovJumpLikelihood = markovJumpLikelihood;
    }

    public int getDimension() {
        return heightLowers.length;
    }

    public double getStatisticValue(int dim) {

        double treeLength = 0;
        double treeDistance = 0;
        double totalMaxDistanceFromRoot = 0;
        double maxDistanceFromRootCumulative = 0; // can only be used when cumulative and not associated with discrete state (not based on the distances on the branches from the root up that point)
        double maxBranchDistanceFromRoot = 0;
        double maxDistanceOverTimeFromRootWA = 0;  // can only be used when cumulative and not associated with discrete state (not based on the distances on the branches from the root up that point)
        double maxBranchDistanceOverTimeFromRootWA = 0;
        List<Double> rates = new ArrayList<Double>();
        List<Double> distances = new ArrayList<Double>();
        List<Double> times = new ArrayList<Double>();
        List<Double> traits = new ArrayList<Double>();
        List<double[]> traits2D = new ArrayList<double[]>();
        //double[] diffusionCoefficients =  null;
        List<Double> diffusionCoefficients = new ArrayList<Double>();
        double waDiffusionCoefficient =  0;

        double lowerHeight = heightLowers[dim];
        double upperHeight = Double.MAX_VALUE;
        if (heightLowers.length == 1){
            upperHeight = heightUpper;
        } else {
            if (dim > 0) {
                if (!cumulative) {
                    upperHeight = heightLowers[dim -1];
                }
            }
        }

        for (AbstractMultivariateTraitLikelihood traitLikelihood : traitLikelihoods) {
            MutableTreeModel tree = traitLikelihood.getTreeModel();
            BranchRateModel branchRates = traitLikelihood.getBranchRateModel();

            String traitName = traitLikelihood.getTraitName();

            for (int i = 0; i < tree.getNodeCount(); i++) {
                NodeRef node = tree.getNode(i);

                if (node != tree.getRoot()) {

                    NodeRef parentNode = tree.getParent(node);

                    boolean testNode = true;
                    if  (branchset.equals(BranchSet.CLADE)){
                        try{
                            testNode = inClade(tree, node, taxonList);
                        } catch (TreeUtils.MissingTaxonException mte) {
                            throw new RuntimeException(mte.toString());
                        }
                    }  else if (branchset.equals(BranchSet.BACKBONE)) {
                        if (backboneTime > 0) {
                            testNode = onAncestralPathTime(tree, node, backboneTime);
                        } else {
                            try{
                                testNode = onAncestralPathTaxa(tree, node, taxonList);
                            } catch (TreeUtils.MissingTaxonException mte) {
                                throw new RuntimeException(mte.toString());
                            }
                        }
                    }

                    if (testNode){

                        if ((tree.getNodeHeight(parentNode) > lowerHeight) && (tree.getNodeHeight(node) < upperHeight)) {

                            double[] trait = traitLikelihood.getTraitForNode(tree, node, traitName);
                            double[] parentTrait = traitLikelihood.getTraitForNode(tree, parentNode, traitName);

                            double[] traitUp = parentTrait;
                            double[] traitLow = trait;

                            double timeUp = tree.getNodeHeight(parentNode);
                            double timeLow = tree.getNodeHeight(node);

                            double rate = (branchRates != null ? branchRates.getBranchRate(tree, node) : 1.0);
//                        System.out.println(rate);
                            MultivariateDiffusionModel diffModel = traitLikelihood.diffusionModel;
                            double[] precision = diffModel.getPrecisionParameter().getParameterValues();

                            History history = null;
                            if (stateString != null) {
                                history = setUpHistory(markovJumpLikelihood.getHistoryForNode(tree, node, SITE), markovJumpLikelihood.getStatesForNode(tree, node)[SITE], markovJumpLikelihood.getStatesForNode(tree, parentNode)[SITE], timeLow, timeUp);
                            }

                            if (tree.getNodeHeight(parentNode) > upperHeight) {
                                timeUp = upperHeight;
                                traitUp = imputeValue(trait, parentTrait, upperHeight, tree.getNodeHeight(node), tree.getNodeHeight(parentNode), precision, rate, trueNoise);
                                if (stateString != null) {
                                    history.truncateUpper(timeUp);
                                }
                            }

                            if (tree.getNodeHeight(node) < lowerHeight) {
                                timeLow = lowerHeight;
                                traitLow = imputeValue(trait, parentTrait, lowerHeight, tree.getNodeHeight(node), tree.getNodeHeight(parentNode), precision, rate, trueNoise);
                                if (stateString != null) {
                                    history.truncateLower(timeLow);
                                }
                            }

                            if (dimension > traitLow.length) {
                                System.err.println("specified trait dimension for continuous trait summary, " + dimension + ", is > dimensionality of trait, " + traitLow.length + ". No trait summarized.");
                            } else {
                                traits.add(traitLow[(dimension - 1)]);
                            }

                            if (traitLow.length == 2){
                                traits2D.add(traitLow);
                            }

                            double time;
                            if (stateString != null) {
                                time = history.getStateTime(stateString);
//                            System.out.println("tine before = "+(timeUp - timeLow)+", time after= "+time);
                            } else {
                                time = timeUp - timeLow;
                            }
                            treeLength += time;
                            times.add(time);

                            //setting up continuous trait values for heights in discrete trait history
                            if (stateString != null) {
                                history.setTraitsforHeights(traitUp, traitLow, precision, rate, trueNoise);
                            }

                            double[] rootTrait = traitLikelihood.getTraitForNode(tree, tree.getRoot(), traitName);
                            double timeFromRoot = (tree.getNodeHeight(tree.getRoot()) - timeLow);

                            if (useGreatCircleDistances && (trait.length == 2)) { // Great Circle distance
                                double distance;
                                if (stateString != null) {
                                    distance = history.getStateGreatCircleDistance(stateString);
                                } else {
                                    distance = getGreatCircleDistance(traitLow, traitUp);
                                }
                                distances.add(distance);

                                if (time > 0) {
                                    treeDistance += distance;
                                    double dc = Math.pow(distance, 2) / (4 * time);
                                    diffusionCoefficients.add(dc);
                                    waDiffusionCoefficient += (dc * time);
                                    rates.add(distance / time);
                                }

                                SphericalPolarCoordinates rootCoord = new SphericalPolarCoordinates(rootTrait[0], rootTrait[1]);
                                double tempDistanceFromRootLow = rootCoord.distance(new SphericalPolarCoordinates(traitUp[0], traitUp[1]));


                                if (tempDistanceFromRootLow > totalMaxDistanceFromRoot) {
                                    totalMaxDistanceFromRoot = tempDistanceFromRootLow;
                                    if (stateString != null) {
                                        double[] stateTimeDistance = getStateTimeAndDistanceFromRoot(tree, node, timeLow, traitLikelihood, traitName, traitLow, precision, branchRates, true);
                                        if (stateTimeDistance[0] > 0) {
                                            maxDistanceFromRootCumulative = tempDistanceFromRootLow * (stateTimeDistance[0] / timeFromRoot);
                                            maxDistanceOverTimeFromRootWA = maxDistanceFromRootCumulative / stateTimeDistance[0];
                                            maxBranchDistanceFromRoot = stateTimeDistance[1];
                                            maxBranchDistanceOverTimeFromRootWA = stateTimeDistance[1] / stateTimeDistance[0];
                                        }
                                    } else {
                                        maxDistanceFromRootCumulative = tempDistanceFromRootLow;
                                        maxDistanceOverTimeFromRootWA = tempDistanceFromRootLow / timeFromRoot;
                                        double[] timeDistance = getTimeAndDistanceFromRoot(tree, node, timeLow, traitLikelihood, traitName, traitLow, true);
                                        maxBranchDistanceFromRoot = timeDistance[1];
                                        maxBranchDistanceOverTimeFromRootWA = timeDistance[1] / timeDistance[0];

                                    }
                                    //distance between traitLow and traitUp for maxDistanceFromRootCumulative
                                    if (timeUp == upperHeight) {
                                        if (time > 0) {
                                            maxDistanceFromRootCumulative = distance;
                                            maxDistanceOverTimeFromRootWA = distance / time;
                                            maxBranchDistanceFromRoot = distance;
                                            maxBranchDistanceOverTimeFromRootWA = distance / time;
                                        }
                                    }
                                }

                            } else {
                                double distance;
                                if (stateString != null) {
                                    distance = history.getStateNativeDistance(stateString);
                                } else {
                                    distance = getNativeDistance(traitLow, traitUp);
                                }
                                distances.add(distance);

                                if (time > 0) {
                                    treeDistance += distance;
                                    double dc = Math.pow(distance, 2) / (4 * time);
                                    diffusionCoefficients.add(dc);
                                    waDiffusionCoefficient += dc * time;
                                    rates.add(distance / time);
                                }

                                double tempDistanceFromRoot = getNativeDistance(traitLow, rootTrait);
                                if (tempDistanceFromRoot > totalMaxDistanceFromRoot) {
                                    totalMaxDistanceFromRoot = tempDistanceFromRoot;
                                    if (stateString != null) {
                                        double[] stateTimeDistance = getStateTimeAndDistanceFromRoot(tree, node, timeLow, traitLikelihood, traitName, traitLow, precision, branchRates, false);
                                        if (stateTimeDistance[0] > 0) {
                                            maxDistanceFromRootCumulative = tempDistanceFromRoot * (stateTimeDistance[0] / timeFromRoot);
                                            maxDistanceOverTimeFromRootWA = maxDistanceFromRootCumulative / stateTimeDistance[0];
                                            maxBranchDistanceFromRoot = stateTimeDistance[1];
                                            maxBranchDistanceOverTimeFromRootWA = stateTimeDistance[1] / stateTimeDistance[0];
                                        }
                                    } else {
                                        maxDistanceFromRootCumulative = tempDistanceFromRoot;
                                        maxDistanceOverTimeFromRootWA = tempDistanceFromRoot / timeFromRoot;
                                        double[] timeDistance = getTimeAndDistanceFromRoot(tree, node, timeLow, traitLikelihood, traitName, traitLow, false);
                                        maxBranchDistanceFromRoot = timeDistance[1];
                                        maxBranchDistanceOverTimeFromRootWA = timeDistance[1] / timeDistance[0];
                                    }
                                    //distance between traitLow and traitUp for maxDistanceFromRootCumulative
                                    if (timeUp == upperHeight) {
                                        if (time > 0) {
                                            maxDistanceFromRootCumulative = distance;
                                            maxDistanceOverTimeFromRootWA = distance / time;
                                            maxBranchDistanceFromRoot = distance;
                                            maxBranchDistanceOverTimeFromRootWA = distance / time;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        if (summaryStat == summaryStatistic.DIFFUSION_RATE){
            if (summaryMode == Mode.AVERAGE) {
                return DiscreteStatistics.mean(toArray(rates));
            } else if (summaryMode == Mode.MEDIAN) {
                return DiscreteStatistics.median(toArray(rates));
            } else if (summaryMode == Mode.COEFFICIENT_OF_VARIATION) {
                final double mean = DiscreteStatistics.mean(toArray(rates));
                return Math.sqrt(DiscreteStatistics.variance(toArray(rates), mean)) / mean;
                //weighted average
            } else {
                return treeDistance / treeLength;
            }
        } else if (summaryStat == summaryStatistic.TRAIT) {
            if (summaryMode == Mode.MEDIAN) {
                return DiscreteStatistics.median(toArray(traits));
            } else if (summaryMode == Mode.COEFFICIENT_OF_VARIATION) {
                // don't compute mean twice
                final double mean = DiscreteStatistics.mean(toArray(traits));
                return Math.sqrt(DiscreteStatistics.variance(toArray(traits), mean)) / mean;
                // default is average. A warning is thrown by the parser when trying to use WEIGHTED_AVERAGE
            } else {
                return DiscreteStatistics.mean(toArray(traits));
            }
        } else if (summaryStat == summaryStatistic.TRAIT2DAREA) {
            double area = getAreaFrom2Dtraits(traits2D, 0.99);
            return  area;
        }  else if (summaryStat == summaryStatistic.DIFFUSION_COEFFICIENT) {
            if (summaryMode == Mode.AVERAGE) {
                return DiscreteStatistics.mean(toArray(diffusionCoefficients));
            } else if (summaryMode == Mode.MEDIAN) {
                return DiscreteStatistics.median(toArray(diffusionCoefficients));
            } else if (summaryMode == Mode.COEFFICIENT_OF_VARIATION) {
                // don't compute mean twice
                final double mean = DiscreteStatistics.mean(toArray(diffusionCoefficients));
                return Math.sqrt(DiscreteStatistics.variance(toArray(diffusionCoefficients), mean)) / mean;
            } else {
                return waDiffusionCoefficient / treeLength;
            }
            //wavefront distance
            //TODO: restrict to non state-specific wavefrontDistance/rate
        }  else if (summaryStat == summaryStatistic.WAVEFRONT_DISTANCE) {
            return maxDistanceFromRootCumulative;
//            return maxBranchDistanceFromRoot;
        } else if (summaryStat == summaryStatistic.WAVEFRONT_DISTANCE_PHYLO) {
            return maxBranchDistanceFromRoot;
            //wavefront rate, only weighted average TODO: extend for average, median, COEFFICIENT_OF_VARIATION?
        }  else if (summaryStat == summaryStatistic.WAVEFRONT_RATE)  {
            return maxDistanceOverTimeFromRootWA;
//            return maxBranchDistanceOverTimeFromRootWA;
        }  else if (summaryStat == summaryStatistic.DIFFUSION_DISTANCE)  {
            return treeDistance;
            //DIFFUSION_TIME
        }  else if (summaryStat == summaryStatistic.DISTANCE_TIME_CORRELATION)  {
            if (summaryMode == Mode.SPEARMAN) {
                return getSpearmanRho(convertDoubles(times),convertDoubles(distances));
            } else if (summaryMode == Mode.R_SQUARED) {
                Regression r = new Regression(convertDoubles(times), convertDoubles(distances));
                return r.getRSquared();
            } else {
                Regression r = new Regression(convertDoubles(times),convertDoubles(distances));
                return r.getCorrelationCoefficient();
            }
         }  else {
            return treeLength;
        }
    }

//    private double getNativeDistance(double[] location1, double[] location2) {
//        return Math.sqrt(Math.pow((location2[0] - location1[0]), 2.0) + Math.pow((location2[1] - location1[1]), 2.0));
//    }

    private double getNativeDistance(double[] location1, double[] location2) {
        int traitDimension = location1.length;
        double sum = 0;
        for (int i = 0; i < traitDimension; i++) {
            sum += Math.pow((location2[i] - location1[i]),2);
        }
        return Math.sqrt(sum);
    }

    public double getGreatCircleDistance(double[] loc1, double[] loc2){
        SphericalPolarCoordinates coord1 = new SphericalPolarCoordinates(loc1[0], loc1[1]);
        SphericalPolarCoordinates coord2 = new SphericalPolarCoordinates(loc2[0], loc2[1]);
        return coord1.distance(coord2);
    }

    private double[] toArray(List<Double> list) {
        double[] returnArray = new double[list.size()];
        for (int i = 0; i < list.size(); i++) {
            returnArray[i] = Double.valueOf(list.get(i).toString());
        }
        return returnArray;
    }

    private double[] imputeValue(double[] nodeValue, double[] parentValue, double time, double nodeHeight, double parentHeight, double[] precisionArray, double rate, boolean trueNoise) {

        final double scaledTimeChild = (time - nodeHeight) * rate;
        final double scaledTimeParent = (parentHeight - time) * rate;
        final double scaledWeightTotal = 1.0 / scaledTimeChild + 1.0 / scaledTimeParent;
        final int dim = nodeValue.length;

        double[][] precision = new double[dim][dim];
        int counter = 0;
        for (int a = 0; a < dim; a++){
            for (int b = 0; b < dim; b++){
                precision[a][b] = precisionArray[counter];
                counter++ ;
            }
        }

        if (scaledTimeChild == 0)
            return nodeValue;

        if (scaledTimeParent == 0)
            return parentValue;

        // Find mean value, weighted average
        double[] mean = new double[dim];
        double[][] scaledPrecision = new double[dim][dim];

        for (int i = 0; i < dim; i++) {
            mean[i] = (nodeValue[i] / scaledTimeChild + parentValue[i] / scaledTimeParent) / scaledWeightTotal;
            if (trueNoise) {
                for (int j = i; j < dim; j++)
                    scaledPrecision[j][i] = scaledPrecision[i][j] = precision[i][j] * scaledWeightTotal;
            }
        }

//        System.out.print(time+"\t"+nodeHeight+"\t"+parentHeight+"\t"+scaledTimeChild+"\t"+scaledTimeParent+"\t"+scaledWeightTotal+"\t"+mean[0]+"\t"+mean[1]+"\t"+scaledPrecision[0][0]+"\t"+scaledPrecision[0][1]+"\t"+scaledPrecision[1][0]+"\t"+scaledPrecision[1][1]);

        if (trueNoise) {
            mean = MultivariateNormalDistribution.nextMultivariateNormalPrecision(mean, scaledPrecision);
        }
//        System.out.println("\t"+mean[0]+"\t"+mean[1]+"\r");

        double[] result = new double[dim];
        for (int i = 0; i < dim; i++)
            result[i] = mean[i];
        return result;
    }

    public static double[] parseVariableLengthDoubleArray(String inString) throws Arguments.ArgumentException {

        List<Double> returnList = new ArrayList<Double>();
        StringTokenizer st = new StringTokenizer(inString, ",");
        while (st.hasMoreTokens()) {
            try {
                returnList.add(Double.parseDouble(st.nextToken()));
            } catch (NumberFormatException e) {
                throw new Arguments.ArgumentException();
            }

        }

        if (returnList.size() > 0) {
            double[] doubleArray = new double[returnList.size()];
            for (int i = 0; i < doubleArray.length; i++)
                doubleArray[i] = returnList.get(i);

            return doubleArray;
        }
        return null;
    }

    @Override
    public String getDimensionName(int dim) {
        if (getDimension() == 1) {
            return getStatisticName();
        } else {
            return getStatisticName() +".height"+ heightLowers[dim];
        }
    }

    public static void reverse(double[] array) {
        if (array == null) {
            return;
        }
        int i = 0;
        int j = array.length - 1;
        double tmp;
        while (j > i) {
            tmp = array[j];
            array[j] = array[i];
            array[i] = tmp;
            j--;
            i++;
        }
    }

    public static double[] extractUnique(double[] array){
        Set<Double> tmp = new LinkedHashSet<Double>();
        for (Double each : array) {
            tmp.add(each);
        }
        double [] output = new double[tmp.size()];
        int i = 0;
        for (Double each : tmp) {
            output[i++] = each;
        }
        return output;
    }

    public History setUpHistory(String historyString, int nodeState, int parentNodeState, double timeLow, double timeUp){
        double[] heights;
        String[] states;
        if (historyString.equals("{}")){
            heights = new double[]{timeUp,timeLow};
            states = new String[]{getState(nodeState)};
//            returnHistory = new History(heights,states);
        } else {
            List<String> returnList = new ArrayList<String>();
            StringTokenizer st = new StringTokenizer(historyString, "},{");
            while (st.hasMoreTokens()) {
                String test = st.nextToken();
//                returnList.add(st.nextToken());
                returnList.add(test);
//                System.out.println(test);
            }

            int numberOfJumps = returnList.size()/NUMBER_OF_HISTORY_ENTRIES;
            String[][] jumpStrings = new String[numberOfJumps][NUMBER_OF_HISTORY_ENTRIES];
            for (int a = 0; a < numberOfJumps; a++){
                jumpStrings[a][0] = returnList.get(a*NUMBER_OF_HISTORY_ENTRIES);
                jumpStrings[a][1] = returnList.get(a*NUMBER_OF_HISTORY_ENTRIES + 1);
                jumpStrings[a][2] = returnList.get(a*NUMBER_OF_HISTORY_ENTRIES + 2);
            }

            //sorting jumpStrings not necessary: jumps are in order of their occurrence
            //fill heights and states
            heights = new double[numberOfJumps+2];
            states = new String[numberOfJumps+1];
            for (int b = 0; b < numberOfJumps; b++){
                states[b] = jumpStrings[b][1];
                heights[b + 1] = Double.valueOf(jumpStrings[b][0]);
            }

            //sanity check
            if (!jumpStrings[0][1].equals(getState(parentNodeState))){
                System.out.println(jumpStrings[0][1]+"\t"+getState(parentNodeState));
                System.err.println("mismatch in jump history and parent node state");
                System.exit(-1);
            }

            //sanity check
            states[numberOfJumps] = jumpStrings[numberOfJumps-1][2];
            if (!jumpStrings[numberOfJumps-1][2].equals(getState(nodeState))){
                System.err.println("mismatch in jump history and node state");
                System.exit(-1);
            }

            heights[0] = timeUp;
            heights[numberOfJumps+1] =  timeLow;

        }
//        System.out.print("\rhistory ");
//        for (int q =0; q < states.length; q++){
//            System.out.print(heights[q] +"\t"+ states[q] +"\t");
//        }
//        System.out.println(heights[states.length]+"\r");
        return new History(heights,states);
    }

    private String getState(int stateInt){
        String returnString = null;
        try{
            returnString = markovJumpLikelihood.formattedState(new int[] {stateInt}).replaceAll("\"","");
        } catch (IndexOutOfBoundsException iobe) {
            System.err.println("no state found for int = "+stateInt+"...");
            System.exit(-1);
        }
        return returnString;
    }

    public double[] getStateTimeAndDistanceFromRoot(MutableTreeModel tree, NodeRef node, double timeLow, AbstractMultivariateTraitLikelihood traitLikelihood, String traitName, double[] traitLow, double[] precision, BranchRateModel branchRates, boolean useGreatCircleDistance){

        NodeRef nodeOfInterest = node;

        double[] timeDistance = new double[]{0,0};

        double[] rootTrait = traitLikelihood.getTraitForNode(tree, tree.getRoot(), traitName);

        int counter = 0;
        while (nodeOfInterest != tree.getRoot()){
            NodeRef parentNode = tree.getParent(nodeOfInterest);
            History history = setUpHistory(markovJumpLikelihood.getHistoryForNode(tree,nodeOfInterest,SITE),markovJumpLikelihood.getStatesForNode(tree,nodeOfInterest)[SITE],markovJumpLikelihood.getStatesForNode(tree,parentNode)[SITE],tree.getNodeHeight(nodeOfInterest),tree.getNodeHeight(parentNode));
            if (counter == 0){
                if (timeLow > tree.getNodeHeight(nodeOfInterest)){
                    history.truncateLower(timeLow);
                }
            }
            double rate = (branchRates != null ? branchRates.getBranchRate(tree, nodeOfInterest) : 1.0);
            double[] parentTrait = traitLikelihood.getTraitForNode(tree, parentNode, traitName);
            double[] nodeTrait = traitLow;
            if (counter > 0){
                nodeTrait = traitLikelihood.getTraitForNode(tree, nodeOfInterest, traitName);
            }
            history.setTraitsforHeights(parentTrait, nodeTrait, precision, rate, trueNoise);

            timeDistance[0] += history.getStateTime(stateString);

            if (useGreatCircleDistance){
                timeDistance[1] += history.getStateDifferenceInGreatCircleDistanceFromRoot(stateString,rootTrait);
            }  else {
                timeDistance[1] += history.getStateDifferenceInNativeDistanceFromRoot(stateString,rootTrait);
            }

            nodeOfInterest = tree.getParent(nodeOfInterest);
            counter++;
        }

        return timeDistance;
    }

    public double[] getTimeAndDistanceFromRoot(MutableTreeModel tree, NodeRef node, double timeLow, AbstractMultivariateTraitLikelihood traitLikelihood, String traitName, double[] traitLow, boolean useGreatCircleDistance){

        NodeRef nodeOfInterest = node;
        double[] timeDistance = new double[]{0,0};

        double[] rootTrait = traitLikelihood.getTraitForNode(tree, tree.getRoot(), traitName);

        int counter = 0;
        while (nodeOfInterest != tree.getRoot()){
            NodeRef parentNode = tree.getParent(nodeOfInterest);
            double[] parentTrait = traitLikelihood.getTraitForNode(tree, parentNode, traitName);
            double[] nodeTrait = traitLow;
            double nodeHeight = timeLow;
            if (counter > 0){
                nodeTrait = traitLikelihood.getTraitForNode(tree, nodeOfInterest, traitName);
                nodeHeight = tree.getNodeHeight(nodeOfInterest);
            }

            timeDistance[0] += tree.getNodeHeight(parentNode) - nodeHeight;

            if (useGreatCircleDistance){
                timeDistance[1] += getGreatCircleDistance(nodeTrait, rootTrait) - getGreatCircleDistance(parentTrait, rootTrait);
            }  else {
                timeDistance[1] += getNativeDistance(nodeTrait,rootTrait) - getNativeDistance(parentTrait,rootTrait);
            }

            nodeOfInterest = tree.getParent(nodeOfInterest);
            counter++;
        }

        return timeDistance;
    }


    public boolean inClade(MutableTreeModel tree, NodeRef node, TaxonList taxonList) throws TreeUtils.MissingTaxonException {

        Set leafSubSet;
        leafSubSet = TreeUtils.getLeavesForTaxa(tree, taxonList);
        NodeRef mrca = TreeUtils.getCommonAncestorNode(tree, leafSubSet);
        Set mrcaLeafSet =  TreeUtils.getDescendantLeaves(tree,mrca);

        Set nodeLeafSet =  TreeUtils.getDescendantLeaves(tree,node);

        if (!nodeLeafSet.isEmpty()){
            nodeLeafSet.removeAll(mrcaLeafSet);
        }

        if (nodeLeafSet.isEmpty()){
            return true;
        }  else {

        }
        return false;
    }
    private static boolean onAncestralPathTaxa(Tree tree, NodeRef node, TaxonList taxonList) throws TreeUtils.MissingTaxonException {

        if (tree.isExternal(node)) return false;

        Set leafSet = TreeUtils.getDescendantLeaves(tree, node);
        int size = leafSet.size();

        Set targetSet = TreeUtils.getLeavesForTaxa(tree, taxonList);
        leafSet.retainAll(targetSet);

        if (leafSet.size() > 0) {

            // if all leaves below are in target then check just above.
            if (leafSet.size() == size) {

                Set superLeafSet = TreeUtils.getDescendantLeaves(tree, tree.getParent(node));
                superLeafSet.removeAll(targetSet);

                // the branch is on ancestral path if the super tree has some non-targets in it
                return (superLeafSet.size() > 0);

            } else return true;

        } else return false;
    }

    //the sum of the branchLength for all the descendent nodes for a particular node should be larger than a user-specified value
    private static boolean onAncestralPathTime(Tree tree, NodeRef node, double time) {

        double maxDescendentTime = 0;

        Set leafSet = TreeUtils.getExternalNodes(tree, node);
        Set nodeSet = TreeUtils.getExternalNodes(tree, node);

        Iterator iter = leafSet.iterator();

        while (iter.hasNext()) {
//            System.out.println("found node set");
            NodeRef currentNode = (NodeRef)iter.next();

            while (tree.getNodeHeight(node) > tree.getNodeHeight(currentNode)) {
//                System.out.println("found node height");
                if (!nodeSet.contains(currentNode)) {
//                    System.out.println("found node");
                    nodeSet.add(currentNode);
                }
                currentNode = tree.getParent(currentNode);
            }
        }

        Iterator nodeIter = nodeSet.iterator();

        while (nodeIter.hasNext()) {
            NodeRef testNode = (NodeRef)nodeIter.next();
            maxDescendentTime += tree.getBranchLength(testNode);
        }

        if (maxDescendentTime > time){
            return true;
        }   else {
            return false;
        }
    }

    private static double getAreaFrom2Dtraits(List<double[]> traits2D,  double hpdValue){

        boolean bandwidthlimit = true;

        double totalArea = 0;
        double[][] y = new double[2][traits2D.size()];
        for (int a=0; a<traits2D.size(); a++){
            double[] trait = traits2D.get(a);
            y[0][a]= trait[0];
            y[1][a]= trait[1];
//            System.err.println(trait[0]+"\t"+trait[1]);
        }

        ContourMaker contourMaker;
        contourMaker = new ContourWithSynder(y[0], y[1], bandwidthlimit);

        ContourPath[] paths = contourMaker.getContourPaths(hpdValue);
        int pathCounter = 1;
        for (ContourPath path : paths) {
            KMLCoordinates coords = new KMLCoordinates(path.getAllX(), path.getAllY());
            Element testElement = new Element("test");
            testElement.addContent(coords.toXML());
            Polygon2D testPolygon = new Polygon2D(testElement);
            totalArea += testPolygon.calculateArea();
//            System.err.println("area: "+testPolygon.calculateArea());
        }

        return totalArea;
    }

    private static double[] convertDoubles(List<Double> doubles) {
        double[] ret = new double[doubles.size()];
        Iterator<Double> iterator = doubles.iterator();
        int i = 0;
        while(iterator.hasNext())
        {
            ret[i] = iterator.next();
            i++;
        }
        return ret;
    }

    private static double getSpearmanRho(double[] data1, double[] data2){

        double data1Ranks[] = new NaturalRanking().rank(data1);
        double data2Ranks[] = new NaturalRanking().rank(data2);

        int counter = 0;
        double d_i = 0;
        while(counter < data1Ranks.length){
            d_i += Math.pow(data1Ranks[counter] -  data2Ranks[counter], 2);
            counter ++;
        }
        return  (1 - (6*d_i)/(data1Ranks.length*(Math.pow(data1Ranks.length,2) - 1)));
    }



//    private int getStateInt(String state){
//        int returnInt = -1;
//        int counter = 0;
//        try{
//            while (returnInt < 0) {
//                if (state.equalsIgnoreCase((markovJumpLikelihood.formattedState(new int[] {counter})).replaceAll("\"",""))) {
//                    returnInt = counter;
//                }
//                counter ++;
//            }
//        } catch (IndexOutOfBoundsException iobe) {
//            int states[] = new int[counter];
//            for (int a = 0; a < states.length; a++){
//                states[a] = a;
//            }
//            System.err.println("state "+state+" not found among "+markovJumpLikelihood.formattedState(states)+ "... ignoring state");
//            System.exit(-1);
//        }
//        return returnInt;
//    }

    enum Mode {
        AVERAGE,
        WEIGHTED_AVERAGE,
        MEDIAN,
        COEFFICIENT_OF_VARIATION,
        SPEARMAN,
        CORRELATION_COEFFICIENT,
        R_SQUARED
    }

    enum summaryStatistic {
        TRAIT,
        TRAIT2DAREA,
        DIFFUSION_TIME,
        DIFFUSION_DISTANCE,
        DIFFUSION_RATE,
        DIFFUSION_COEFFICIENT,
        WAVEFRONT_DISTANCE,
        WAVEFRONT_DISTANCE_PHYLO,
        WAVEFRONT_RATE,
        DISTANCE_TIME_CORRELATION
    }

    enum BranchSet {
        ALL,
        CLADE,
        BACKBONE, //TODO: to implement
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return CONTINUOUS_DIFFUSION_STATISTIC;
        }

        @Override
        public String[] getParserNames() {
            return new String[]{getParserName(), DIFFUSION_RATE_STATISTIC, TREE_DISPERSION_STATISTIC};
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String name = xo.getAttribute(NAME, xo.getId());

            boolean greatCircleDistances = xo.getAttribute(USE_GREATCIRCLEDISTANCES, false); // Default value is false
            Mode statMode;
            String mode = xo.getAttribute(MODE, WEIGHTED_AVERAGE);

            if (mode.equals(AVERAGE)) {
                statMode = Mode.AVERAGE;
            } else if (mode.equals(MEDIAN)) {
                statMode = Mode.MEDIAN;
            } else if (mode.equals(COEFFICIENT_OF_VARIATION)) {
                statMode = Mode.COEFFICIENT_OF_VARIATION;
            } else if (mode.equals(WEIGHTED_AVERAGE)) {
                statMode = Mode.WEIGHTED_AVERAGE;
            } else if (mode.equals(SPEARMAN)) {
                statMode = Mode.SPEARMAN;
            } else if (mode.equals(CORRELATION_COEFFICIENT)) {
                statMode = Mode.CORRELATION_COEFFICIENT;
            } else if (mode.equals(R_SQUARED)) {
                statMode = Mode.R_SQUARED;
            } else {
                System.err.println("Unknown mode: "+mode+". Reverting to weighted average for "+name);
                statMode = Mode.WEIGHTED_AVERAGE;
            }

            final double upperHeight = xo.getAttribute(HEIGHT_UPPER, Double.MAX_VALUE);
            final double lowerHeight = xo.getAttribute(HEIGHT_LOWER, 0.0);

            double[] lowerHeights = null;
            if (xo.hasAttribute(HEIGHT_LOWER_SERIE)){
                String lowerHeightsString = xo.getStringAttribute(HEIGHT_LOWER_SERIE);
                try {
                    lowerHeights = parseVariableLengthDoubleArray(lowerHeightsString);
                } catch (Arguments.ArgumentException e) {
                    System.err.println(name+": error reading " + HEIGHT_LOWER_SERIE);
                    System.exit(1);
                }
            }

            boolean cumulative = xo.getAttribute(CUMULATIVE, false);

            boolean trueNoise = xo.getAttribute(NOISE, false); // Default value is false

//            boolean diffCoeff = xo.getAttribute(BOOLEAN_DC_OPTION, false); // Default value is false
            summaryStatistic summaryStat;
            String statistic = xo.getAttribute(STATISTIC, DIFFUSION_RATE);
            int dimension = 1;
            if (statistic.equals(DIFFUSION_RATE)) {
                summaryStat = summaryStatistic.DIFFUSION_RATE;
                if (mode.equals(SPEARMAN) || mode.equals(R_SQUARED) || mode.equals(CORRELATION_COEFFICIENT)){
                    System.err.println(name+": mode = "+mode+" ignored for "+DIFFUSION_TIME+", reverting to weighted average mode");
                    statMode = Mode.WEIGHTED_AVERAGE;
                }
            } else if (statistic.equals(DIFFUSION_TIME)) {
                summaryStat = summaryStatistic.DIFFUSION_TIME;
                if (!mode.equals(WEIGHTED_AVERAGE)) {
                    System.err.println(name+": mode = "+mode+" ignored for "+DIFFUSION_TIME);
                }
            } else if (statistic.equals(DIFFUSION_DISTANCE)) {
                summaryStat = summaryStatistic.DIFFUSION_DISTANCE;
                if (!mode.equals(WEIGHTED_AVERAGE)) {
                    System.err.println(name+": mode = "+mode+" ignored for "+DIFFUSION_DISTANCE);
                }
            } else if (statistic.equals(DISTANCE_TIME_CORRELATION)) {
                summaryStat = summaryStatistic.DISTANCE_TIME_CORRELATION;
                if (mode.equals(AVERAGE) || mode.equals(WEIGHTED_AVERAGE) || mode.equals(COEFFICIENT_OF_VARIATION) || mode.equals(MEDIAN)){
                    System.err.println(name+": mode = "+mode+" ignored for "+DISTANCE_TIME_CORRELATION+", reverting to correlation coefficient mode");
                    statMode = Mode.CORRELATION_COEFFICIENT;
                }
            } else if (statistic.equals(WAVEFRONT_DISTANCE)) {
                summaryStat = summaryStatistic.WAVEFRONT_DISTANCE;
                if (!mode.equals(WEIGHTED_AVERAGE)) {
                    System.err.println(name + ": mode = " + mode + " ignored for " + WAVEFRONT_DISTANCE);
                }
            } else if (statistic.equals(WAVEFRONT_DISTANCE_PHYLO)) {
                summaryStat = summaryStatistic.WAVEFRONT_DISTANCE_PHYLO;
                if (!mode.equals(WEIGHTED_AVERAGE)) {
                    System.err.println(name + ": mode = " + mode + " ignored for " + WAVEFRONT_DISTANCE);
                }
            } else if (statistic.equals(TRAIT)) {
                summaryStat = summaryStatistic.TRAIT;
                if (mode.equals(WEIGHTED_AVERAGE)) {
                    System.err.println(name + ": mode = " + mode + " ignored for " + TRAIT + ", resorting to " + AVERAGE);
                    statMode = Mode.AVERAGE;
                }
                if (upperHeight < Double.MAX_VALUE) {
                    System.err.println(name + ": only " + HEIGHT_LOWER + " or " + HEIGHT_LOWER_SERIE + " are relevant for " + TRAIT);
                }
                dimension = xo.getAttribute(DIMENSION, 1);
                if (dimension == 0) {
                    System.err.println(name + ": trait dimensions start from 1. Setting dimension to 1");
                    dimension = 1;
                }
                if (cumulative) {
                    System.err.println(name + ": " + CUMULATIVE + " is ignored for " + TRAIT);
                }
                if (greatCircleDistances) {
                    System.err.println(name + ": " + USE_GREATCIRCLEDISTANCES + " is ignored for " + TRAIT);
                }
            } else if (statistic.equals(TRAIT2DAREA)) {
                summaryStat = summaryStatistic.TRAIT2DAREA;
                dimension = xo.getAttribute(DIMENSION, 2);
                if (dimension != 2){
                    System.err.println(name + ": trait dimension ("+dimension+") is not 2. Cannot calculate 2D area for the traits, 0's will be returned");
                }
          } else if (statistic.equals(WAVEFRONT_RATE)) {
                summaryStat = summaryStatistic.WAVEFRONT_RATE;
            } else if (statistic.equals(DIFFUSION_COEFFICIENT)) {
                summaryStat = summaryStatistic.DIFFUSION_COEFFICIENT;
            } else if (statistic.equals(DISTANCE_TIME_CORRELATION)) {
                summaryStat = summaryStatistic.DISTANCE_TIME_CORRELATION;
            } else {
                System.err.println(name+": unknown statistic: "+statistic+". Reverting to diffusion rate.");
                summaryStat = summaryStatistic.DIFFUSION_RATE;
            }

            BranchSet branchset;
            String branchMode = xo.getAttribute(BRANCHSET, ALL);
            if (branchMode.equals(CLADE)) {
                branchset = BranchSet.CLADE;
            } else if (branchMode.equals(BACKBONE)) {
                branchset = BranchSet.BACKBONE;
            } else if (branchMode.equals(ALL)) {
                branchset = BranchSet.ALL;
            } else {
                System.err.println(name+": unknown branchset: "+branchMode+". Reverting to all branches.");
                branchset = BranchSet.ALL;
            }

            TaxonList taxonList = null;
            double backboneTime = 0;
            if  (branchset.equals(BranchSet.CLADE)){
                taxonList = (TaxonList) xo.getChild(TaxonList.class);
                if (taxonList==null){
                    System.err.println("empty taxon list in continuousDiffusionStatistic despite 'clade' branchSet attribute");
                }
            } else if (branchset.equals(BranchSet.BACKBONE)){
                taxonList = (TaxonList) xo.getChild(TaxonList.class);
                if (xo.hasAttribute(BACKBONE_TIME)){
                    backboneTime = xo.getAttribute(BACKBONE_TIME, 0.0);
                    if (taxonList!=null){
                        System.err.println("both backbone time and taxon list provided for backbone definition in continuousDiffusionStatistic. Ignoring taxon list...");
                    }
                }  else if (taxonList==null){
                    System.err.println("empty taxon list and no backboneTime in continuousDiffusionStatistic despite 'backbone' branchSet attribute. Ignoring 'backbone' branchSet...");
                }
            } else if (branchset.equals(BranchSet.ALL)){
                taxonList = (TaxonList) xo.getChild(TaxonList.class);
                if (taxonList!=null){
                    System.err.println("taxon list provided in continuousDiffusionStatistic but no 'clade' or 'backbone' branchSet attribute?? Ignoring taxon list...");
                }
                if (xo.hasAttribute(BACKBONE_TIME)){
                    System.err.println("backoneTime provided in continuousDiffusionStatistic but no 'backbone' branchSet attribute?? Ignoring backboneTime list...");
                }
            }


            String stateString = null;
            if (xo.hasAttribute(DISCRETE_STATE)){
                stateString = xo.getStringAttribute(DISCRETE_STATE);
            }

            List<AbstractMultivariateTraitLikelihood> traitLikelihoods = new ArrayList<AbstractMultivariateTraitLikelihood>();
            MarkovJumpsBeagleTreeLikelihood mjtl = null;

            for (int i = 0; i < xo.getChildCount(); i++) {
//                System.err.println("child is = "+xo.getChildName(i));
                if (xo.getChild(i) instanceof AbstractMultivariateTraitLikelihood) {
                    AbstractMultivariateTraitLikelihood amtl = (AbstractMultivariateTraitLikelihood) xo.getChild(i);
                    traitLikelihoods.add(amtl);
                }
                if (xo.getChild(i) instanceof MarkovJumpsBeagleTreeLikelihood) {
                    mjtl = (MarkovJumpsBeagleTreeLikelihood) xo.getChild(i);
                }
            }

            if (stateString == null && mjtl != null) {
                System.err.println(name+": markovJumpsTreeLikelihood specified for state-specific summaries but no state string.. ignoring markovJumpsTreeLikelihood");
                mjtl = null;
            }  else if (stateString != null && mjtl == null){
                System.err.println(name+": state number provided for state-specific summaries but no markovJumpsTreeLikelihood specified.. ignoring state");
                stateString = null;
            }  else if (stateString != null && mjtl != null) {
                if (statistic.equals(TRAIT)){
                    System.err.println(name+": ignoring state-specific summary (for "+stateString+") for " + TRAIT+", resorting to overall summary");
                } else {
                    int stateInt = -1;
                    int counter = 0;
                    try{
                        while (stateInt < 0) {
                            if (stateString.equalsIgnoreCase((mjtl.formattedState(new int[] {counter})).replaceAll("\"",""))) {
                                stateInt = counter;
                                System.out.println(name+": summarizing continuous diffusion statistic for state "+mjtl.formattedState(new int[] {counter}));
                            }
                            counter ++;
                        }
                    } catch (IndexOutOfBoundsException iobe) {
                        int states[] = new int[counter];
                        for (int a = 0; a < states.length; a++){
                            states[a] = a;
                        }
                        System.err.println(name+": state "+stateString+" not found among "+mjtl.formattedState(states)+ "... ignoring state");
                        mjtl = null;
                        stateString = null;
                    }
                }
            }



            return new ContinuousDiffusionStatistic(name, traitLikelihoods, greatCircleDistances, statMode, summaryStat, upperHeight, lowerHeight, lowerHeights, cumulative, trueNoise, dimension, taxonList, branchset, backboneTime, stateString, mjtl);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A statistic that returns the average of the branch diffusion rates";
        }

        public Class getReturnType() {
            return TreeStatistic.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newStringRule(NAME, true),
                AttributeRule.newBooleanRule(USE_GREATCIRCLEDISTANCES, true),
                AttributeRule.newStringRule(MODE, true),
                AttributeRule.newStringRule(STATISTIC,true),
                AttributeRule.newStringRule(DISCRETE_STATE,true),
                AttributeRule.newDoubleRule(HEIGHT_UPPER, true),
                AttributeRule.newDoubleRule(HEIGHT_LOWER, true),
                AttributeRule.newStringRule(HEIGHT_LOWER_SERIE, true),
                AttributeRule.newDoubleRule(DIMENSION, true),
                AttributeRule.newBooleanRule(CUMULATIVE, true),
                AttributeRule.newBooleanRule(NOISE, true),
                AttributeRule.newStringRule(BRANCHSET, true),
                new ElementRule(TaxonList.class,true),
                new ElementRule(AbstractMultivariateTraitLikelihood.class, 1, Integer.MAX_VALUE),
                new ElementRule(MarkovJumpsBeagleTreeLikelihood.class, true)
        };
    };

    private boolean useGreatCircleDistances;
    private List<AbstractMultivariateTraitLikelihood> traitLikelihoods;
    private MarkovJumpsBeagleTreeLikelihood markovJumpLikelihood;
    //    private int stateInt;
    private String stateString;
    private Mode summaryMode;
    private summaryStatistic summaryStat;
    private double heightUpper;
    private double[] heightLowers;
    private boolean cumulative;
    private boolean trueNoise;
    private int dimension;
    private TaxonList taxonList;
    private BranchSet branchset;
    private double backboneTime;

    private class History {

        private double[] historyHeights;
        private String[] historyStates;
        private double[][] historyTraits;

        public History(double historyHeights[], String historyStates[]) {
            this.historyHeights = historyHeights;
            this.historyStates = historyStates;
        }

        public void truncateUpper(double time) {
            int cutFrom = -1;
            for (int a = 0; a < (historyHeights.length - 1); a++) {
                if ((time < historyHeights[a]) && (time > historyHeights[a+1])) {
                    cutFrom = a;
                }
            }

            if (cutFrom < 0){
                System.err.println("no upper truncation of discrete trait history on branch possible");
                System.exit(0);
            }

            double[] tempHeights = new double[historyHeights.length - cutFrom];
            String[] tempStates = new String[historyStates.length - cutFrom];

            tempHeights = Arrays.copyOfRange(historyHeights, cutFrom, historyHeights.length);
            tempHeights[0] = time;

            tempStates = Arrays.copyOfRange(historyStates, cutFrom, historyStates.length);

            historyHeights = tempHeights;
            historyStates =  tempStates;
        }

        public void truncateLower (double time) {
            int cutTo = -1;

            for (int a = (historyHeights.length - 1); a > 0; a--) {
                if ((time > historyHeights[a]) && (time < historyHeights[a - 1])) {
                    cutTo = a;
                }
            }

            if (cutTo < 0){
                System.err.println("no lower truncation of discrete trait history on branch possible");
                System.exit(0);
            }

            double[] tempHeights = new double[cutTo + 1];
            String[] tempStates = new String[cutTo];

            tempHeights = Arrays.copyOfRange(historyHeights, 0, cutTo + 1);
            tempHeights[(tempHeights.length-1)] = time;

            tempStates = Arrays.copyOfRange(historyStates, 0, cutTo);

            historyHeights = tempHeights;
            historyStates =  tempStates;
        }

        public double getStateTime(String state){
            double time = 0;
            for (int x = 0; x < historyStates.length; x++){
                if (state.equals(historyStates[x])){
                    time += (historyHeights[x] - historyHeights[x+1]);
                }
            }
            return time;
        }

        private void setTraitsforHeights(double[] traitUp,  double[] traitLow, double[] precisionArray, double rate, boolean trueNoise){
            historyTraits = new double[historyHeights.length][2];
            for (int x = 0; x < historyHeights.length; x++){
                if(x == 0){
                    historyTraits[x] = traitUp;
                } else if (x == (historyTraits.length -1)) {
                    historyTraits[x] = traitLow;
                } else {
                    historyTraits[x] = imputeValue(traitUp, traitLow, historyHeights[x], historyHeights[(historyHeights.length -1)], historyHeights[0], precisionArray, rate, trueNoise);
                }

            }
        }

        public double getStateGreatCircleDistance(String state){
            double distance = 0;
            for (int x = 0; x < historyStates.length; x++){
                if (state.equals(historyStates[x])){
                    distance += getGreatCircleDistance(historyTraits[x],historyTraits[x+1]);
                }
            }
            return distance;
        }

        public double getStateDifferenceInGreatCircleDistanceFromRoot(String state, double[] rootTrait){
            double distance = 0;
            for (int x = 0; x < historyStates.length; x++){
                if (state.equals(historyStates[x])){
                    distance += (getGreatCircleDistance(historyTraits[x+1],rootTrait) - getGreatCircleDistance(historyTraits[x],rootTrait));
                }
            }
            return distance;
        }

        public double getStateNativeDistance(String state){
            double distance = 0;
            for (int x = 0; x < historyStates.length; x++){
                if (state.equals(historyStates[x])){
                    distance += getNativeDistance(historyTraits[x],historyTraits[x+1]);
                }
            }
            return distance;
        }

        public double getStateDifferenceInNativeDistanceFromRoot(String state, double[] rootTrait){
            double distance = 0;
            for (int x = 0; x < historyStates.length; x++){
                if (state.equals(historyStates[x])){
                    distance += (getNativeDistance(historyTraits[x+1],rootTrait) - getNativeDistance(historyTraits[x],rootTrait));
                }
            }
            return distance;
        }

    }

}

