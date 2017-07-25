/*
 * DiffusionRateStatistic.java
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

import dr.app.util.Arguments;
import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeStatistic;
import dr.geo.math.SphericalPolarCoordinates;
import dr.inference.model.Statistic;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.stats.DiscreteStatistics;
import dr.xml.*;

import java.util.*;

/**
 * @author Marc Suchard
 * @author Philippe Lemey
 * @author Andrew Rambaut
 */
@Deprecated
public class DiffusionRateStatistic extends Statistic.Abstract {

    public static final String DIFFUSION_RATE_STATISTIC = "diffusionRateStatistic";
    public static final String TREE_DISPERSION_STATISTIC = "treeDispersionStatistic";
    public static final String BOOLEAN_DIS_OPTION = "greatCircleDistance";
    public static final String MODE = "mode";
    public static final String MEDIAN = "median";
    public static final String AVERAGE = "average";  // average over all branches
    public static final String WEIGHTED_AVERAGE = "weightedAverage"; // weighted average (=total distance/total time)
    public static final String COEFFICIENT_OF_VARIATION = "coefficientOfVariation"; // weighted average (=total distance/total time)
    public static final String STATISTIC = "statistic";
    public static final String DIFFUSION_RATE = "diffusionRate"; // weighted average (=total distance/total time)
    public static final String WAVEFRONT_DISTANCE = "wavefrontDistance"; // weighted average (=total distance/total time)
    public static final String WAVEFRONT_RATE = "wavefrontRate"; // weighted average (=total distance/total time)
    public static final String DIFFUSION_COEFFICIENT = "diffusionCoefficient";
    //    public static final String DIFFUSIONCOEFFICIENT = "diffusionCoefficient"; // weighted average (=total distance/total time)
    //    public static final String BOOLEAN_DC_OPTION = "diffusionCoefficient";
    public static final String HEIGHT_UPPER = "heightUpper";
    public static final String HEIGHT_LOWER = "heightLower";
    public static final String HEIGHT_LOWER_SERIE = "heightLowerSerie";
    public static final String CUMULATIVE = "cumulative";

    public DiffusionRateStatistic(String name, List<AbstractMultivariateTraitLikelihood> traitLikelihoods,
                                        boolean option, Mode mode,
                                        summaryStatistic statistic, double heightUpper, double heightLower,
                                        double[] lowerHeights, boolean cumulative) {
        super(name);
        this.traitLikelihoods = traitLikelihoods;
        this.useGreatCircleDistances = option;
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
    }

    public int getDimension() {
        return heightLowers.length;
    }

    public double getStatisticValue(int dim) {

        String traitName = traitLikelihoods.get(0).getTraitName();

        double treelength = 0;
        double treeDistance = 0;
        double maxDistanceFromRoot = 0;
        double maxDistanceOverTimeFromRoot = 0;

        //double[] rates =  null;
        List<Double> rates = new ArrayList<Double>();
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

//        System.out.println("dim = "+dim+", heightLower = "+lowerHeight+", heightUpper = "+upperHeight);

        for (AbstractMultivariateTraitLikelihood traitLikelihood : traitLikelihoods) {
            MutableTreeModel tree = traitLikelihood.getTreeModel();
            BranchRateModel branchRates = traitLikelihood.getBranchRateModel();

            for (int i = 0; i < tree.getNodeCount(); i++) {
                NodeRef node = tree.getNode(i);

                if (node != tree.getRoot()) {

                    NodeRef parentNode = tree.getParent(node);
                    if ((tree.getNodeHeight(parentNode) > lowerHeight) && (tree.getNodeHeight(node) < upperHeight)) {

                        double[] trait = traitLikelihood.getTraitForNode(tree, node, traitName);
                        double[] parentTrait = traitLikelihood.getTraitForNode(tree, parentNode, traitName);

                        double[] traitUp = parentTrait;
                        double[] traitLow = trait;

                        double timeUp = tree.getNodeHeight(parentNode);
                        double timeLow = tree.getNodeHeight(node);

                        double rate = (branchRates != null ? branchRates.getBranchRate(tree, node) : 1.0);

                        MultivariateDiffusionModel diffModel = traitLikelihood.diffusionModel;
                        double[] precision = diffModel.getPrecisionParameter().getParameterValues();

                        if (tree.getNodeHeight(parentNode) > upperHeight) {
                            timeUp = upperHeight;
                            //TODO: implement TrueNoise??
                            traitUp = imputeValue(trait, parentTrait, upperHeight, tree.getNodeHeight(node), tree.getNodeHeight(parentNode), precision, rate, false);
                        }

                        if (tree.getNodeHeight(node) < lowerHeight) {
                            timeLow = lowerHeight;
                            traitLow = imputeValue(trait, parentTrait, lowerHeight, tree.getNodeHeight(node), tree.getNodeHeight(parentNode), precision, rate, false);
                        }

                        double time = timeUp - timeLow;
                        treelength += time;

                        double[] rootTrait = traitLikelihood.getTraitForNode(tree, tree.getRoot(), traitName);

                        if (useGreatCircleDistances && (trait.length == 2)) { // Great Circle distance
                            SphericalPolarCoordinates coord1 = new SphericalPolarCoordinates(traitLow[0], traitLow[1]);
                            SphericalPolarCoordinates coord2 = new SphericalPolarCoordinates(traitUp[0], traitUp[1]);
                            double distance = coord1.distance(coord2);
                            treeDistance += distance;
                            double dc = Math.pow(distance,2)/(4*time);
                            diffusionCoefficients.add(dc);
                            waDiffusionCoefficient +=  dc*time;
                            rates.add(distance/time);

                            SphericalPolarCoordinates rootCoord = new SphericalPolarCoordinates(rootTrait[0], rootTrait[1]);
                            double tempDistanceFromRoot = rootCoord.distance(coord2);
                            if (tempDistanceFromRoot > maxDistanceFromRoot){
                                maxDistanceFromRoot = tempDistanceFromRoot;
                                maxDistanceOverTimeFromRoot = tempDistanceFromRoot/(tree.getNodeHeight(tree.getRoot()) - timeLow);
                                //distance between traitLow and traitUp for maxDistanceFromRoot
                                if (timeUp == upperHeight) {
                                    maxDistanceFromRoot = distance;
                                    maxDistanceOverTimeFromRoot = distance/time;
                                }
                            }


                        } else {
                            double distance = getNativeDistance(traitLow, traitUp);
                            treeDistance += distance;
                            double dc = Math.pow(distance,2)/(4*time);
                            diffusionCoefficients.add(dc);
                            waDiffusionCoefficient += dc*time;
                            rates.add(distance/time);

                            double tempDistanceFromRoot = getNativeDistance(traitLow, rootTrait);
                            if (tempDistanceFromRoot > maxDistanceFromRoot){
                                maxDistanceFromRoot = tempDistanceFromRoot;
                                maxDistanceOverTimeFromRoot = tempDistanceFromRoot/(tree.getNodeHeight(tree.getRoot()) - timeLow);
                                //distance between traitLow and traitUp for maxDistanceFromRoot
                                if (timeUp == upperHeight) {
                                    maxDistanceFromRoot = distance;
                                    maxDistanceOverTimeFromRoot = distance/time;
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
                // don't compute mean twice
                final double mean = DiscreteStatistics.mean(toArray(rates));
                return Math.sqrt(DiscreteStatistics.variance(toArray(rates), mean)) / mean;
            } else {
                return treeDistance / treelength;
            }
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
                return waDiffusionCoefficient/treelength;
            }
        }  else if (summaryStat == summaryStatistic.WAVEFRONT_DISTANCE) {
            return maxDistanceFromRoot;
        }  else {
            return maxDistanceOverTimeFromRoot;
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
//            System.out.println(sum);
        }
        return Math.sqrt(sum);
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

    enum Mode {
        AVERAGE,
        WEIGHTED_AVERAGE,
        MEDIAN,
        COEFFICIENT_OF_VARIATION
    }

    enum summaryStatistic {
        DIFFUSION_RATE,
        DIFFUSION_COEFFICIENT,
        WAVEFRONT_DISTANCE,
        WAVEFRONT_RATE,
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return DIFFUSION_RATE_STATISTIC;
        }

        @Override
        public String[] getParserNames() {
            return new String[]{getParserName(), TREE_DISPERSION_STATISTIC};
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String name = xo.getAttribute(NAME, xo.getId());

            boolean option = xo.getAttribute(BOOLEAN_DIS_OPTION, true); // Default value is true
            Mode averageMode;
            String mode = xo.getAttribute(MODE, WEIGHTED_AVERAGE);

            if (mode.equals(AVERAGE)) {
                averageMode = Mode.AVERAGE;
            } else if (mode.equals(MEDIAN)) {
                averageMode = Mode.MEDIAN;
            } else if (mode.equals(COEFFICIENT_OF_VARIATION)) {
                averageMode = Mode.COEFFICIENT_OF_VARIATION;
            } else if (mode.equals(WEIGHTED_AVERAGE)) {
                averageMode = Mode.WEIGHTED_AVERAGE;
            } else {
                System.err.println("Unknown mode: "+mode+". Reverting to weighted average");
                averageMode = Mode.WEIGHTED_AVERAGE;
            }

//            boolean diffCoeff = xo.getAttribute(BOOLEAN_DC_OPTION, false); // Default value is false
            summaryStatistic summaryStat;
            String statistic = xo.getAttribute(STATISTIC, DIFFUSION_RATE);
            if (statistic.equals(DIFFUSION_RATE)) {
                summaryStat = summaryStatistic.DIFFUSION_RATE;
            } else if (statistic.equals(WAVEFRONT_DISTANCE)) {
                summaryStat = summaryStatistic.WAVEFRONT_DISTANCE;
            } else if (statistic.equals(WAVEFRONT_RATE)) {
                summaryStat = summaryStatistic.WAVEFRONT_RATE;
            } else if (statistic.equals(DIFFUSION_COEFFICIENT)) {
                summaryStat = summaryStatistic.DIFFUSION_COEFFICIENT;
            } else {
                System.err.println("Unknown statistic: "+statistic+". Reverting to diffusion rate");
                summaryStat = summaryStatistic.DIFFUSION_COEFFICIENT;
            }

            final double upperHeight = xo.getAttribute(HEIGHT_UPPER, Double.MAX_VALUE);
            final double lowerHeight = xo.getAttribute(HEIGHT_LOWER, 0.0);

            double[] lowerHeights = null;
            if (xo.hasAttribute(HEIGHT_LOWER_SERIE)){
                String lowerHeightsString = xo.getStringAttribute(HEIGHT_LOWER_SERIE);
                try {
                    lowerHeights = parseVariableLengthDoubleArray(lowerHeightsString);
                } catch (Arguments.ArgumentException e) {
                    System.err.println("Error reading " + HEIGHT_LOWER_SERIE);
                    System.exit(1);
                }
            }

            boolean cumulative = xo.getAttribute(CUMULATIVE, false);

            List<AbstractMultivariateTraitLikelihood> traitLikelihoods = new ArrayList<AbstractMultivariateTraitLikelihood>();

            for (int i = 0; i < xo.getChildCount(); i++) {
                if (xo.getChild(i) instanceof AbstractMultivariateTraitLikelihood) {
                    AbstractMultivariateTraitLikelihood amtl = (AbstractMultivariateTraitLikelihood) xo.getChild(i);
                    traitLikelihoods.add(amtl);
                }
            }

            return new DiffusionRateStatistic(name, traitLikelihoods, option, averageMode, summaryStat, upperHeight, lowerHeight, lowerHeights, cumulative);
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
                AttributeRule.newBooleanRule(BOOLEAN_DIS_OPTION, true),
                AttributeRule.newStringRule(MODE, true),
                AttributeRule.newStringRule(STATISTIC,true),
                AttributeRule.newDoubleRule(HEIGHT_UPPER, true),
                AttributeRule.newDoubleRule(HEIGHT_LOWER, true),
                AttributeRule.newStringRule(HEIGHT_LOWER_SERIE,true),
                AttributeRule.newBooleanRule(CUMULATIVE, true),
                new ElementRule(AbstractMultivariateTraitLikelihood.class, 1, Integer.MAX_VALUE),
        };
    };

    private boolean useGreatCircleDistances;
    private List<AbstractMultivariateTraitLikelihood> traitLikelihoods;
    private Mode summaryMode;
    private summaryStatistic summaryStat;
    private double heightUpper;
    private double[] heightLowers;
    private boolean cumulative;

}

