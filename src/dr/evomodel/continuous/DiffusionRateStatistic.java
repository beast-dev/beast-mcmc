/*
 * DiffusionRateStatistic.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.tree.MultivariateTraitTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.DiscretizedBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeStatistic;
import dr.geo.math.SphericalPolarCoordinates;
import dr.inference.model.Statistic;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.stats.DiscreteStatistics;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc Suchard
 * @author Philippe Lemey
 * @author Andrew Rambaut
 */
public class DiffusionRateStatistic extends Statistic.Abstract {

    public static final String DIFFUSION_RATE_STATISTIC = "diffusionRateStatistic";
    public static final String TREE_DISPERSION_STATISTIC = "treeDispersionStatistic";
    public static final String BOOLEAN_DIS_OPTION = "greatCircleDistance";
    public static final String MODE = "mode";
    public static final String MEDIAN = "median";
    public static final String AVERAGE = "average";  // average over all branches
    public static final String WEIGHTEDAVERAGE = "weightedAverage"; // weighted average (=total distance/total time)
    public static final String COEFFICIENT_OF_VARIATION = "coefficientOfVariation"; // weighted average (=total distance/total time)
//    public static final String DIFFUSIONCOEFFICIENT = "diffusionCoefficient"; // weighted average (=total distance/total time)
    public static final String BOOLEAN_DC_OPTION = "diffusionCoefficient";
    public static final String HEIGHT_UPPER = "heightUpper";
    public static final String HEIGHT_LOWER = "heightLower";

    public DiffusionRateStatistic(String name, List<AbstractMultivariateTraitLikelihood> traitLikelihoods,
                                  boolean option, Mode mode, boolean diffusionCoefficient, double heightUpper, double heightLower, DiscretizedBranchRates branchRates) {
        super(name);
        this.traitLikelihoods = traitLikelihoods;
        this.useGreatCircleDistances = option;
        summaryMode =  mode;
        this.diffusionCoefficient = diffusionCoefficient;
        this.heightUpper = heightUpper;
        this.heightLower = heightLower;
        this.branchRates = branchRates;

    }

    public int getDimension() {
        return 1;
    }

    public double getStatisticValue(int dim) {

        String traitName = traitLikelihoods.get(0).getTraitName();

        double treelength = 0;
        double treeDistance = 0;

        //double[] rates =  null;
        List<Double> rates = new ArrayList<Double>();
        //double[] diffusionCoefficients =  null;
        List<Double> diffusionCoefficients = new ArrayList<Double>();
        double waDiffusionCoefficient =  0;

        for (AbstractMultivariateTraitLikelihood traitLikelihood : traitLikelihoods) {
            MultivariateTraitTree tree = traitLikelihood.getTreeModel();

            for (int i = 0; i < tree.getNodeCount(); i++) {
                NodeRef node = tree.getNode(i);

                if (node != tree.getRoot()) {

                    NodeRef parentNode = tree.getParent(node);
                    if ((tree.getNodeHeight(parentNode) > heightLower) && (tree.getNodeHeight(node) < heightUpper)) {

                        double[] trait = traitLikelihood.getTraitForNode(tree, node, traitName);
                        double[] parentTrait = traitLikelihood.getTraitForNode(tree, parentNode, traitName);

                        double[] traitUp = parentTrait;
                        double[] traitLow = trait;

                        double timeUp = tree.getNodeHeight(parentNode);
                        double timeLow = tree.getNodeHeight(node);

                        double rate = 1;
                        if (branchRates != null){
                            rate = branchRates.getBranchRate(tree,node);
                        }

                        MultivariateDiffusionModel diffModel = traitLikelihoods.get(0).diffusionModel;
                        double[] precision = diffModel.getPrecisionParameter().getParameterValues();

                        if (tree.getNodeHeight(parentNode) > heightUpper) {
                            timeUp = heightUpper;
                            //TODO: implement TrueNoise
                            traitUp = imputeValue(trait, parentTrait, heightUpper, tree.getNodeHeight(node), tree.getNodeHeight(parentNode), precision, rate, false);
                        }

                        if (tree.getNodeHeight(node) < heightLower) {
                            timeLow = heightLower;
                            traitLow = imputeValue(trait, parentTrait, heightLower, tree.getNodeHeight(node), tree.getNodeHeight(parentNode), precision, rate, false);
                        }

                        double time = timeUp - timeLow;
                        treelength += time;

                        if (useGreatCircleDistances && (trait.length == 2)) { // Great Circle distance
                            SphericalPolarCoordinates coord1 = new SphericalPolarCoordinates(traitLow[0], traitLow[1]);
                            SphericalPolarCoordinates coord2 = new SphericalPolarCoordinates(traitUp[0], traitUp[1]);
                            double distance = coord1.distance(coord2);
                            treeDistance += distance;
                            double dc = Math.pow(distance,2)/(4*time);
                            diffusionCoefficients.add(dc);
                            waDiffusionCoefficient +=  dc*time;
                            rates.add(distance/time);
                        } else {
                            double distance = getNativeDistance(trait, parentTrait);
                            treeDistance += distance;
                            double dc = Math.pow(distance,2)/(4*time);
                            diffusionCoefficients.add(dc);
                            waDiffusionCoefficient += dc*time;
                            rates.add(distance/time);
                        }
                    }
                }
            }
        }

        if (!diffusionCoefficient){
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
        }  else {
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


    enum Mode {
        AVERAGE,
        WEIGHTED_AVERAGE,
        MEDIAN,
        COEFFICIENT_OF_VARIATION
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

            boolean option = xo.getAttribute(BOOLEAN_DIS_OPTION, false); // Default value is false
            Mode averageMode;
            String mode = xo.getAttribute(MODE, WEIGHTEDAVERAGE);

            if (mode.equals(AVERAGE)) {
                averageMode = Mode.AVERAGE;
            } else if (mode.equals(MEDIAN)) {
                averageMode = Mode.MEDIAN;
            } else if (mode.equals(COEFFICIENT_OF_VARIATION)) {
                averageMode = Mode.COEFFICIENT_OF_VARIATION;
            } else {
                averageMode = Mode.WEIGHTED_AVERAGE;
            }

            boolean diffCoeff = xo.getAttribute(BOOLEAN_DC_OPTION, false); // Default value is false

            final double upperHeight = xo.hasAttribute(HEIGHT_UPPER) ? xo.getDoubleAttribute(HEIGHT_UPPER) : Double.MAX_VALUE;
            final double lowerHeight = xo.hasAttribute(HEIGHT_LOWER) ? xo.getDoubleAttribute(HEIGHT_LOWER) : 0;

            List<AbstractMultivariateTraitLikelihood> traitLikelihoods = new ArrayList<AbstractMultivariateTraitLikelihood>();
            DiscretizedBranchRates branchRates = null;

            for (int i = 0; i < xo.getChildCount(); i++) {
                if (xo.getChild(i) instanceof AbstractMultivariateTraitLikelihood) {
                    traitLikelihoods.add((AbstractMultivariateTraitLikelihood) xo.getChild(i));
                }
                if (xo.getChild(i) instanceof DiscretizedBranchRates) {
                    branchRates = (DiscretizedBranchRates) xo.getChild(i);
                }
            }

            return new DiffusionRateStatistic(name, traitLikelihoods, option, averageMode, diffCoeff, upperHeight, lowerHeight, branchRates);
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
                AttributeRule.newBooleanRule(BOOLEAN_DC_OPTION, true),
                AttributeRule.newStringRule(MODE, true),
                AttributeRule.newDoubleRule(HEIGHT_UPPER, true),
                AttributeRule.newDoubleRule(HEIGHT_LOWER, true),
                new ElementRule(MultivariateTraitTree.class),
                new ElementRule(DiscretizedBranchRates.class),
                new ElementRule(AbstractMultivariateTraitLikelihood.class, 1, Integer.MAX_VALUE),
        };
    };

    private boolean useGreatCircleDistances;
    private List<AbstractMultivariateTraitLikelihood> traitLikelihoods;
    private Mode summaryMode;
    private boolean diffusionCoefficient;
    private double heightUpper;
    private double heightLower;
    private DiscretizedBranchRates branchRates;
}
