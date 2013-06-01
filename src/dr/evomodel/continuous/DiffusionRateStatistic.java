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
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeStatistic;
import dr.geo.math.SphericalPolarCoordinates;
import dr.inference.model.Statistic;
import dr.stats.DiscreteStatistics;
import dr.util.HeapSort;
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

    public DiffusionRateStatistic(String name, TreeModel tree, List<AbstractMultivariateTraitLikelihood> traitLikelihoods,
                                  boolean option, Mode mode, boolean diffusionCoefficient) {
        super(name);
        this.traitLikelihoods = traitLikelihoods;
        this.useGreatCircleDistances = option;
        summaryMode =  mode;
        this.diffusionCoefficient = diffusionCoefficient;
    }

    public int getDimension() {
        return 1;
    }

    public double getStatisticValue(int dim) {

        String traitName = traitLikelihoods.get(0).getTraitName();
        double treelength = 0;
        double treeDistance = 0;

        double[] rates =  null;
        double[] diffusionCoefficients =  null;
        double waDiffusionCoefficient =  0;

        for (AbstractMultivariateTraitLikelihood traitLikelihood : traitLikelihoods) {
            MultivariateTraitTree tree = traitLikelihood.getTreeModel();

            rates = new double[(tree.getNodeCount() -1)];
            diffusionCoefficients = new double[(tree.getNodeCount() -1)];

            int counter = 0;

            for (int i = 0; i < tree.getNodeCount(); i++) {
                NodeRef node = tree.getNode(i);
                double[] trait = traitLikelihood.getTraitForNode(tree, node, traitName);

                if (node != tree.getRoot()) {

                    double[] parentTrait = traitLikelihood.getTraitForNode(tree, tree.getParent(node), traitName);
                    double time = tree.getBranchLength(node);
                    treelength += time;

                    if (useGreatCircleDistances && (trait.length == 2)) { // Great Circle distance
                        SphericalPolarCoordinates coord1 = new SphericalPolarCoordinates(trait[0], trait[1]);
                        SphericalPolarCoordinates coord2 = new SphericalPolarCoordinates(parentTrait[0], parentTrait[1]);
                        double distance = coord1.distance(coord2);
                        treeDistance += distance;
                        diffusionCoefficients[counter] = Math.pow(distance,2)/(4*time);
                        waDiffusionCoefficient +=  diffusionCoefficients[counter]*time;
                        rates[counter] = distance/time;
                    } else {
                        double distance = getNativeDistance(trait, parentTrait);
                        treeDistance += distance;
                        diffusionCoefficients[counter] = Math.pow(distance,2)/(4*time);
                        waDiffusionCoefficient += diffusionCoefficients[counter]*time;
                        rates[counter] = distance/time;
                    }
                    counter ++;
                }
            }
        }
        if (!diffusionCoefficient){
            if (summaryMode == Mode.AVERAGE) {
                return DiscreteStatistics.mean(rates);
            } else if (summaryMode == Mode.MEDIAN) {
                return DiscreteStatistics.median(rates);
            } else if (summaryMode == Mode.COEFFICIENT_OF_VARIATION) {
                // don't compute mean twice
                final double mean = DiscreteStatistics.mean(rates);
                return Math.sqrt(DiscreteStatistics.variance(rates, mean)) / mean;
            } else {
                return treeDistance / treelength;
            }
        }  else {
            if (summaryMode == Mode.AVERAGE) {
                return DiscreteStatistics.mean(diffusionCoefficients);
            } else if (summaryMode == Mode.MEDIAN) {
                return DiscreteStatistics.median(diffusionCoefficients);
            } else if (summaryMode == Mode.COEFFICIENT_OF_VARIATION) {
                // don't compute mean twice
                final double mean = DiscreteStatistics.mean(diffusionCoefficients);
                return Math.sqrt(DiscreteStatistics.variance(diffusionCoefficients, mean)) / mean;
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
            TreeModel tree = (TreeModel) xo.getChild(Tree.class);

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

            List<AbstractMultivariateTraitLikelihood> traitLikelihoods = new ArrayList<AbstractMultivariateTraitLikelihood>();

            for (int i = 0; i < xo.getChildCount(); i++) {
                if (xo.getChild(i) instanceof AbstractMultivariateTraitLikelihood) {
                    traitLikelihoods.add((AbstractMultivariateTraitLikelihood) xo.getChild(i));
                }
            }

            return new DiffusionRateStatistic(name, tree, traitLikelihoods, option, averageMode, diffCoeff);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A statistic that returns the average of the branch rates";
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
                new ElementRule(MultivariateTraitTree.class),
                new ElementRule(AbstractMultivariateTraitLikelihood.class, 1, Integer.MAX_VALUE),
        };
    };

    private boolean useGreatCircleDistances;
    private List<AbstractMultivariateTraitLikelihood> traitLikelihoods;
    private Mode summaryMode;
    private boolean diffusionCoefficient;
}
