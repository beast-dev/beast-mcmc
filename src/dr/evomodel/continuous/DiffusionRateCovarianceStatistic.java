/*
 * DiffusionRateCovarianceStatistic.java
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
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeStatistic;
import dr.geo.math.SphericalPolarCoordinates;
import dr.inference.model.Statistic;
import dr.stats.DiscreteStatistics;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc Suchard
 * @author Philippe Lemey
 * @author Andrew Rambaut
 */
public class DiffusionRateCovarianceStatistic extends Statistic.Abstract {

    public static final String DIFFUSION_RATE_COVARIANCE_STATISTIC = "diffusionRateCovarianceStatistic";
    public static final String TREE_DISPERSION_COVARIANCE_STATISTIC = "treeDispersionCovarianceStatistic";
    public static final String BOOLEAN_DIS_OPTION = "greatCircleDistance";
    public static final String BOOLEAN_DC_OPTION = "diffusionCoefficient";

    public DiffusionRateCovarianceStatistic(String name, TreeModel tree, List<AbstractMultivariateTraitLikelihood> traitLikelihoods,
                                            boolean option, boolean diffusionCoefficient) {
        super(name);
        this.traitLikelihoods = traitLikelihoods;
        useGreatCircleDistances = option;
        this.diffusionCoefficient = diffusionCoefficient;
        int n = tree.getExternalNodeCount();
        childRate = new double[2 * n - 4];
        parentRate = new double[childRate.length];
    }

    public int getDimension() {
        return 1;
    }

    /**
     * @return whatever Philippe wants
     */
    public double getStatisticValue(int dim) {

        String traitName = traitLikelihoods.get(0).getTraitName();

        for (AbstractMultivariateTraitLikelihood traitLikelihood : traitLikelihoods) {
            MutableTreeModel tree = traitLikelihood.getTreeModel();

            int counter = 0;
            int index = 0;

            for (int i = 0; i < tree.getNodeCount(); i++) {

                NodeRef child = tree.getNode(i);
                NodeRef parent = tree.getParent(child);
                if (parent != null & !tree.isRoot(parent)) {
                    double[] childTrait = traitLikelihood.getTraitForNode(tree, child, traitName);
                    double[] parentTrait = traitLikelihood.getTraitForNode(tree, parent, traitName);
                    double childTime = tree.getBranchLength(child);
                    double parentTime = tree.getBranchLength(parent);

                    NodeRef grandParent = tree.getParent(parent);
                    double[] grandParentTrait = traitLikelihood.getTraitForNode(tree, grandParent, traitName);

                    if (useGreatCircleDistances && (childTrait.length == 2)) { // Great Circle distance
                        SphericalPolarCoordinates coord1 = new SphericalPolarCoordinates(childTrait[0], childTrait[1]);
                        SphericalPolarCoordinates coord2 = new SphericalPolarCoordinates(parentTrait[0], parentTrait[1]);
                        double childDistance = coord1.distance(coord2);
                        SphericalPolarCoordinates coord3 = new SphericalPolarCoordinates(grandParentTrait[0], grandParentTrait[1]);
                        double parentDistance = coord2.distance(coord3);

                        if (!diffusionCoefficient){
                            childRate[index] = childDistance/childTime;
                            parentRate[index] = parentDistance/parentTime;
                        } else {
                            childRate[index] = Math.pow(childDistance,2)/(4*childTime);
                            parentRate[index] = Math.pow(parentDistance,2)/(4*parentTime);
                        }

                    } else {
                        double childDistance = getNativeDistance(childTrait, parentTrait);
                        double parentDistance = getNativeDistance(parentTrait, grandParentTrait);

                        if (!diffusionCoefficient){
                            childRate[index] = childDistance/childTime;
                            parentRate[index] = parentDistance/parentTime;
                        } else {
                            childRate[index] = Math.pow(childDistance,2)/(4*childTime);
                            parentRate[index] = Math.pow(parentDistance,2)/(4*parentTime);
                        }

                    }
                    index += 1;
                }

            }
        }

        return DiscreteStatistics.covariance(childRate, parentRate);
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

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return DIFFUSION_RATE_COVARIANCE_STATISTIC;
        }

        @Override
        public String[] getParserNames() {
            return new String[]{getParserName(), TREE_DISPERSION_COVARIANCE_STATISTIC};
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String name = xo.getAttribute(NAME, xo.getId());
            TreeModel tree = (TreeModel) xo.getChild(Tree.class);

            boolean option = xo.getAttribute(BOOLEAN_DIS_OPTION, false); // Default value is false

            boolean diffCoeff = xo.getAttribute(BOOLEAN_DC_OPTION, false); // Default value is false

            List<AbstractMultivariateTraitLikelihood> traitLikelihoods = new ArrayList<AbstractMultivariateTraitLikelihood>();

            for (int i = 0; i < xo.getChildCount(); i++) {
                if (xo.getChild(i) instanceof AbstractMultivariateTraitLikelihood) {
                    traitLikelihoods.add((AbstractMultivariateTraitLikelihood) xo.getChild(i));
                }
            }

            return new DiffusionRateCovarianceStatistic(name, tree, traitLikelihoods, option, diffCoeff);
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
                new ElementRule(MutableTreeModel.class),
                new ElementRule(AbstractMultivariateTraitLikelihood.class, 1, Integer.MAX_VALUE),
        };
    };

    private boolean useGreatCircleDistances;
    private List<AbstractMultivariateTraitLikelihood> traitLikelihoods;
    private double[] childRate = null;
    private double[] parentRate = null;
    private boolean diffusionCoefficient;
}
