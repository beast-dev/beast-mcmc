/*
 * PrecisionMatrixGibbsOperator.java
 *
 * Copyright (C) 2002-2007 Alexei Drummond and Andrew Rambaut
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

package dr.evomodel.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.continuous.MultivariateTraitLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.geo.GeoSpatialDistribution;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.MatrixParameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.math.distributions.MultivariateDistribution;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.xml.*;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Marc Suchard
 */
public class TraitGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    public static final String GIBBS_OPERATOR = "traitGibbsOperator";
    public static final String INTERNAL_ONLY = "onlyInternalNodes";
    public static final String NODE_PRIOR = "nodePrior";
    public static final String NODE_LABEL = "taxon";

    private final TreeModel treeModel;
    private final MatrixParameter precisionMatrixParameter;
    private final MultivariateTraitLikelihood traitModel;
    private final int dim;
    private final String traitName;

    private Map<NodeRef, GeoSpatialDistribution> nodePrior;
    private boolean onlyInternalNodes = true;

    public TraitGibbsOperator(MultivariateTraitLikelihood traitModel, boolean onlyInternalNodes) {
        super();
        this.traitModel = traitModel;
        this.treeModel = traitModel.getTreeModel();
        this.precisionMatrixParameter = (MatrixParameter) traitModel.getDiffusionModel().getPrecisionParameter();
        this.traitName = traitModel.getTraitName();
        this.onlyInternalNodes = onlyInternalNodes;
        this.dim = treeModel.getMultivariateNodeTrait(treeModel.getRoot(), traitName).length;
    }

    public void setNodePrior(NodeRef node, GeoSpatialDistribution distribution) {
        if (nodePrior == null)
            nodePrior = new HashMap<NodeRef, GeoSpatialDistribution>();

        nodePrior.put(node, distribution);
    }


    public int getStepCount() {
        return 1;
    }


    public double doOperation() throws OperatorFailedException {

        double[][] precision = precisionMatrixParameter.getParameterAsMatrix();

        NodeRef node = treeModel.getRoot();
        while (node == treeModel.getRoot()) {
            if (onlyInternalNodes)
                node = treeModel.getInternalNode(MathUtils.nextInt(
                        treeModel.getInternalNodeCount()));
            else
                node = treeModel.getNode(MathUtils.nextInt(
                        treeModel.getNodeCount()));
        } // select any internal (or internal/external) node but the root.
        // TODO Should use Gibbs update on root given MVN prior

        NodeRef parent = treeModel.getParent(node);

        double[] weightedAverage = new double[dim];

        double weight = 1.0 / traitModel.getRescaledBranchLength(node);

        double[] trait = treeModel.getMultivariateNodeTrait(parent, traitName);

        for (int i = 0; i < dim; i++)
            weightedAverage[i] = trait[i] * weight;

        double weightTotal = weight;
        for (int j = 0; j < treeModel.getChildCount(node); j++) {
            NodeRef child = treeModel.getChild(node, j);
            trait = treeModel.getMultivariateNodeTrait(child, traitName);
            weight = 1.0 / traitModel.getRescaledBranchLength(child);

            for (int i = 0; i < dim; i++)
                weightedAverage[i] += trait[i] * weight;

            weightTotal += weight;
        }

        for (int i = 0; i < dim; i++) {
            weightedAverage[i] /= weightTotal;
            for (int j = i; j < dim; j++)
                precision[j][i] = precision[i][j] *= weightTotal;
        }

        boolean priorExists = nodePrior != null && nodePrior.containsKey(node);

        double[] draw;

        do {
            draw = MultivariateNormalDistribution.nextMultivariateNormalPrecision(
                    weightedAverage, precision);
        } while (priorExists &&  // There is a prior for this node
                (nodePrior.get(node)).logPdf(draw) == Double.NEGATIVE_INFINITY); // And draw is invalid under prior
        // TODO Currently only works for flat/truncated priors, make work for MVN

        treeModel.setMultivariateTrait(node, traitName, draw);

        return 0;  //To change body of implemented methods use File | Settings | File Templates.

    }

    public String getPerformanceSuggestion() {
        return null;
    }

    public String getOperatorName() {
        return GIBBS_OPERATOR;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return GIBBS_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);
            boolean onlyInternalNodes = xo.getAttribute(INTERNAL_ONLY, true);
            MultivariateTraitLikelihood traitModel = (MultivariateTraitLikelihood) xo.getChild(MultivariateTraitLikelihood.class);

            TraitGibbsOperator operator = new TraitGibbsOperator(traitModel, onlyInternalNodes);
            operator.setWeight(weight);

            // Get node priors
            for (int i = 0; i < xo.getChildCount(); i++) {
                 if (xo.getChild(i) instanceof MultivariateDistributionLikelihood) {
                    MultivariateDistribution dist = ((MultivariateDistributionLikelihood) xo.getChild(i)).getDistribution();
                    if (dist instanceof GeoSpatialDistribution) {
                        GeoSpatialDistribution prior = (GeoSpatialDistribution) dist;
                        String nodeLabel = prior.getLabel();
                        TreeModel treeModel = traitModel.getTreeModel();

                        // Get taxon node from tree
                        int index = treeModel.getTaxonIndex(nodeLabel);
                        if (index == -1) {
                            throw new XMLParseException("taxon '" + nodeLabel + "' not found for geoSpatialDistribution element in traitGibbsOperator element");
                        }
                        NodeRef node = treeModel.getExternalNode(index);
                        operator.setNodePrior(node, prior);
                        System.err.println("Adding truncated prior for "+node);
                    }
                }
            }
            return operator;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a multivariate Gibbs operator on traits for possible all nodes.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(WEIGHT),
                AttributeRule.newBooleanRule(INTERNAL_ONLY, true),
                new ElementRule(MultivariateTraitLikelihood.class),
//                new ElementRule(NODE_PRIOR, new XMLSyntaxRule[] {
//                        AttributeRule.newStringRule(NODE_LABEL),
//                        new ElementRule(MultivariateDistributionLikelihood.class),
//                }),
                new ElementRule(MultivariateDistributionLikelihood.class, 0, Integer.MAX_VALUE),
        };

    };

}