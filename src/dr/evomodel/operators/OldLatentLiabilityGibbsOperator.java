/*
 * OldLatentLiabilityGibbsOperator.java
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

package dr.evomodel.operators;

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evomodel.continuous.AbstractMultivariateTraitLikelihood;
import dr.evomodel.continuous.BinaryLatentLiabilityLikelihood;
import dr.evomodel.continuous.IntegratedMultivariateTraitLikelihood;
import dr.inference.model.MatrixParameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.util.Citable;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Marc A. Suchard
 * @author Joe Felsenstein
 */
public class OldLatentLiabilityGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    public static final String GIBBS_OPERATOR = "oldLatentLiabilityGibbsOperator";
//    public static final String INTERNAL_ONLY = "onlyInternalNodes";
//    public static final String TIP_WITH_PRIORS_ONLY = "onlyTipsWithPriors";
//    public static final String NODE_PRIOR = "nodePrior";
//    public static final String NODE_LABEL = "taxon";
//    public static final String ROOT_PRIOR = "rootPrior";

    private final MutableTreeModel treeModel;
    private final MatrixParameter precisionMatrixParameter;
    private final IntegratedMultivariateTraitLikelihood traitModel;
    private final BinaryLatentLiabilityLikelihood liabilityLikelihood;

    private final int dim;
    private final String traitName;

//    private Map<Taxon, GeoSpatialDistribution> nodeGeoSpatialPrior;
//    private Map<Taxon, MultivariateNormalDistribution> nodeMVNPrior;
//    private GeoSpatialCollectionModel parameterPrior = null;

//    private boolean onlyInternalNodes = true;
//    private boolean onlyTipsWithPriors = true;
//    private boolean sampleRoot = false;
//    private double[] rootPriorMean;
//    private double[][] rootPriorPrecision;

    private final int maxTries = 10000;

    public OldLatentLiabilityGibbsOperator(IntegratedMultivariateTraitLikelihood traitModel,
                                           BinaryLatentLiabilityLikelihood liabilityLikelihood) {
        super();
        this.traitModel = traitModel;
        this.liabilityLikelihood = liabilityLikelihood;

        this.treeModel = traitModel.getTreeModel();
        this.precisionMatrixParameter = (MatrixParameter) traitModel.getDiffusionModel().getPrecisionParameter();
        this.traitName = traitModel.getTraitName();
        this.dim = traitModel.getDimTrait();

        StringBuilder sb = new StringBuilder();
        sb.append("Using a latent trait Gibbs operator.  Please cite:");
        sb.append(Citable.Utils.getCitationString(liabilityLikelihood));
        Logger.getLogger("dr.evomodel.continuous").info(sb.toString());
    }

//    public void setRootPrior(MultivariateNormalDistribution rootPrior) {
//        rootPriorMean = rootPrior.getMean();
//        rootPriorPrecision = rootPrior.getScaleMatrix();
//        sampleRoot = true;
//    }

//    public void setTaxonPrior(Taxon taxon, MultivariateDistribution distribution) {
//
//        if (distribution instanceof GeoSpatialDistribution) {
//            if (nodeGeoSpatialPrior == null) {
//                nodeGeoSpatialPrior = new HashMap<Taxon, GeoSpatialDistribution>();
//            }
//            nodeGeoSpatialPrior.put(taxon, (GeoSpatialDistribution)distribution);
//
//        } else if (distribution instanceof MultivariateNormalDistribution) {
//            if (nodeMVNPrior == null) {
//                nodeMVNPrior = new HashMap<Taxon, MultivariateNormalDistribution>();
//            }
//            nodeMVNPrior.put(taxon, (MultivariateNormalDistribution)distribution);
//        } else {
//            throw new RuntimeException("Only flat/truncated geospatial and multivariate normal distributions allowed");
//        }
//    }
//
//    public void setParameterPrior(GeoSpatialCollectionModel distribution) {
//        parameterPrior = distribution;
//    }

    public int getStepCount() {
        return 1;
    }

//    private boolean nodeGeoSpatialPriorExists(NodeRef node) {
//        return nodeGeoSpatialPrior != null && nodeGeoSpatialPrior.containsKey(treeModel.getNodeTaxon(node));
//    }
//
//    private boolean nodeMVNPriorExists(NodeRef node) {
//        return nodeMVNPrior != null && nodeMVNPrior.containsKey(treeModel.getNodeTaxon(node));
//    }

    public double doOperation() {

        traitModel.redrawAncestralStates();

        NodeRef node = treeModel.getNode(MathUtils.nextInt(treeModel.getExternalNodeCount()));
        int tip = node.getNumber();

        // Draw truncated MVN using rejection sampling
        do {
            // Nothing
        } while (!liabilityLikelihood.validTraitForTip(tip));

//        NodeRef node = null;
//        final NodeRef root = treeModel.getRoot();

//        while (node == null) {
//            if (onlyInternalNodes)
//                node = treeModel.getInternalNode(MathUtils.nextInt(
//                        treeModel.getInternalNodeCount()));
//            else {
//                node = treeModel.getNode(MathUtils.nextInt(
//                        treeModel.getNodeCount()));
//                if (onlyTipsWithPriors &&
//                    (treeModel.getChildCount(node) == 0) && // Is a tip
//                    !nodeGeoSpatialPriorExists(node)) { // Does not have a prior
//                    node = null;
//                }
//            }
//            if (!sampleRoot && node == root)
//                node = null;
//        } // select any internal (or internal/external) node
//
//        final double[] initialValue = treeModel.getMultivariateNodeTrait(node,traitName);
//
//        MeanPrecision mp;
//
//        if (node != root)
//            mp = operateNotRoot(node);
//        else
//            mp = operateRoot(node);
//
//        final Taxon taxon = treeModel.getNodeTaxon(node);
//
//        final boolean nodePriorExists = nodeGeoSpatialPriorExists(node);
//
//        int count = 0;
//
//        final boolean parameterPriorExists = parameterPrior != null;
//
//        double[] draw;
//
//        do {
//            do {
//                if (count > maxTries)  {
//                    treeModel.setMultivariateTrait(node,traitName,initialValue);
//                    throw new OperatorFailedException("Truncated Gibbs is stuck!");
//                }
//
//                draw = MultivariateNormalDistribution.nextMultivariateNormalPrecision(
//                        mp.mean, mp.precision);
//                count++;
//
//            } while (nodePriorExists &&  // There is a prior for this node
//                    (nodeGeoSpatialPrior.get(taxon)).logPdf(draw) == Double.NEGATIVE_INFINITY); // And draw is invalid under prior
//
//            treeModel.setMultivariateTrait(node, traitName, draw);
//
//        } while (parameterPriorExists &&
//                (parameterPrior.getLogLikelihood() == Double.NEGATIVE_INFINITY));

        return 0;
    }

//    private MeanPrecision operateNotRoot(NodeRef node) {
//
//        double[][] precision = precisionMatrixParameter.getParameterAsMatrix();
//
//        NodeRef parent = treeModel.getParent(node);
//
//        double[] mean = new double[dim];
//
//        double weight = 1.0 / traitModel.getRescaledBranchLength(node);
//
//        double[] trait = treeModel.getMultivariateNodeTrait(parent, traitName);
//
//        for (int i = 0; i < dim; i++)
//            mean[i] = trait[i] * weight;
//
//        double weightTotal = weight;
//        for (int j = 0; j < treeModel.getChildCount(node); j++) {
//            NodeRef child = treeModel.getChild(node, j);
//            trait = treeModel.getMultivariateNodeTrait(child, traitName);
//            weight = 1.0 / traitModel.getRescaledBranchLength(child);
//
//            for (int i = 0; i < dim; i++)
//                mean[i] += trait[i] * weight;
//
//            weightTotal += weight;
//        }
//
//        for (int i = 0; i < dim; i++) {
//            mean[i] /= weightTotal;
//            for (int j = i; j < dim; j++)
//                precision[j][i] = precision[i][j] *= weightTotal;
//        }
//
////        if (nodeMVNPriorExists(node)) {
////            throw new RuntimeException("Still trying to implement multivariate normal taxon priors");
////        }
//
//        return new MeanPrecision(mean,precision);
//    }

    class MeanPrecision {
        final double[] mean;
        final double[][] precision;

        MeanPrecision(double[] mean, double[][] precision) {
            this.mean = mean;
            this.precision = precision;
        }
    }

//    private MeanPrecision operateRoot(NodeRef node) {
//
//        double[] trait;
//        double weightTotal = 0.0;
//
//        double[] weightedAverage = new double[dim];
//
//        double[][] precision = precisionMatrixParameter.getParameterAsMatrix();
//
//        for (int k = 0; k < treeModel.getChildCount(node); k++) {
//            NodeRef child = treeModel.getChild(node, k);
//            trait = treeModel.getMultivariateNodeTrait(child, traitName);
//            final double weight = 1.0 / traitModel.getRescaledBranchLength(child);
//
//            for (int i = 0; i < dim; i++) {
//                for (int j=0; j<dim; j++)
//                    weightedAverage[i] += precision[i][j] * weight * trait[j];
//            }
//
//            weightTotal += weight;
//        }
//
//        for (int i=0; i<dim; i++) {
//            for (int j=0; j<dim; j++) {
//                weightedAverage[i] += rootPriorPrecision[i][j] * rootPriorMean[j];
//                precision[i][j]  = precision[i][j] * weightTotal + rootPriorPrecision[i][j];
//            }
//        }
//
//        double[][] variance = new SymmetricMatrix(precision).inverse().toComponents();
//
//        trait = new double[dim];
//        for (int i=0; i<dim; i++) {
//            for (int j=0; j<dim; j++)
//                trait[i] += variance[i][j] * weightedAverage[j];
//        }
//
//        return new MeanPrecision(trait,precision);
//    }

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
//            boolean onlyInternalNodes = xo.getAttribute(INTERNAL_ONLY, true);
//            boolean onlyTipsWithPriors = xo.getAttribute(TIP_WITH_PRIORS_ONLY, true);
//            boolean onlyInternalNodes = true;
//            boolean onlyTipsWithPriors = true;

            IntegratedMultivariateTraitLikelihood traitModel = (IntegratedMultivariateTraitLikelihood)
                    xo.getChild(AbstractMultivariateTraitLikelihood.class);

            BinaryLatentLiabilityLikelihood liabilityLikelihood = (BinaryLatentLiabilityLikelihood)
                    xo.getChild(BinaryLatentLiabilityLikelihood.class);

            OldLatentLiabilityGibbsOperator operator = new OldLatentLiabilityGibbsOperator(traitModel, liabilityLikelihood);
            operator.setWeight(weight);

            return operator;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a multivariate Gibbs operator on traits for tip nodes under a latent liability model.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(WEIGHT),
                new ElementRule(IntegratedMultivariateTraitLikelihood.class),
                new ElementRule(BinaryLatentLiabilityLikelihood.class),
        };

    };

}