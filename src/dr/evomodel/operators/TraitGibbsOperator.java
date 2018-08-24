/*
 * TraitGibbsOperator.java
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
import dr.evolution.util.Taxon;
import dr.evomodel.continuous.AbstractMultivariateTraitLikelihood;
import dr.evomodel.continuous.SampledMultivariateTraitLikelihood;
import dr.geo.GeoSpatialCollectionModel;
import dr.geo.GeoSpatialDistribution;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.MatrixParameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.math.distributions.MultivariateDistribution;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.xml.*;

import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;

/**
 * @author Marc Suchard
 */
public class TraitGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    public static final String GIBBS_OPERATOR = "traitGibbsOperator";
    public static final String INTERNAL_ONLY = "onlyInternalNodes";
    public static final String TIP_WITH_PRIORS_ONLY = "onlyTipsWithPriors";
    public static final String NODE_PRIOR = "nodePrior";
    public static final String NODE_LABEL = "taxon";
    public static final String ROOT_PRIOR = "rootPrior";

    private final MutableTreeModel treeModel;
    private final MatrixParameter precisionMatrixParameter;
    private final SampledMultivariateTraitLikelihood traitModel;
    private final int dim;
    private final String traitName;

    private Map<Taxon, GeoSpatialDistribution> nodeGeoSpatialPrior;
    private Map<Taxon, MultivariateNormalDistribution> nodeMVNPrior;
    private GeoSpatialCollectionModel parameterPrior = null;

    private boolean onlyInternalNodes = true;
    private boolean onlyTipsWithPriors = true;
    private boolean sampleRoot = false;
    private double[] rootPriorMean;
    private double[][] rootPriorPrecision;
    private final int maxTries = 10000;

    public TraitGibbsOperator(SampledMultivariateTraitLikelihood traitModel, boolean onlyInternalNodes,
                              boolean onlyTipsWithPriors) {
        super();
        this.traitModel = traitModel;
        this.treeModel = traitModel.getTreeModel();
        this.precisionMatrixParameter = (MatrixParameter) traitModel.getDiffusionModel().getPrecisionParameter();
        this.traitName = traitModel.getTraitName();
        this.onlyInternalNodes = onlyInternalNodes;
        this.onlyTipsWithPriors = onlyTipsWithPriors;
        this.dim = treeModel.getMultivariateNodeTrait(treeModel.getRoot(), traitName).length;
        Logger.getLogger("dr.evomodel").info("Using *NEW* trait Gibbs operator");
    }

    public void setRootPrior(MultivariateNormalDistribution rootPrior) {
        rootPriorMean = rootPrior.getMean();
        rootPriorPrecision = rootPrior.getScaleMatrix();
        sampleRoot = true;
    }

    public void setTaxonPrior(Taxon taxon, MultivariateDistribution distribution) {

        if (distribution instanceof GeoSpatialDistribution) {
            if (nodeGeoSpatialPrior == null) {
                nodeGeoSpatialPrior = new HashMap<Taxon, GeoSpatialDistribution>();
            }
            nodeGeoSpatialPrior.put(taxon, (GeoSpatialDistribution) distribution);

        } else if (distribution instanceof MultivariateNormalDistribution) {
            if (nodeMVNPrior == null) {
                nodeMVNPrior = new HashMap<Taxon, MultivariateNormalDistribution>();
            }
            nodeMVNPrior.put(taxon, (MultivariateNormalDistribution) distribution);
        } else {
            throw new RuntimeException("Only flat/truncated geospatial and multivariate normal distributions allowed");
        }
    }

    public void setParameterPrior(GeoSpatialCollectionModel distribution) {
        parameterPrior = distribution;
    }

    public int getStepCount() {
        return 1;
    }

    private boolean nodeGeoSpatialPriorExists(NodeRef node) {
        return nodeGeoSpatialPrior != null && nodeGeoSpatialPrior.containsKey(treeModel.getNodeTaxon(node));
    }

    private boolean nodeMVNPriorExists(NodeRef node) {
        return nodeMVNPrior != null && nodeMVNPrior.containsKey(treeModel.getNodeTaxon(node));
    }

    public double doOperation() {

        NodeRef node = null;
        final NodeRef root = treeModel.getRoot();

        while (node == null) {
            if (onlyInternalNodes)
                node = treeModel.getInternalNode(MathUtils.nextInt(
                        treeModel.getInternalNodeCount()));
            else {
                node = treeModel.getNode(MathUtils.nextInt(
                        treeModel.getNodeCount()));
                if (onlyTipsWithPriors &&
                        (treeModel.getChildCount(node) == 0) && // Is a tip
                        !nodeGeoSpatialPriorExists(node)) { // Does not have a prior
                    node = null;
                }
            }
            if (!sampleRoot && node == root)
                node = null;
        } // select any internal (or internal/external) node

        final double[] initialValue = treeModel.getMultivariateNodeTrait(node, traitName);

        MeanPrecision mp;

        if (node != root)
            mp = operateNotRoot(node);
        else
            mp = operateRoot(node);

        final Taxon taxon = treeModel.getNodeTaxon(node);

//        final boolean nodePriorExists = nodeGeoSpatialPrior != null && nodeGeoSpatialPrior.containsKey(taxon);
        final boolean nodePriorExists = nodeGeoSpatialPriorExists(node);

//        if (!onlyInternalNodes) {
//            final boolean isTip = (treeModel.getChildCount(node) == 0);
//            if (!nodePriorExists && isTip)
//                System.err.println("Warning: sampling taxon '"+treeModel.getNodeTaxon(node).getId()
//                        +"' tip trait without a prior!!!");
//        }

        int count = 0;

        final boolean parameterPriorExists = parameterPrior != null;

        double[] draw;

        do {
            do {
                if (count > maxTries) {
                    treeModel.setMultivariateTrait(node, traitName, initialValue);  // TODO Add to MTT interface
                    throw new RuntimeException("Truncated Gibbs is stuck!");
                }

                draw = MultivariateNormalDistribution.nextMultivariateNormalPrecision(
                        mp.mean, mp.precision);
                count++;

            } while (nodePriorExists &&  // There is a prior for this node
                    (nodeGeoSpatialPrior.get(taxon)).logPdf(draw) == Double.NEGATIVE_INFINITY); // And draw is invalid under prior
            // TODO Currently only works for flat/truncated priors, make work for MVN

            treeModel.setMultivariateTrait(node, traitName, draw);

        } while (parameterPriorExists &&
                (parameterPrior.getLogLikelihood() == Double.NEGATIVE_INFINITY));

        return 0;
    }

    private MeanPrecision operateNotRoot(NodeRef node) {

        double[][] precision = precisionMatrixParameter.getParameterAsMatrix();

        NodeRef parent = treeModel.getParent(node);

        double[] mean = new double[dim];

        double weight = 1.0 / traitModel.getRescaledBranchLengthForPrecision(node);

        double[] trait = treeModel.getMultivariateNodeTrait(parent, traitName);

        for (int i = 0; i < dim; i++)
            mean[i] = trait[i] * weight;

        double weightTotal = weight;
        for (int j = 0; j < treeModel.getChildCount(node); j++) {
            NodeRef child = treeModel.getChild(node, j);
            trait = treeModel.getMultivariateNodeTrait(child, traitName);
            weight = 1.0 / traitModel.getRescaledBranchLengthForPrecision(child);

            for (int i = 0; i < dim; i++)
                mean[i] += trait[i] * weight;

            weightTotal += weight;
        }

        for (int i = 0; i < dim; i++) {
            mean[i] /= weightTotal;
            for (int j = i; j < dim; j++)
                precision[j][i] = precision[i][j] *= weightTotal;
        }

        if (nodeMVNPriorExists(node)) {
            throw new RuntimeException("Still trying to implement multivariate normal taxon priors");
        }

        return new MeanPrecision(mean, precision);
    }

    class MeanPrecision {
        final double[] mean;
        final double[][] precision;

        MeanPrecision(double[] mean, double[][] precision) {
            this.mean = mean;
            this.precision = precision;
        }
    }

    private MeanPrecision operateRoot(NodeRef node) {

        double[] trait;
        double weightTotal = 0.0;

        double[] weightedAverage = new double[dim];

        double[][] precision = precisionMatrixParameter.getParameterAsMatrix();

        for (int k = 0; k < treeModel.getChildCount(node); k++) {
            NodeRef child = treeModel.getChild(node, k);
            trait = treeModel.getMultivariateNodeTrait(child, traitName);
            final double weight = 1.0 / traitModel.getRescaledBranchLengthForPrecision(child);

            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++)
                    weightedAverage[i] += precision[i][j] * weight * trait[j];
            }

            weightTotal += weight;
        }

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                weightedAverage[i] += rootPriorPrecision[i][j] * rootPriorMean[j];
                precision[i][j] = precision[i][j] * weightTotal + rootPriorPrecision[i][j];
            }
        }

        double[][] variance = new SymmetricMatrix(precision).inverse().toComponents();

        trait = new double[dim];
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++)
                trait[i] += variance[i][j] * weightedAverage[j];
        }

        return new MeanPrecision(trait, precision);
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

        private final String[] names = {GIBBS_OPERATOR, "internalTraitGibbsOperator"};

        public String[] getParserNames() {
            return names;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);
            boolean onlyInternalNodes = xo.getAttribute(INTERNAL_ONLY, true);
            boolean onlyTipsWithPriors = xo.getAttribute(TIP_WITH_PRIORS_ONLY, true);
            SampledMultivariateTraitLikelihood traitModel = (SampledMultivariateTraitLikelihood) xo.getChild(AbstractMultivariateTraitLikelihood.class);

            TraitGibbsOperator operator = new TraitGibbsOperator(traitModel, onlyInternalNodes, onlyTipsWithPriors);
            operator.setWeight(weight);

            // Get root prior

            XMLObject cxo = xo.getChild(ROOT_PRIOR);
            if (cxo != null) {

                MultivariateDistributionLikelihood rootPrior = (MultivariateDistributionLikelihood) cxo.getChild(MultivariateDistributionLikelihood.class);
                if (!(rootPrior.getDistribution() instanceof MultivariateDistribution))
                    throw new XMLParseException("Only multivariate normal priors allowed for Gibbs sampling the root trait");
                operator.setRootPrior((MultivariateNormalDistribution) rootPrior.getDistribution());
            }


            // Get node priors
            for (int i = 0; i < xo.getChildCount(); i++) {
                if (xo.getChild(i) instanceof MultivariateDistributionLikelihood) {
                    MultivariateDistribution dist = ((MultivariateDistributionLikelihood) xo.getChild(i)).getDistribution();
                    if (dist instanceof GeoSpatialDistribution) {
                        GeoSpatialDistribution prior = (GeoSpatialDistribution) dist;
                        String nodeLabel = prior.getLabel();
                        Taxon taxon = getTaxon(traitModel.getTreeModel(), nodeLabel);
                        operator.setTaxonPrior(taxon, prior);
                        System.err.println("Adding truncated prior for taxon '" + taxon + "'");
                    }
                }
            }

            GeoSpatialCollectionModel collectionModel = (GeoSpatialCollectionModel) xo.getChild(GeoSpatialCollectionModel.class);
            if (collectionModel != null) {
                operator.setParameterPrior(collectionModel);
                System.err.println("Adding truncated prior '" + collectionModel.getId() +
                        "' for parameter '" + collectionModel.getParameter().getId() + "'");
            }

            return operator;
        }

        private Taxon getTaxon(MutableTreeModel treeModel, String taxonLabel) throws XMLParseException {
            // Get taxon node from tree
            int index = treeModel.getTaxonIndex(taxonLabel);
            if (index == -1) {
                throw new XMLParseException("Taxon '" + taxonLabel + "' not found for geoSpatialDistribution element in traitGibbsOperator element");
            }
            return treeModel.getTaxon(index);
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
                AttributeRule.newBooleanRule(TIP_WITH_PRIORS_ONLY, true),
                new ElementRule(SampledMultivariateTraitLikelihood.class),
//                new ElementRule(NODE_PRIOR, new XMLSyntaxRule[] {
//                        AttributeRule.newStringRule(NODE_LABEL),
//                        new ElementRule(MultivariateDistributionLikelihood.class),
//                }),
                new ElementRule(MultivariateDistributionLikelihood.class, 0, Integer.MAX_VALUE),
                new ElementRule(ROOT_PRIOR,
                        new XMLSyntaxRule[]{
                                new ElementRule(MultivariateDistributionLikelihood.class)
                        }, true),
                new ElementRule(GeoSpatialCollectionModel.class, true),
        };

    };

}