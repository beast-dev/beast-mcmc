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
import dr.evolution.util.Taxon;
import dr.evomodel.continuous.AbstractMultivariateTraitLikelihood;
import dr.evomodel.continuous.SampledMultivariateTraitLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.geo.GeoSpatialDistribution;
import dr.geo.GeoSpatialCollectionModel;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.MatrixParameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.math.distributions.MultivariateDistribution;
import dr.math.distributions.MultivariateNormalDistribution;
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
    public static final String NODE_PRIOR = "nodePrior";
    public static final String NODE_LABEL = "taxon";
    public static final String ROOT_PRIOR = "rootPrior";

    private final TreeModel treeModel;
    private final MatrixParameter precisionMatrixParameter;
    private final SampledMultivariateTraitLikelihood traitModel;
    private final int dim;
    private final String traitName;

    private Map<Taxon, GeoSpatialDistribution> nodePrior;
    private GeoSpatialCollectionModel parameterPrior = null;

    private boolean onlyInternalNodes = true;
    private boolean sampleRoot = false;
    private double[] rootPriorMean;
    private double[][] rootPriorPrecision;
    private final int maxTries = 10000;

    public TraitGibbsOperator(SampledMultivariateTraitLikelihood traitModel, boolean onlyInternalNodes) {
        super();
        this.traitModel = traitModel;
        this.treeModel = traitModel.getTreeModel();
        this.precisionMatrixParameter = (MatrixParameter) traitModel.getDiffusionModel().getPrecisionParameter();
        this.traitName = traitModel.getTraitName();
        this.onlyInternalNodes = onlyInternalNodes;
        this.dim = treeModel.getMultivariateNodeTrait(treeModel.getRoot(), traitName).length;
        Logger.getLogger("dr.evomodel").info("Using *NEW* trait Gibbs operator");
    }

    public void setRootPrior(MultivariateNormalDistribution rootPrior) {
        rootPriorMean = rootPrior.getMean();
        rootPriorPrecision = rootPrior.getScaleMatrix();
        sampleRoot = true;
    }

    public void setTaxonPrior(Taxon taxon, GeoSpatialDistribution distribution) {
        if (nodePrior == null)
            nodePrior = new HashMap<Taxon, GeoSpatialDistribution>();

        nodePrior.put(taxon, distribution);
    }

    public void setParameterPrior(GeoSpatialCollectionModel distribution) {
        parameterPrior = distribution;
    }

    public int getStepCount() {
        return 1;
    }


    public double doOperation() throws OperatorFailedException {

        NodeRef node = null;
        final NodeRef root = treeModel.getRoot();

        while (node == null) {
            if (onlyInternalNodes)
                node = treeModel.getInternalNode(MathUtils.nextInt(
                        treeModel.getInternalNodeCount()));
            else
                node = treeModel.getNode(MathUtils.nextInt(
                        treeModel.getNodeCount()));
            if (!sampleRoot && node == root)
                node = null;
        } // select any internal (or internal/external) node

        final double[] initialValue = treeModel.getMultivariateNodeTrait(node,traitName);

        MeanPrecision mp;

        if (node != root)
            mp = operateNotRoot(node);
        else
            mp = operateRoot(node);

        final Taxon taxon = treeModel.getNodeTaxon(node);

        final boolean nodePriorExists = nodePrior != null && nodePrior.containsKey(taxon);

        if (!onlyInternalNodes) {
            final boolean isTip = (treeModel.getChildCount(node) == 0);
            if (!nodePriorExists && isTip)
                System.err.println("Warning: sampling taxon '"+treeModel.getNodeTaxon(node).getId()
                        +"' tip trait without a prior!!!");
        }

        int count = 0;

        final boolean parameterPriorExists = parameterPrior != null;
       
        double[] draw;

        do {
            do {
                if (count > maxTries)  {
                    treeModel.setMultivariateTrait(node,traitName,initialValue);
                    throw new OperatorFailedException("Truncated Gibbs is stuck!");
                }

                draw = MultivariateNormalDistribution.nextMultivariateNormalPrecision(
                        mp.mean, mp.precision);
                count++;

            } while (nodePriorExists &&  // There is a prior for this node
                    (nodePrior.get(taxon)).logPdf(draw) == Double.NEGATIVE_INFINITY); // And draw is invalid under prior
            // TODO Currently only works for flat/truncated priors, make work for MVN

            treeModel.setMultivariateTrait(node, traitName, draw);

        } while (parameterPriorExists &&
                (parameterPrior.getLogLikelihood() == Double.NEGATIVE_INFINITY));

        return 0;  //To change body of implemented methods use File | Settings | File Templates.

    }

    private MeanPrecision operateNotRoot(NodeRef node) {

        double[][] precision = precisionMatrixParameter.getParameterAsMatrix();

        NodeRef parent = treeModel.getParent(node);

        double[] mean = new double[dim];

        double weight = 1.0 / traitModel.getRescaledBranchLength(node);

        double[] trait = treeModel.getMultivariateNodeTrait(parent, traitName);

        for (int i = 0; i < dim; i++)
            mean[i] = trait[i] * weight;

        double weightTotal = weight;
        for (int j = 0; j < treeModel.getChildCount(node); j++) {
            NodeRef child = treeModel.getChild(node, j);
            trait = treeModel.getMultivariateNodeTrait(child, traitName);
            weight = 1.0 / traitModel.getRescaledBranchLength(child);

            for (int i = 0; i < dim; i++)
                mean[i] += trait[i] * weight;

            weightTotal += weight;
        }

        for (int i = 0; i < dim; i++) {
            mean[i] /= weightTotal;
            for (int j = i; j < dim; j++)
                precision[j][i] = precision[i][j] *= weightTotal;
        }
        return new MeanPrecision(mean,precision);
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

        double[] trait = null;
        double weightTotal = 0.0;

        double[] weightedAverage = new double[dim];

        double[][] precision = precisionMatrixParameter.getParameterAsMatrix();
        
        for (int k = 0; k < treeModel.getChildCount(node); k++) {
            NodeRef child = treeModel.getChild(node, k);
            trait = treeModel.getMultivariateNodeTrait(child, traitName);
            final double weight = 1.0 / traitModel.getRescaledBranchLength(child);

            for (int i = 0; i < dim; i++) {
                for (int j=0; j<dim; j++)
                    weightedAverage[i] += precision[i][j] * weight * trait[j];
            }

            weightTotal += weight;
        }

        for (int i=0; i<dim; i++) {
            for (int j=0; j<dim; j++) {
                weightedAverage[i] += rootPriorPrecision[i][j] * rootPriorMean[j];
                precision[i][j]  = precision[i][j] * weightTotal + rootPriorPrecision[i][j];
            }
        }

        double[][] variance = new SymmetricMatrix(precision).inverse().toComponents();

        for (int i=0; i<dim; i++) {
            // todo: (FIXME) if code is correct in using the arbitrary last trait from 2 loops above that
            // todo: (FIXME) requires at least a comment.
            assert trait != null;
            
            trait[i] = 0;
            for (int j=0; j<dim; j++)
                trait[i] += variance[i][j] * weightedAverage[j];
        }

        return new MeanPrecision(trait,precision);
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

        private final String[] names = { GIBBS_OPERATOR, "internalTraitGibbsOperator" };

        public String[] getParserNames() { return names; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
    
            double weight = xo.getDoubleAttribute(WEIGHT);
            boolean onlyInternalNodes = xo.getAttribute(INTERNAL_ONLY, true);
            SampledMultivariateTraitLikelihood traitModel = (SampledMultivariateTraitLikelihood) xo.getChild(AbstractMultivariateTraitLikelihood.class);

            TraitGibbsOperator operator = new TraitGibbsOperator(traitModel, onlyInternalNodes);
            operator.setWeight(weight);

            // Get root prior

            XMLObject cxo = xo.getChild(ROOT_PRIOR);
            if (cxo != null) {

                MultivariateDistributionLikelihood rootPrior = (MultivariateDistributionLikelihood) cxo.getChild(MultivariateDistributionLikelihood.class);
                if( !(rootPrior.getDistribution() instanceof MultivariateDistribution))
                    throw new XMLParseException("Only multivariate normal priors allowed for Gibbs sampling the root trait");
               operator.setRootPrior((MultivariateNormalDistribution)rootPrior.getDistribution());
            }


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
                        operator.setTaxonPrior(treeModel.getTaxon(index),prior);
                        System.err.println("Adding truncated prior for taxon '"+treeModel.getTaxon(index)+"'");
                    }
                }
            }

            GeoSpatialCollectionModel collectionModel = (GeoSpatialCollectionModel) xo.getChild(GeoSpatialCollectionModel.class);
            if (collectionModel != null) {
                operator.setParameterPrior(collectionModel);
                System.err.println("Adding truncated prior '"+collectionModel.getId()+
                        "' for parameter '"+collectionModel.getParameter().getId()+"'");
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
                new ElementRule(GeoSpatialCollectionModel.class,true),
        };

    };

}