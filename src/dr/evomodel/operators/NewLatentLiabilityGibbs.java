


/*
 * NewLatentLiabilityGibbs.java
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


/*
 * Gabriela Cybis
 *
 *   Gibbs operator for latent variable in latent liability model
 */

package dr.evomodel.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.continuous.LatentTruncation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.preorder.WrappedMeanPrecision;
import dr.evomodel.treedatalikelihood.preorder.WrappedTipFullConditionalDistributionDelegate;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.*;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class NewLatentLiabilityGibbs extends SimpleMCMCOperator {

    private static final String NEW_LATENT_LIABILITY_GIBBS_OPERATOR = "newlatentLiabilityGibbsOperator";

    private final LatentTruncation latentLiability;
    private final CompoundParameter tipTraitParameter;
    private final TreeTrait<List<WrappedMeanPrecision>> fullConditionalDensity;

    private final Tree treeModel;
    private final int dim;

    private Parameter mask;
    private int numFixed = 0;
    private int numUpdate = 0;
    private int[] doUpdate;
    private int[] dontUpdate;


    public NewLatentLiabilityGibbs(
            TreeDataLikelihood treeDataLikelihood,
            LatentTruncation LatentLiability, CompoundParameter tipTraitParameter, Parameter mask,
            double weight, String traitName) {
        super();

        this.latentLiability = LatentLiability;
        this.tipTraitParameter = tipTraitParameter;
        this.treeModel = treeDataLikelihood.getTree();
        ContinuousDataLikelihoodDelegate likelihoodDelegate = (ContinuousDataLikelihoodDelegate) treeDataLikelihood.getDataLikelihoodDelegate();
        this.dim = likelihoodDelegate.getTraitDim();
        String fcdName = WrappedTipFullConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(fcdName) == null) {
            likelihoodDelegate.addWrappedFullConditionalDensityTrait(traitName);
        }
        this.fullConditionalDensity = castTreeTrait(treeDataLikelihood.getTreeTrait(fcdName));
        this.mask = mask;

        dontUpdate = new int[dim];
        doUpdate = new int[dim];

        if (mask != null) {
            for (int i = 0; i < dim; i++) {
                if (mask.getParameterValue(i) == 0.0) {
                    dontUpdate[numFixed] = i;
                    numFixed++;
                } else {
                    doUpdate[numUpdate] = i;
                    numUpdate++;
                }

            }
        }
        setWeight(weight);
    }


    public int getStepCount() {
        return 1;
    }

    public double doOperation() {

        final int pos = MathUtils.nextInt(treeModel.getExternalNodeCount());
        NodeRef node = treeModel.getExternalNode(pos);

        final List<WrappedMeanPrecision> allStatistics = fullConditionalDensity.getTrait(treeModel, node);
        final WrappedMeanPrecision statistic = allStatistics.get(0);

        double logq = sampleNode2(node, statistic);
        tipTraitParameter.fireParameterChangedEvent();
        return logq;
    }

    private double[] getNodeTrait(NodeRef node) {
        int index = node.getNumber();
        return tipTraitParameter.getParameter(index).getParameterValues();
    }

//    public double getNodeTrait(NodeRef node, int entry) {
//        int index = node.getNumber();
//        double traitValue = tipTraitParameter.getParameter(index).getParameterValue(entry);
//        return traitValue;
//    }

    private void setNodeTrait(NodeRef node, double[] traitValue) {
        int index = node.getNumber();
        for (int i = 0; i < dim; i++) {
            tipTraitParameter.getParameter(index).setParameterValue(i, traitValue[i]);
        }
    }

//    private void setNodeTrait(NodeRef node, int entry, double traitValue) {
//        int index = node.getNumber();
//
//        tipTraitParameter.getParameter(index).setParameterValue(entry, traitValue);
//
//    }

    private double sampleNode2(NodeRef node, WrappedMeanPrecision statistics) {

        final int thisNumber = node.getNumber();

        ReadableVector mean = statistics.getMean();
        ReadableMatrix thisP = statistics.getPrecision();

        //todo: add mask; pass the joint conditional mean and precision for both discrete and continuous traits.
        addMaskIfNeeded(mean, thisP);

        double[] meanVector = new double[mean.getDim()];
        double[][] precisionMatrix = new double[mean.getDim()][mean.getDim()];

        for (int i = 0; i < mean.getDim(); ++i) {
            meanVector[i] = mean.get(i);
        }

        for (int i = 0; i < mean.getDim(); ++i) {
            for (int j = 0; j < mean.getDim(); ++j) {
                precisionMatrix[i][j] = thisP.get(i, j);
            }
        }

        double[] oldValue = getNodeTrait(node);
        double[] value = new double[oldValue.length];

        int attempt = 0;
        boolean validTip = false;

        //	  Sample it all traits together as one multivariate normal

        MultivariateNormalDistribution distribution = new MultivariateNormalDistribution(meanVector, precisionMatrix);

        final int max = 10000;

        while (!validTip & attempt < max) {

            value = distribution.nextMultivariateNormal();

            setNodeTrait(node, value);

            if (latentLiability.validTraitForTip(thisNumber)) {
                validTip = true;
            }
            attempt++;
        }

        if (attempt == max) {
            return Double.NEGATIVE_INFINITY;
        }

        double pOld = distribution.logPdf(oldValue);

        double pNew = distribution.logPdf(value);

        return pOld - pNew;
    }

    private void addMaskIfNeeded(ReadableVector jointMean, ReadableMatrix jointPrecision) {

        if (mask != null) {

            List<Integer> missingIndex = new ArrayList<Integer>();
            List<Integer> notmissingIndex = new ArrayList<Integer>();

            for(int i = 0; i < dim; ++i){
                if (mask.getParameterValue(i) == 1.0) {
                    missingIndex.add(i);
                } else {
                    notmissingIndex.add(i);
                }
            }
        }
    }

    @SuppressWarnings("unchecked")
    private TreeTrait<List<WrappedMeanPrecision>> castTreeTrait(TreeTrait trait) {
        return trait;
    }

    public String getPerformanceSuggestion() {
        return null;
    }


    public String getOperatorName() {
        return NEW_LATENT_LIABILITY_GIBBS_OPERATOR;
    }


    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        private final static String MASK = "mask";


        public String getParserName() {
            return NEW_LATENT_LIABILITY_GIBBS_OPERATOR;
        }


        public Object parseXMLObject(XMLObject xo) throws XMLParseException {


            if (xo.getChildCount() < 3) {
                throw new XMLParseException(
                        "Element with id = '" + xo.getName() + "' should contain:\n" +
                                "\t 1 conjugate multivariateTraitLikelihood, 1 latentLiabilityLikelihood and one parameter \n"
                );
            }

            double weight = xo.getDoubleAttribute(WEIGHT);

            TreeDataLikelihood traitModel = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
            LatentTruncation LLModel = (LatentTruncation) xo.getChild(LatentTruncation.class);
            CompoundParameter tipTraitParameter = (CompoundParameter) xo.getChild(CompoundParameter.class);

            Parameter mask = null;

            if (xo.hasChildNamed(MASK)) {
                mask = (Parameter) xo.getElementFirstChild(MASK);
            }

            return new NewLatentLiabilityGibbs(traitModel, LLModel, tipTraitParameter, mask, weight, "latent");
        }

        public String getParserDescription() {
            return "This element returns a gibbs sampler on tip latent trais for latent liability model.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(WEIGHT),
                new ElementRule(TreeDataLikelihood.class, "The model for the latent random variables"),
                new ElementRule(LatentTruncation.class, "The model that links latent and observed variables"),
                new ElementRule(MASK, dr.inference.model.Parameter.class, "Mask: 1 for latent variables that should be sampled", true),
                new ElementRule(CompoundParameter.class, "The parameter of tip locations from the tree")

        };
    };
}




