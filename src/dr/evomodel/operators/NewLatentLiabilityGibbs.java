


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
import dr.evomodel.treedatalikelihood.preorder.ConditionalVarianceAndTransform;
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

import java.util.List;

public class NewLatentLiabilityGibbs extends SimpleMCMCOperator {

    public static final String NEW_LATENT_LIABILITY_GIBBS_OPERATOR = "newlatentLiabilityGibbsOperator";
    public static final String TREE_MODEL = "treeModel";


    private final LatentTruncation latentLiability;

    private final CompoundParameter tipTraitParameter;

    private final TreeTrait<List<WrappedMeanPrecision>> fullConditionalDensity;


    private final Tree treeModel;
    private final int dim;

    public double[][] postMeans;
    public double[][] preMeans;
    public double[] preP;
    public double[] postP;
    private Parameter mask;
    private boolean hasMask = false;
    private int numFixed = 0;
    private int numUpdate = 0;
    private int[] doUpdate;
    private int[] dontUpdate;


    public NewLatentLiabilityGibbs(
            TreeDataLikelihood traitModel,
            LatentTruncation LatentLiability, CompoundParameter tipTraitParameter, Parameter mask,
            double weight, String traitName) {
        super();

        this.latentLiability = LatentLiability;
        this.tipTraitParameter = tipTraitParameter;
        this.treeModel = traitModel.getTree();
        ContinuousDataLikelihoodDelegate likelihoodDelegate = (ContinuousDataLikelihoodDelegate) traitModel.getDataLikelihoodDelegate();
        this.dim = likelihoodDelegate.getTraitDim();
        String fcdName = WrappedTipFullConditionalDistributionDelegate.getName(traitName);
        if (traitModel.getTreeTrait(fcdName) == null) {
            likelihoodDelegate.addWrappedFullConditionalDensityTrait(traitName);
        }
        this.fullConditionalDensity = castTreeTrait(traitModel.getTreeTrait(fcdName));
        this.mask = mask;
        if (mask != null) {
            hasMask = true;
        }

        postMeans = new double[treeModel.getNodeCount()][dim];
        preMeans = new double[treeModel.getNodeCount()][dim];
        preP = new double[treeModel.getNodeCount()];
        postP = new double[treeModel.getNodeCount()];

        dontUpdate = new int[dim];
        doUpdate = new int[dim];

        if (hasMask) {
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

        final List<WrappedMeanPrecision> allStatistics = fullConditionalDensity.getTrait(treeModel, treeModel.getExternalNode(pos));

        final WrappedMeanPrecision statistic;

        statistic = allStatistics.get(0);

        NodeRef node = treeModel.getExternalNode(pos);

        double logq = sampleNode2(node, statistic);

        tipTraitParameter.fireParameterChangedEvent();

        return logq;
    }

    public double[] getNodeTrait(NodeRef node) {
        int index = node.getNumber();
        double[] traitValue = tipTraitParameter.getParameter(index).getParameterValues();
        return traitValue;
    }

    public double getNodeTrait(NodeRef node, int entry) {
        int index = node.getNumber();
        double traitValue = tipTraitParameter.getParameter(index).getParameterValue(entry);
        return traitValue;
    }

    public void setNodeTrait(NodeRef node, double[] traitValue) {
        int index = node.getNumber();
        for (int i = 0; i < dim; i++) {

            tipTraitParameter.getParameter(index).setParameterValue(i, traitValue[i]);
        }
    }

    public void setNodeTrait(NodeRef node, int entry, double traitValue) {
        int index = node.getNumber();

        tipTraitParameter.getParameter(index).setParameterValue(entry, traitValue);

    }

    public double sampleNode2(NodeRef node, WrappedMeanPrecision statistics) {

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
            // TODO Automatically reject?
            return Double.NEGATIVE_INFINITY;
        }

        double pOld = distribution.logPdf(oldValue);

        double pNew = distribution.logPdf(value);

        double logq = pOld - pNew;

        return logq;
    }

    private void addMaskIfNeeded(ReadableVector jointmean, ReadableMatrix jointPrecision) {

        if(mask == null){
            return;
        } else {
            //todo
//            ConditionalVarianceAndTransform cVarianceJoint = new ConditionalVarianceAndTransform(
//                    new Matrix(jointGraphVariance), cMissingJoint, cNotMissingJoint);

        }
    }

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

        public final static String MASK = "mask";


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




