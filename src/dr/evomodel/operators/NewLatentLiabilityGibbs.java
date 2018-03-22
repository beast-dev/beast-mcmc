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

import java.util.ArrayList;
import java.util.List;

public class NewLatentLiabilityGibbs extends SimpleMCMCOperator {

    private static final String NEW_LATENT_LIABILITY_GIBBS_OPERATOR = "newlatentLiabilityGibbsOperator";
    private static final String MAX_ATTEMPTS = "numAttempts";

    private final LatentTruncation latentLiability;
    private final CompoundParameter tipTraitParameter;
    private final TreeTrait<List<WrappedMeanPrecision>> fullConditionalDensity;

    private int maxAttempts;

    private final Tree treeModel;
    private final int dim;

    private Parameter mask;
    private MaskIndices maskIndices;

    private double[] fcdMean;
    private double[][] fcdPrecision;
    private double[][] fcdVaraince;
    private double[] maskedMean;
    private double[][] maskedPrecision;


    public NewLatentLiabilityGibbs(
            TreeDataLikelihood treeDataLikelihood,
            LatentTruncation LatentLiability, CompoundParameter tipTraitParameter, Parameter mask,
            double weight, String traitName, int maxAttempts) {
        super();

        this.latentLiability = LatentLiability;
        this.tipTraitParameter = tipTraitParameter;
        this.treeModel = treeDataLikelihood.getTree();
        ContinuousDataLikelihoodDelegate likelihoodDelegate = (ContinuousDataLikelihoodDelegate) treeDataLikelihood
                .getDataLikelihoodDelegate();
        this.dim = likelihoodDelegate.getTraitDim();
        String fcdName = WrappedTipFullConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(fcdName) == null) {
            likelihoodDelegate.addWrappedFullConditionalDensityTrait(traitName);
        }
        this.fullConditionalDensity = castTreeTrait(treeDataLikelihood.getTreeTrait(fcdName));
        this.mask = mask;
        setupmask();
        this.fcdMean = new double[dim];
        this.fcdVaraince = new double[dim][dim];
        this.fcdPrecision = new double[dim][dim];

        this.maxAttempts = maxAttempts;

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

        double logq = sampleNode(node, statistic);
        tipTraitParameter.fireParameterChangedEvent();
        return logq;
    }

    private double[] getNodeTrait(NodeRef node) {
        int index = node.getNumber();
        return tipTraitParameter.getParameter(index).getParameterValues();
    }

    private void setNodeTrait(NodeRef node, double[] traitValue) {

        int index = node.getNumber();
        if (maskIndices == null) {
            for (int i = 0; i < dim; i++) {
                tipTraitParameter.getParameter(index).setParameterValue(i, traitValue[i]);
            }
        } else {
            int j = 0;
            for (int i : maskIndices.missingIndices) {
                tipTraitParameter.getParameter(index).setParameterValue(i, traitValue[j]);
                j++;
            }
        }
    }

    private double sampleNode(NodeRef node, WrappedMeanPrecision statistics) {

        final int thisNumber = node.getNumber();

        ReadableVector mean = statistics.getMean();
        ReadableMatrix thisP = statistics.getPrecision();
        double precisionScalar = statistics.getPrecisionScalar();

        for (int i = 0; i < mean.getDim(); ++i) {
            fcdMean[i] = mean.get(i);
        }

        for (int i = 0; i < mean.getDim(); ++i) {
            for (int j = 0; j < mean.getDim(); ++j) {
                fcdPrecision[i][j] = thisP.get(i, j) * precisionScalar;
            }
        }

        fcdVaraince = (new Matrix(fcdPrecision)).inverse().toComponents();

        MultivariateNormalDistribution fullDistribution = new MultivariateNormalDistribution(fcdMean, fcdPrecision);
        MultivariateNormalDistribution drawDistribution;

        if (maskIndices != null) {
            addMaskIfNeeded();
            drawDistribution = new MultivariateNormalDistribution(maskedMean, maskedPrecision);
        } else {
            drawDistribution = fullDistribution;
        }

        double[] oldValue = getNodeTrait(node);

        int attempt = 0;
        boolean validTip = false;

        while (!validTip & attempt < maxAttempts) {

            setNodeTrait(node, drawDistribution.nextMultivariateNormal());

            if (latentLiability.validTraitForTip(thisNumber)) {
                validTip = true;
            }
            attempt++;
        }

        if (attempt == maxAttempts) {
            return Double.NEGATIVE_INFINITY;
        }

        double[] newValue = getNodeTrait(node);

        return fullDistribution.logPdf(oldValue) - fullDistribution.logPdf(newValue);
    }

    private void addMaskIfNeeded() {

        ConditionalVarianceAndTransform cVarianceJoint = new ConditionalVarianceAndTransform(
                new Matrix(fcdVaraince), maskIndices.missingIndices, maskIndices.notMissingIndex);

        maskedPrecision = cVarianceJoint.getConditionalVariance().inverse().toComponents();
        maskedMean = cVarianceJoint.getConditionalMean(tipTraitParameter.getParameterValues(), 0, fcdMean, 0);
    }

    private int[] convertListToArray(List<Integer> listResult) { //todo this shouldn't be here...
        int[] result = new int[listResult.size()];
        int i = 0;
        for (int num : listResult) {
            result[i++] = num;
        }
        return result;
    }

    private void setupmask() {

        if (mask != null) {

            List<Integer> missingIndex = new ArrayList<Integer>();
            List<Integer> notmissingIndex = new ArrayList<Integer>();

            for (int i = 0; i < dim; ++i) {

                if (mask.getParameterValue(i) == 1.0) {
                    missingIndex.add(i);
                } else {
                    notmissingIndex.add(i);
                }
            }

            int[] cMissingJoint = convertListToArray(missingIndex);
            int[] cNotMissingJoint = convertListToArray(notmissingIndex);

            maskIndices = new MaskIndices(cMissingJoint, cNotMissingJoint);
        }
    }

    protected class MaskIndices {

        final int[] missingIndices;
        final int[] notMissingIndex;

        private MaskIndices(int[] latentindex, int[] contindex) {
            this.missingIndices = latentindex;
            this.notMissingIndex = contindex;
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
                                "\t 1 conjugate multivariateTraitLikelihood, 1 latentLiabilityLikelihood and one " +
                                "parameter \n"
                );
            }

            double weight = xo.getDoubleAttribute(WEIGHT);

            TreeDataLikelihood traitModel = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
            LatentTruncation LLModel = (LatentTruncation) xo.getChild(LatentTruncation.class);
            CompoundParameter tipTraitParameter = (CompoundParameter) xo.getChild(CompoundParameter.class);
            int numAttempts = xo.getAttribute(MAX_ATTEMPTS, 100000);

            Parameter mask = null;

            if (xo.hasChildNamed(MASK)) {
                mask = (Parameter) xo.getElementFirstChild(MASK);
            }

            return new NewLatentLiabilityGibbs(traitModel, LLModel, tipTraitParameter, mask, weight, "latent",
                    numAttempts);
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
                new ElementRule(MASK, dr.inference.model.Parameter.class, "Mask: 1 for latent variables that should " +
                        "be sampled", true),
                new ElementRule(CompoundParameter.class, "The parameter of tip locations from the tree")

        };
    };
}




