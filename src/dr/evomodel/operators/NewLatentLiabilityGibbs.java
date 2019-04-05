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
import dr.evomodel.treedatalikelihood.continuous.RepeatedMeasuresTraitDataModel;
import dr.evomodel.treedatalikelihood.preorder.ConditionalPrecisionAndTransform;
import dr.evomodel.treedatalikelihood.preorder.WrappedNormalSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.WrappedTipFullConditionalDistributionDelegate;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.*;
import dr.math.matrixAlgebra.missingData.MissingOps;
import dr.xml.*;
import org.ejml.data.D1Matrix64F;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.ArrayList;
import java.util.List;

import static org.ejml.alg.dense.mult.MatrixVectorMult.mult;
import static org.ejml.alg.dense.mult.MatrixVectorMult.multAdd;

public class NewLatentLiabilityGibbs extends SimpleMCMCOperator {

    private static final String NEW_LATENT_LIABILITY_GIBBS_OPERATOR = "newlatentLiabilityGibbsOperator";
    private static final String MAX_ATTEMPTS = "numAttempts";
    private static final String MISSING_BY_COLUMN = "missingByColumn";
    private static final String FORCE_ALL_MISSING = "forceAllMissing";

    private final LatentTruncation latentLiability;
    private final CompoundParameter tipTraitParameter;
    private final TreeTrait<List<WrappedNormalSufficientStatistics>> fullConditionalDensity;
    private final RepeatedMeasuresTraitDataModel repeatedMeasuresModel;

    private int maxAttempts;

    private final Tree treeModel;
    private final int dim;

    private Parameter mask;
    //    private MaskIndices maskIndices;
    private final MaskIndicesDelegate maskDelegate;
    private final Boolean missingByColumn;
    private final int[] needSampling;

    private double[] fcdMean;
    private double[][] fcdPrecision;
    private double[][] fcdVaraince;
    private double[] maskedMean;
    private double[][] maskedPrecision;


    public NewLatentLiabilityGibbs(
            TreeDataLikelihood treeDataLikelihood,
            LatentTruncation LatentLiability, CompoundParameter tipTraitParameter,
            RepeatedMeasuresTraitDataModel repeatedMeasuresModel, Parameter mask,
            double weight, String traitName, int maxAttempts, boolean missingByColumn) {
        super();

        this.latentLiability = LatentLiability;
        this.tipTraitParameter = tipTraitParameter;
        this.treeModel = treeDataLikelihood.getTree();
        this.repeatedMeasuresModel = repeatedMeasuresModel;
        ContinuousDataLikelihoodDelegate likelihoodDelegate = (ContinuousDataLikelihoodDelegate) treeDataLikelihood
                .getDataLikelihoodDelegate();
        this.dim = likelihoodDelegate.getTraitDim();
        String fcdName = WrappedTipFullConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(fcdName) == null) {
            likelihoodDelegate.addWrappedFullConditionalDensityTrait(traitName);
        }
        this.fullConditionalDensity = castTreeTrait(treeDataLikelihood.getTreeTrait(fcdName));

        this.missingByColumn = missingByColumn;
        this.mask = mask;
        this.maskDelegate = new MaskIndicesDelegate();
        this.needSampling = setupNeedSampling();
//        setupMaskDelegate();
//        setupmask();
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

        final int m = needSampling.length;
        final int pos = MathUtils.nextInt(m);


        NodeRef node = treeModel.getExternalNode(needSampling[pos]);

        final List<WrappedNormalSufficientStatistics> allStatistics = fullConditionalDensity.getTrait(treeModel, node);
        final WrappedNormalSufficientStatistics statistic = allStatistics.get(0);

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
        if (mask == null) {
            Parameter tip = tipTraitParameter.getParameter(index);
            for (int i = 0; i < dim; i++) {
                tip.setParameterValueQuietly(i, traitValue[i]);
            }
            tip.fireParameterChangedEvent(-1, Parameter.ChangeType.ALL_VALUES_CHANGED);
        } else {
            int j = 0;
            Parameter tip = tipTraitParameter.getParameter(index);

            for (int i : maskDelegate.getLatentIndices(node)) {
                tip.setParameterValueQuietly(i, traitValue[j]);
                j++;
            }
            tip.fireParameterChangedEvent(-1, Parameter.ChangeType.ALL_VALUES_CHANGED);
        }
    }

    private double sampleNode(NodeRef node, WrappedNormalSufficientStatistics statistics) {

        final int thisNumber = node.getNumber();
        final int[] obsInds = maskDelegate.getObservedIndices(node);
        final int obsDim = obsInds.length;

        if (obsDim == dim) {
            return 0;
        }

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

        if (repeatedMeasuresModel != null) {
            //TODO: preallocate memory for all these matrices/vectors
            DenseMatrix64F Q = new DenseMatrix64F(fcdPrecision); //storing original fcd precision
            double[] tipPartial = repeatedMeasuresModel.getTipPartial(thisNumber, false);
            int offset = dim;
            DenseMatrix64F P = MissingOps.wrap(tipPartial, offset, dim, dim);

            DenseMatrix64F preOrderMean = new DenseMatrix64F(dim, 1);
            for (int i = 0; i < dim; i++) {
                preOrderMean.set(i, 0, fcdMean[i]);
            }
            DenseMatrix64F dataMean = new DenseMatrix64F(dim, 1);
            for (int i = 0; i < dim; i++) {
                dataMean.set(i, 0, tipPartial[i]);
            }

            D1Matrix64F bufferMean = new DenseMatrix64F(dim, 1);

            mult(Q, preOrderMean, bufferMean); //bufferMean = Q * preOrderMean
            multAdd(P, dataMean, bufferMean); //bufferMean = Q * preOderMean + P * dataMean


            CommonOps.addEquals(P, Q); //P = P + Q
            DenseMatrix64F V = new DenseMatrix64F(dim, dim);
            CommonOps.invert(P, V); //V = inv(P + Q)
            mult(V, bufferMean, dataMean); // dataMean = inv(P + Q) * (Q * preOderMean + P * dataMean)

            for (int i = 0; i < dim; i++) {
                fcdMean[i] = dataMean.get(i);
                for (int j = 0; j < dim; j++) {
                    fcdPrecision[i][j] = P.get(i, j);
                }

            }
        }

        MultivariateNormalDistribution fullDistribution = new MultivariateNormalDistribution(fcdMean, fcdPrecision); //TODO: should this not be declared until 'else' statement?
        MultivariateNormalDistribution drawDistribution;

        if (mask != null && obsDim > 0) {
            addMaskOnContiuousTraitsPrecisionSpace(thisNumber);
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

//    private void addMaskOnContiuousTraits(int nodeNumber) {
//
//        double[] currentValues = new double[dim];
//
//        for (int i = 0; i < currentValues.length; ++i) {
//            currentValues[i] = tipTraitParameter.getParameterValues()[nodeNumber * dim + i];
//        }
//
//        ConditionalVarianceAndTransform cVarianceJoint = new ConditionalVarianceAndTransform(
//                new Matrix(fcdVaraince), maskIndices.discreteIndices, maskIndices.continuousIndex);
//
//        maskedPrecision = cVarianceJoint.getConditionalVariance().inverse().toComponents();
//        maskedMean = cVarianceJoint.getConditionalMean(currentValues, 0, fcdMean, 0);
//    }

    private void addMaskOnContiuousTraitsPrecisionSpace(int nodeNumber) {

        double[] currentValues = new double[dim];

        for (int i = 0; i < currentValues.length; ++i) {
            currentValues[i] = tipTraitParameter.getParameterValues()[nodeNumber * dim + i];
        }


        ConditionalPrecisionAndTransform cVarianceJoint = new ConditionalPrecisionAndTransform(
                new Matrix(fcdPrecision), maskDelegate.getLatentIndices(nodeNumber), maskDelegate.getObservedIndices(nodeNumber));

        maskedPrecision = cVarianceJoint.getConditionalPrecision().toComponents();
        maskedMean = cVarianceJoint.getConditionalMean(currentValues, 0, fcdMean, 0);
    }

    private int[] convertListToArray(List<Integer> listResult) { //todo this shouldn't be here...
        int[] result = new int[listResult.size()];
        int i = 0;
        for (int num : listResult) {
            result[i++] = num;
        }
        return result;
    }

//    private void setupmask() {
//
//        if (mask != null) {
//
//            List<Integer> missingIndex = new ArrayList<Integer>();
//            List<Integer> notmissingIndex = new ArrayList<Integer>();
//
//            for (int i = 0; i < dim; ++i) {
//
//                if (mask.getParameterValue(i) == 1.0) {
//                    missingIndex.add(i);
//                } else {
//                    notmissingIndex.add(i);
//                }
//            }
//
//            int[] cMissingJoint = convertListToArray(missingIndex);
//            int[] cNotMissingJoint = convertListToArray(notmissingIndex);
//
//            maskIndices = new MaskIndices(cMissingJoint, cNotMissingJoint);
//        }
//    }

    protected class MaskIndices {

        final int[] discreteIndices;
        final int[] continuousIndex;

        private MaskIndices(int[] latentindex, int[] contindex) {
            this.discreteIndices = latentindex;
            this.continuousIndex = contindex;
        }
    }


    private class MaskIndicesDelegate {

        int[] latentColumns = null;
        int[] observedColumns = null;

        private MaskIndicesDelegate() {
            if (mask != null) {
                if (missingByColumn) {

                    List<Integer> missingIndex = new ArrayList<Integer>();
                    List<Integer> notmissingIndex = new ArrayList<Integer>();

                    for (int i = 0; i < dim; ++i) {

                        if (mask.getParameterValue(i) == 1.0) {
                            missingIndex.add(i);
                        } else {
                            notmissingIndex.add(i);
                        }
                    }

                    this.latentColumns = convertListToArray(missingIndex);
                    this.observedColumns = convertListToArray(notmissingIndex);
                }
            }


        }

        private int[] getLatentIndices(NodeRef nodeRef) {
            return getLatentIndices(nodeRef.getNumber());
        }

        private int[] getLatentIndices(int nodeNumber) {
            if (missingByColumn) {
                return this.latentColumns;
            } else {

                int offset = dim * nodeNumber;
                List<Integer> latentArray = new ArrayList<Integer>();

                for (int i = offset; i < offset + dim; i++) {
                    if (mask.getParameterValue(i) == 1.0) {
                        latentArray.add(i - offset);
                    }
                }

                return convertListToArray(latentArray);

            }
        }

        private int[] getObservedIndices(int nodeNumber) {
            if (missingByColumn) {
                return this.observedColumns;
            } else {

                int offset = dim * nodeNumber;
                List<Integer> obsArray = new ArrayList<Integer>();

                for (int i = offset; i < offset + dim; i++) {
                    if (mask.getParameterValue(i) == 0.0) {
                        obsArray.add(i - offset);
                    }
                }

                return convertListToArray(obsArray);
            }
        }

        private int[] getObservedIndices(NodeRef nodeRef) {
            return getObservedIndices(nodeRef.getNumber());
        }


    }

    private int[] setupNeedSampling() {
        int n = treeModel.getExternalNodeCount();
        List<Integer> sampleList = new ArrayList<Integer>();

        for (int i = 0; i < n; i++) {

            int obsDim = maskDelegate.getObservedIndices(i).length;

            if (obsDim != dim) {
                sampleList.add(i);
            }
        }

        return convertListToArray(sampleList);
    }

    @SuppressWarnings("unchecked")
    private TreeTrait<List<WrappedNormalSufficientStatistics>> castTreeTrait(TreeTrait trait) {
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
        private final static String PARTIALS_PROVIDER = "partialsProvider";


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

            boolean missingByColumn = xo.getAttribute(MISSING_BY_COLUMN, true);


            Parameter mask = null;

            if (xo.hasChildNamed(MASK)) {
                mask = (Parameter) xo.getElementFirstChild(MASK);
            }

            RepeatedMeasuresTraitDataModel repeatedMeasuresModel = null;
            if (xo.hasChildNamed(PARTIALS_PROVIDER)) {
                repeatedMeasuresModel = (RepeatedMeasuresTraitDataModel) xo.getElementFirstChild(PARTIALS_PROVIDER);
            }

            if (xo.getAttribute(FORCE_ALL_MISSING, false)) {
                int dim = traitModel.getDataLikelihoodDelegate().getTraitDim();
                mask = new Parameter.Default(dim);
                for (int i = 0; i < dim; i++) {
                    mask.setParameterValue(i, 1.0);
                }
                missingByColumn = true;

            }

            return new NewLatentLiabilityGibbs(traitModel, LLModel, tipTraitParameter, repeatedMeasuresModel, mask, weight, "latent",
                    numAttempts, missingByColumn);
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
                AttributeRule.newBooleanRule(MISSING_BY_COLUMN, true),
//                AttributeRule.newBooleanRule(FORCE_ALL_MISSING, false),
                new ElementRule(TreeDataLikelihood.class, "The model for the latent random variables"),
                new ElementRule(LatentTruncation.class, "The model that links latent and observed variables"),
                new ElementRule(MASK, dr.inference.model.Parameter.class, "Mask: 1 for latent variables that should " +
                        "be sampled", true),
                new ElementRule(CompoundParameter.class, "The parameter of tip locations from the tree"),
                new ElementRule(PARTIALS_PROVIDER, RepeatedMeasuresTraitDataModel.class,
                        "Provides information about model extensions", true)

        };
    };
}




