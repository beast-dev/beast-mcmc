/*
 * SubstitutionModelCrossProductDelegate.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.treedatalikelihood.discrete;

import beagle.Beagle;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.preorder.AbstractBeagleGradientDelegate;
import dr.evomodel.treedatalikelihood.preorder.AdjointMethods;
import dr.evomodel.treedatalikelihood.preorder.DiscretePartialsType;

import dr.evomodel.treedatalikelihood.discrete.discretetreedataLikelihood.ComplexBlockKernelUtils;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.Arrays;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class SpectralBeagleCrossProductDelegate extends AbstractBeagleGradientDelegate {

    private final String name;
    private final Tree tree;
    private final BranchRateModel branchRateModel;
    private final BranchModel branchModel;
    private final int stateCount;
    private final int substitutionModelCount;

    private static final String GRADIENT_TRAIT_NAME = "substitutionModelCrossProductGradient";
    private static final String USE_COMPLEX_BLOCK_KERNEL_UTILS_PROPERTY =
            "dr.evomodel.treedatalikelihood.discrete.SpectralBeagleCrossProductDelegate.useComplexBlockKernelUtils";
    private static final boolean DEFAULT_USE_COMPLEX_BLOCK_KERNEL_UTILS =
            Boolean.parseBoolean(System.getProperty(USE_COMPLEX_BLOCK_KERNEL_UTILS_PROPERTY, "true"));

    private final double[] postOrderPartials;
    private final double[] preOrderPartials;
    private final double[] intermediate;
    private final double[] eigenBasisAccum;
    private final double[] midBuffer;
    private final ComplexBlockKernelUtils.ComplexKernelPlan kernelPlan;
    private final boolean useComplexBlockKernelUtils;
    private int[] postBufferIndices;
    private int[] preBufferIndices;
    private int[] matrixBufferIndices;
    private double[] branchLengths;

    private static final boolean BEAGLE_OVERRIDE = true;

    public SpectralBeagleCrossProductDelegate(String name,
                                              Tree tree,
                                              BeagleDataLikelihoodDelegate likelihoodDelegate,
                                              BranchRateModel branchRateModel,
                                              int stateCount) {

        this(name, tree, likelihoodDelegate, branchRateModel, stateCount, DEFAULT_USE_COMPLEX_BLOCK_KERNEL_UTILS);
    }

    public SpectralBeagleCrossProductDelegate(String name,
                                              Tree tree,
                                              BeagleDataLikelihoodDelegate likelihoodDelegate,
                                              BranchRateModel branchRateModel,
                                              int stateCount,
                                              boolean useComplexBlockKernelUtils) {

        super(name, tree, likelihoodDelegate);
        this.name = name;
        this.tree = tree;
        this.stateCount = stateCount;
        this.branchRateModel = branchRateModel;
        this.branchModel = likelihoodDelegate.getBranchModel();
        this.substitutionModelCount = branchModel.getSubstitutionModels().size();
        this.useComplexBlockKernelUtils = useComplexBlockKernelUtils;

        this.eigenBasisAccum = new double[stateCount * stateCount];
        this.midBuffer = new double[stateCount * stateCount];

        final int partialLength = categoryCount * patternCount * stateCount;
        this.postOrderPartials = new double[partialLength];
        this.preOrderPartials = new double[partialLength];
        this.intermediate = new double[partialLength];
        this.kernelPlan = useComplexBlockKernelUtils ?
                new ComplexBlockKernelUtils.ComplexKernelPlan(stateCount) : null;
    }

    @Override
    protected DiscretePartialsType getPreOrderType() {
        return DiscretePartialsType.TOP;
    }

    private double getBranchLength(NodeRef node) {

        final double branchRate;
        synchronized (branchRateModel) {
            branchRate = branchRateModel.getBranchRate(tree, node);
        }

        double parentHeight = tree.getNodeHeight(tree.getParent(node));
        double nodeHeight = tree.getNodeHeight(node);

        return branchRate * (parentHeight - nodeHeight);
    }

    @Override
    protected int getGradientLength() {
        return stateCount * stateCount * substitutionModelCount;
    }

    private int getRootPostOrderBuffer() {
        NodeRef root = tree.getRoot();
        return getPostOrderPartialIndex(root.getNumber());
    }

    private int coverWholeTree(int[] postBufferIndices,
                               int[] preBufferIndices,
                               int[] matrixBufferIndices,
                               double[] branchLengths) {
        int u = 0;
        for (int nodeNum = 0; nodeNum < tree.getNodeCount(); nodeNum++) {
            NodeRef node = tree.getNode(nodeNum);
            if (!tree.isRoot(tree.getNode(nodeNum))) {
                postBufferIndices[u] = getPostOrderPartialIndex(nodeNum);
                preBufferIndices[u]  = getPreOrderPartialIndex(nodeNum);
                matrixBufferIndices[u] = likelihoodDelegate.getEvolutionaryProcessDelegate().getMatrixIndex(nodeNum);
                branchLengths[u] = getBranchLength(node);
                u++;
            }
        }
        return tree.getNodeCount() - 1;
    }

    private int coverPartialTree(int modelNumber,
                                 int[] postBufferIndices,
                                 int[] preBufferIndices,
                                 double[] branchLengths) {
        int u = 0;
        for (int nodeNum = 0; nodeNum < tree.getNodeCount(); nodeNum++) {
            NodeRef node = tree.getNode(nodeNum);
            if (!tree.isRoot(tree.getNode(nodeNum))) {

                BranchModel.Mapping mapping = branchModel.getBranchModelMapping(node);
                int[] order = mapping.getOrder();
                double[] weights = mapping.getWeights();

                for (int k = 0; k < order.length; ++k) {
                    if (order[k] == modelNumber) {
                        postBufferIndices[u] = getPostOrderPartialIndex(nodeNum);
                        preBufferIndices[u]  = getPreOrderPartialIndex(nodeNum);
                        branchLengths[u] = getBranchLength(node) * relativeWeight(k, weights);
                        u++;
                    }
                }
            }
        }
        return u;
    }

    private double relativeWeight(int k, double[] weights) {
        double sum = 0.0;
        for (double w : weights) {
            sum += w;
        }
        return weights[k] / sum;
    }

    @Override
    protected void getNodeDerivatives(Tree tree, double[] first, double[] second) {

        final int length = stateCount * stateCount;

        assert (first.length >= length * substitutionModelCount);
        assert (second == null || second.length >= stateCount * stateCount * substitutionModelCount);

        if (second != null) {
            throw new RuntimeException("Not yet implemented");
        }

        ensureBranchBuffers(tree.getNodeCount() - 1);

        List<SubstitutionModel> substitutionModels = likelihoodDelegate.getBranchModel().getSubstitutionModels();

        if (substitutionModelCount == 1) {

            EigenDecomposition ed = substitutionModels.get(0).getEigenDecomposition();
            // ted.getInverseEigenVectors() is R^T for the original eigensystem;
            // this rotates branch-top pre-order messages into the kernel basis.
            EigenDecomposition ted = ed.transpose();

            if (useComplexBlockKernelUtils && !BEAGLE_OVERRIDE) {
                ComplexBlockKernelUtils.fillStructure(kernelPlan, ed);
            }
            Arrays.fill(eigenBasisAccum, 0, eigenBasisAccum.length, 0.0);

            final double[] patternWeights  = patternList.getPatternWeights();
            final double[] categoryWeights = siteRateModel.getCategoryProportions();
            final double[] categoryRates   = siteRateModel.getCategoryRates();

            double[] buffer = new double[first.length];

            int count = coverWholeTree(postBufferIndices, preBufferIndices, matrixBufferIndices, branchLengths);
            calculateCrossProductDifferentials(0, count,
                    postBufferIndices, preBufferIndices, matrixBufferIndices, ed, ted, branchLengths,
                    patternWeights, categoryWeights, categoryRates, buffer);

//            System.err.println("beagle-output : " + new WrappedVector.Raw(buffer));
//            System.err.println("java-output   : " + new WrappedVector.Raw(eigenBasisAccum));
//            System.exit(-1);

            if (BEAGLE_OVERRIDE) {
                rotateIntoOutput(buffer, ed.getEigenVectors(), ed.getInverseEigenVectors(), first, 0);
            } else {
                rotateIntoOutput(eigenBasisAccum, ed.getEigenVectors(), ed.getInverseEigenVectors(), first, 0);
            }
        } else {

            throw new RuntimeException("Not yet implemented");

//            double[] buffer = new double[length];
//
//            for (int i = 0; i < substitutionModelCount; ++i) {
//
//                Arrays.fill(buffer, 0, buffer.length, 0.0);
//                int count = coverPartialTree(i, postBufferIndices, preBufferIndices,
//                        branchLengths);
//                beagle.calculateCrossProductDifferentials(postBufferIndices, preBufferIndices,
//                        new int[] { 0 }, new int[] { 0 },
//                        branchLengths,
//                        count,
//                        buffer, null);
//
//                System.arraycopy(buffer, 0, first, i * length, length);
//            }
        }
        // TOOD handle `firstSquared` and `second`
    }

    private void ensureBranchBuffers(int branchCount) {
        if (postBufferIndices == null || postBufferIndices.length != branchCount) {
            postBufferIndices = new int[branchCount];
            preBufferIndices = new int[branchCount];
            matrixBufferIndices = new int[branchCount];
            branchLengths = new double[branchCount];
        }
    }

    private void calculateCrossProductDifferentials(int start,
                                                    int end,
                                                    int[] postBufferIndices,
                                                    int[] preBufferIndices,
                                                    int[] matrixBufferIndices,
                                                    EigenDecomposition ed,
                                                    EigenDecomposition ted,
                                                    double[] branchLengths,
                                                    double[] patternWeights,
                                                    double[] categoryWeights,
                                                    double[] categoryRates,
                                                    double[] first) {

        if (BEAGLE_OVERRIDE) {
            beagle.calculateAdjointCrossProductDifferentials(
                    postBufferIndices, preBufferIndices,
                    matrixBufferIndices,
                    new int[]{0}, new int[]{0}, getRootPostOrderBuffer(), 0,
                    branchLengths.length, first, null);

        } else {
            for (int i = start; i < end; ++i) {
                calculateCrossProductDifferentials(postBufferIndices[i], preBufferIndices[i],
                        ed, ted, branchLengths[i], patternWeights, categoryWeights, categoryRates);
            }
        }
    }

    private void calculateCrossProductDifferentials(int postBufferIndex,
                                                    int preBufferIndex,
                                                    EigenDecomposition ed,
                                                    EigenDecomposition ted,
                                                    double branchLength,
                                                    double[] patternWeights,
                                                    double[] categoryWeights,
                                                    double[] categoryRates) {

        getRotatedPartial(postBufferIndex,
                ed.getInverseEigenVectors(),
                postOrderPartials, intermediate);
        getRotatedPartial(preBufferIndex,
                ted.getInverseEigenVectors(),
                preOrderPartials, intermediate);

        // TODO off load this to BEAGLE (after profiling)

        final double[] eigenValues = ed.getEigenValues();
        int offset = 0;
        for (int i = 0; i < categoryCount; ++i) {
            final double wc = categoryWeights[i];
            final double tc = branchLength * (categoryRates == null ? 1.0 : categoryRates[i]);
            if (useComplexBlockKernelUtils) {
                ComplexBlockKernelUtils.fillTransitionCoefficients(kernelPlan, ed, tc);
            }

            for (int j = 0; j < patternCount; ++j) {
                final double wp = patternWeights[j];
                if (wp == 0.0) {
                    offset += stateCount;
                    continue;
                }

                final double denom = branchLikelihoodInRotatedBasis(eigenValues, tc, offset);
                if (denom <= 0.0 || Double.isNaN(denom) || Double.isInfinite(denom)) {
                    offset += stateCount;
                    continue;
                }

                final double scale = (wp * wc) / denom;

//                System.err.println("java likelihood - " + denom + " with scale " + scale);
//                System.err.println("branchLength: " + branchLength);

                accumulateEigenBasisGradient(eigenValues, tc, offset, scale);
                offset += stateCount;
            }
        }

//        System.err.println("post-order-rotated: " + new WrappedVector.Raw(postOrderPartials));
//        System.err.println("pre-order-rotated : " + new WrappedVector.Raw(preOrderPartials));
//        System.err.println("eigen-basis-accum : " + new WrappedVector.Raw(eigenBasisAccum));
    }

    private double branchLikelihoodInRotatedBasis(double[] eigenValues,
                                                  double time,
                                                  int offset) {
        if (useComplexBlockKernelUtils) {
            return ComplexBlockKernelUtils.blockDiagonalTransitionInnerProduct(
                    kernelPlan, preOrderPartials, offset, postOrderPartials, offset);
        }

        return AdjointMethods.branchLikelihoodInRotatedBasis(
                eigenValues, time, preOrderPartials, offset, postOrderPartials, offset, stateCount);
    }

    private void accumulateEigenBasisGradient(double[] eigenValues,
                                              double time,
                                              int offset,
                                              double scale) {
        if (useComplexBlockKernelUtils) {
            ComplexBlockKernelUtils.applyTimeDependentCoefficientsToOuterProduct(
                    kernelPlan, time, preOrderPartials, offset, postOrderPartials, offset, scale, eigenBasisAccum);
        } else {
            AdjointMethods.accumulateEigenBasisGradientForOuterProduct(
                    eigenValues, time, preOrderPartials, offset, postOrderPartials, offset,
                    scale, eigenBasisAccum, stateCount);
        }
    }

    private void getRotatedPartial(int partialIndex,
                                   double[] rotation,
                                   double[] out,
                                   double[] intermediate) {

        // TODO does not work for STATES
        beagle.getPartials(partialIndex, Beagle.NONE, intermediate);

        int offset = 0;
        for (int i = 0; i < categoryCount; ++i) {
            for (int j = 0; j < patternCount; ++j) {

                multiplyMatrixVector(rotation, intermediate, offset, out, offset, stateCount);

                offset += stateCount;
            }
        }
    }

    @Override
    protected String getGradientTraitName() {
        return GRADIENT_TRAIT_NAME + "." + name;
    }

    public static String getName(String name) {
        return GRADIENT_TRAIT_NAME + "." + name;
    }

    @Override
    protected void constructTraits(Helper treeTraitHelper) {

        treeTraitHelper.addTrait(new TreeTrait.DA() {
            @Override
            public String getTraitName() {
                return getGradientTraitName();
            }

            @Override
            public Intent getIntent() {
                return Intent.WHOLE_TREE;
            }

            @Override
            public double[] getTrait(Tree tree, NodeRef node) {
                return getGradient(node);
            }
        });
    }

    /**
     * out[modelOffset .. modelOffset+K*K) = R^{-T} * eigenBasisAccum * R^T
     *
     * Step 1: mid = R^{-T} * accum -> mid[i,col] = sum_k ievc[k*K+i] * accum[k*K+col]
     * Step 2: out = mid * R^T      -> out[row,col] = sum_j mid[row*K+j] * evec[col*K+j]
     */
    private void rotateIntoOutput(double[] accum, double[] evec, double[] ievc,
                                   double[] out, int modelOffset) {
        final int K = stateCount;

        for (int row = 0; row < K; row++) {
            final int rowOff = row * K;
            Arrays.fill(midBuffer, rowOff, rowOff + K, 0.0);
            for (int k = 0; k < K; k++) {
                final double a = ievc[k * K + row];
                final int accumOff = k * K;
                for (int col = 0; col < K; col++) {
                    midBuffer[rowOff + col] += a * accum[accumOff + col];
                }
            }
        }

        for (int row = 0; row < K; row++) {
            final int rowOff = row * K;
            final int outRowOff = modelOffset + rowOff;
            for (int col = 0; col < K; col++) {
                double sum = 0.0;
                final int colOff = col * K;
                for (int j = 0; j < K; j++) {
                    sum += midBuffer[rowOff + j] * evec[colOff + j];
                }
                out[outRowOff + col] = sum;
            }
        }
    }

    private static void multiplyMatrixVector(double[] matrix, double[] vector, int vectorOffset, double[] out, int outOffset, int dim) {
        for (int i = 0; i < dim; i++) {
            double sum = 0.0;
            final int rowBase = i * dim;
            for (int j = 0; j < dim; j++) {
                sum += matrix[rowBase + j] * vector[vectorOffset + j];
            }
            out[outOffset + i] = sum;
        }
    }
}
