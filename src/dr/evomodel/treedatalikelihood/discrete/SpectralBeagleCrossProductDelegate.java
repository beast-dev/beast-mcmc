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
import dr.evomodel.treedatalikelihood.preorder.DiscretePartialsType;

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

    public SpectralBeagleCrossProductDelegate(String name,
                                              Tree tree,
                                              BeagleDataLikelihoodDelegate likelihoodDelegate,
                                              BranchRateModel branchRateModel,
                                              int stateCount) {

        super(name, tree, likelihoodDelegate);
        this.name = name;
        this.tree = tree;
        this.stateCount = stateCount;
        this.branchRateModel = branchRateModel;
        this.branchModel = likelihoodDelegate.getBranchModel();
        this.substitutionModelCount = branchModel.getSubstitutionModels().size();

        this.tmp1 = new double[categoryCount * patternCount * stateCount];
        this.tmp2 = new double[categoryCount * patternCount * stateCount];
        this.tmp3 = new double[categoryCount * patternCount * stateCount];
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

    private int coverWholeTree(int[] postBufferIndices,
                               int[] preBufferIndices,
                               double[] branchLengths) {
        int u = 0;
        for (int nodeNum = 0; nodeNum < tree.getNodeCount(); nodeNum++) {
            NodeRef node = tree.getNode(nodeNum);
            if (!tree.isRoot(tree.getNode(nodeNum))) {
                postBufferIndices[u] = getPostOrderPartialIndex(nodeNum);
                preBufferIndices[u]  = getPreOrderPartialIndex(nodeNum);
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

        final int[] postBufferIndices = new int[tree.getNodeCount() - 1];
        final int[] preBufferIndices = new int[tree.getNodeCount() - 1];
        final double[] branchLengths = new double[tree.getNodeCount() - 1];

        List<SubstitutionModel> substitutionModels = likelihoodDelegate.getBranchModel().getSubstitutionModels();

        if (substitutionModelCount == 1) {

            EigenDecomposition ed = substitutionModels.get(0).getEigenDecomposition();
            EigenDecomposition ted = ed.transpose();

            Arrays.fill(first, 0, first.length, 0.0);

            int count = coverWholeTree(postBufferIndices, preBufferIndices, branchLengths);
            calculateCrossProductDifferentials(postBufferIndices, preBufferIndices,
                    ed, ted,
                    branchLengths,
                    count,
                    first, null);
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

    private void calculateCrossProductDifferentials(int[] postBufferIndices,
                                                    int[] preBufferIndices,
                                                    EigenDecomposition ed,
                                                    EigenDecomposition ted,
                                                    double[] branchLengths,
                                                    int count,
                                                    double[] first,
                                                    double[] second) {
        for (int i = 0; i < count; ++i) {
            calculateCrossProductDifferentials(postBufferIndices[i], preBufferIndices[i],
                    ed, ted, branchLengths[i], first, second);
        }
    }

    private final double[] tmp1;
    private final double[] tmp2;
    private final double[] tmp3;

    private void calculateCrossProductDifferentials(int postBufferIndex,
                                                    int preBufferIndex,
                                                    EigenDecomposition ed,
                                                    EigenDecomposition ted,
                                                    double branchLength,
                                                    double[] first,
                                                    double[] second) {

        double[] postOrderPartials = tmp1;
        double[] preOrderPartials = tmp2;
        double[] intermediate = tmp3;

        getRotatedPartial(postBufferIndex,
                ed.getInverseEigenVectors(),
                postOrderPartials, intermediate);
        getRotatedPartial(preBufferIndex,
                ted.getInverseEigenVectors(),
                preOrderPartials, intermediate);

        double likelihood = blockDiagonalWeightedInnerProduct(
                preOrderPartials, 0,
                postOrderPartials, 0,
                ed.getEigenValues(), branchLength,
                stateCount);

        System.err.println("like = " + Math.log(likelihood));

    }

    private double blockDiagonalWeightedInnerProduct(double[] lhs, final int offsetLhs,
                                                     double[] rhs, final int offsetRhs,
                                                     double[] eval,
                                                     final double t,
                                                     final int dim) {
        double sum = 0.0;

        for (int i = 0; i < dim; ++i) {
            double a = eval[i];
            double expat = Math.exp(a * t);
            double b = eval[dim + i];
            if (b == 0.0) {
                sum += lhs[offsetLhs + i] * expat * rhs[offsetRhs + i];
            } else {
                double expatcosbt = expat * Math.cos(b * t);
                double expatsinbt = expat * Math.sin(b * t);

                double x0 = expatcosbt * rhs[offsetRhs + i] + expatsinbt * rhs[offsetRhs + i + 1];
                double x1 = expatcosbt * rhs[offsetRhs + i + 1] - expatsinbt * rhs[offsetRhs + 1];
                sum += lhs[offsetLhs + i] * x0 + lhs[offsetLhs + i + 1] * x1;

                ++i;
            }
        }

        return sum;
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
