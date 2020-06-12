/*
 * MultivariateTraitDebugUtilities.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.BranchRates;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.math.KroneckerOperation;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;

import java.util.HashSet;
import java.util.Set;

/**
 * @author Marc A. Suchard
 */
public class MultivariateTraitDebugUtilities {

    public static double getLengthToRoot(final Tree tree, final BranchRates branchRates, final NodeRef nodeRef) {

        double length = 0;

        if (!tree.isRoot(nodeRef)) {
            NodeRef parent = tree.getParent(nodeRef);
            double rateScalar = 1.0;
            if (branchRates != null) {
                rateScalar = branchRates.getBranchRate(tree, nodeRef);
            }
            length += rateScalar * tree.getBranchLength(nodeRef) + getLengthToRoot(tree, branchRates, parent);
        }

        return length;
    }

    private static NodeRef findMRCA(final Tree tree, final int iTip, final int jTip) {
        Set<String> leafNames = new HashSet<String>();
        leafNames.add(tree.getTaxonId(iTip));
        leafNames.add(tree.getTaxonId(jTip));
        return TreeUtils.getCommonAncestorNode(tree, leafNames);
    }

    public static void insertPrecision(Tree tree, NodeRef parent, NodeRef child, double[][] matrix, double normalization) {
        matrix[parent.getNumber()][child.getNumber()] =
                matrix[child.getNumber()][parent.getNumber()] = -1.0 / (tree.getBranchLength(child) * normalization);
        recurseGraph(tree, child, matrix, normalization);

    }

    public static void recurseGraph(Tree tree, NodeRef node, double[][] matrix, double normalization) {
        if (!tree.isExternal(node)) {
            insertPrecision(tree, node, tree.getChild(node, 0), matrix, normalization);
            insertPrecision(tree, node, tree.getChild(node, 1), matrix, normalization);
        }
    }

    public static double[][] getGraphVariance(final Tree tree, final BranchRates branchRates, final double normalization, final double priorSampleSize) {

        final int nodeCount = tree.getNodeCount();
//        final double[][] variance = new double[nodeCount + 1][nodeCount + 1];
        final double[][] variance = new double[nodeCount][nodeCount];

        for (int i = 0; i < nodeCount; i++) {

            // Fill in diagonal
            double marginalTime = getLengthToRoot(tree, branchRates, tree.getNode(i)) * normalization;
            variance[i][i] = marginalTime;

            // Fill in upper right triangle,
            for (int j = i + 1; j < nodeCount; j++) {
                NodeRef mrca = TreeUtils.getCommonAncestorSafely(tree, tree.getNode(i), tree.getNode(j));
                variance[i][j] = getLengthToRoot(tree, branchRates, mrca) * normalization;
            }
        }

        makeSymmetric(variance);

        addPrior(variance, priorSampleSize);

        return variance;
    }

    @Deprecated
    public static double[][] getTreeVarianceOld(final Tree tree, final BranchRates branchRates, final double normalization, final double priorSampleSize) {
        final int tipCount = tree.getExternalNodeCount();

        double[][] variance = new double[tipCount][tipCount];

        for (int i = 0; i < tipCount; i++) {

            // Fill in diagonal
            double marginalTime = getLengthToRoot(tree, branchRates, tree.getExternalNode(i)) * normalization;
            variance[i][i] = marginalTime;

            // Fill in upper right triangle,

            for (int j = i + 1; j < tipCount; j++) {
                NodeRef mrca = findMRCA(tree, i, j);
                variance[i][j] = getLengthToRoot(tree, branchRates, mrca) * normalization;
            }
        }

        makeSymmetric(variance);

        addPrior(variance, priorSampleSize);

        return variance;
    }

    public static double[][] getTreeVariance(final Tree tree, final BranchRates branchRates, final double normalization, final double priorSampleSize) {

        final int tipCount = tree.getExternalNodeCount();

        double[][] variance = new double[tipCount][tipCount];

        NodeRef rootNode = tree.getRoot();

        recursiveTreeVariance(variance, rootNode, tree, branchRates, normalization, priorSampleSize);

        makeSymmetric(variance);
        addPrior(variance, priorSampleSize);

        return variance;

    }


    public static PostOrderBranchStats recursiveTreeVariance(double[][] treeVariance, NodeRef node, final Tree tree, final BranchRates branchRates, final double normalization, final double priorSampleSize) {
        if (tree.isExternal(node)) {
            return new PostOrderBranchStats(new int[]{node.getNumber()}, getScaledBranchLength(node, tree, branchRates, normalization));
        }

        PostOrderBranchStats pobsL = recursiveTreeVariance(treeVariance, tree.getChild(node, 0), tree, branchRates, normalization, priorSampleSize);
        PostOrderBranchStats pobsR = recursiveTreeVariance(treeVariance, tree.getChild(node, 1), tree, branchRates, normalization, priorSampleSize);

        accumulateBranchLengths(treeVariance, pobsL);
        accumulateBranchLengths(treeVariance, pobsR);


        int nL = pobsL.dims.length;
        int nR = pobsR.dims.length;
        int n = nL + nR;

        int[] dims = new int[n];
        System.arraycopy(pobsL.dims, 0, dims, 0, nL);
        System.arraycopy(pobsR.dims, 0, dims, nL, nR);

        double scaledBranchLength = getScaledBranchLength(node, tree, branchRates, normalization);

        return new PostOrderBranchStats(dims, scaledBranchLength);

    }

    private static double getScaledBranchLength(NodeRef node, Tree tree, BranchRates branchRates, double normalization) {
        double branchLength = tree.getBranchLength(node);
        if (branchRates != null) {
            branchLength *= branchRates.getBranchRate(tree, node);
        }
        branchLength *= normalization;
        return branchLength;
    }

    private static void accumulateBranchLengths(double[][] treeVariance, PostOrderBranchStats pobs) {
        for (int i : pobs.dims) {
            for (int j : pobs.dims) {
                treeVariance[i][j] += pobs.branchLength;
            }
        }
    }

    private static class PostOrderBranchStats {
        private final int[] dims;
        private final double branchLength;

        PostOrderBranchStats(int[] dims, double branchLength) {
            this.dims = dims;
            this.branchLength = branchLength;
        }

    }


    private static void makeSymmetric(final double[][] variance) {
        for (int i = 0; i < variance.length; i++) {
            for (int j = i + 1; j < variance[i].length; j++) {
                variance[j][i] = variance[i][j];
            }
        }
    }

    private static void addPrior(final double[][] variance, double priorSampleSize) {
        if (!Double.isInfinite(priorSampleSize)) {
            for (int i = 0; i < variance.length; ++i) {
                for (int j = 0; j < variance[i].length; ++j) {
                    variance[i][j] += 1.0 / priorSampleSize;
                }
            }
        }
    }

    public static double[][] getTreeDrift(Tree tree, double[] priorMean, ContinuousDiffusionIntegrator cdi, DiffusionProcessDelegate diffusion) {

        final int dim = cdi.getDimTrait();
        final double[][] drift = new double[tree.getExternalNodeCount()][dim];

        for (int tip = 0; tip < tree.getExternalNodeCount(); ++tip) {
            drift[tip] = diffusion.getAccumulativeDrift(tree.getExternalNode(tip), priorMean, cdi, dim);
        }

        return drift;
    }

    public static double[][] getGraphDrift(Tree tree, ContinuousDiffusionIntegrator cdi, DiffusionProcessDelegate diffusion) {

        final int dim = cdi.getDimTrait();
        final double[][] drift = new double[tree.getNodeCount()][dim];

        final double[] temp = new double[dim];

        for (int node = 0; node < tree.getNodeCount(); ++node) {
            drift[node] = diffusion.getAccumulativeDrift(tree.getNode(node), temp, cdi, dim);
        }

        return drift;
    }

    public static Matrix getJointVarianceFactor(final double priorSampleSize, double[][] treeVariance, double[][] treeSharedLengths,
                                                double[][] loadingsVariance, double[][] diffusionVariance,
                                                DiffusionProcessDelegate diffusionProcessDelegate,
                                                Matrix Lt) {
        if (!diffusionProcessDelegate.hasActualization()) {

            double[][] jointVariance = diffusionProcessDelegate.getJointVariance(priorSampleSize, treeVariance, treeVariance, loadingsVariance);
            Matrix loadingsFactorsVariance = new Matrix(jointVariance);
            return loadingsFactorsVariance;

        } else {

            double[][] jointVariance = diffusionProcessDelegate.getJointVariance(priorSampleSize, treeVariance, treeSharedLengths, diffusionVariance);

            Matrix jointVarianceMatrix = new Matrix(jointVariance);

            double[][] identity = KroneckerOperation.makeIdentityMatrixArray(treeSharedLengths[0].length);
            Matrix loadingsTree = new Matrix(KroneckerOperation.product(identity, Lt.toComponents()));

            Matrix loadingsFactorsVariance = null;
            try {
                loadingsFactorsVariance = loadingsTree.product(jointVarianceMatrix.product(loadingsTree.transpose()));
            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();
            }
            return loadingsFactorsVariance;
        }
    }

    private enum Accumulator {

        OFF_DIAGONAL {
            @Override
            BranchCumulant accumulate(BranchCumulant cumulant0, double length0,
                                      BranchCumulant cumulant1, double length1) {
                return new BranchCumulant(
                        cumulant0.nTaxa + cumulant1.nTaxa,
                        cumulant0.sharedLength + cumulant1.sharedLength +
                                (cumulant0.nTaxa - 1) * cumulant0.nTaxa * length0 +
                                (cumulant1.nTaxa - 1) * cumulant1.nTaxa * length1
                );
            }
        },
        DIAGONAL {
            @Override
            BranchCumulant accumulate(BranchCumulant cumulant0, double length0,
                                      BranchCumulant cumulant1, double length1) {
                return new BranchCumulant(
                        cumulant0.nTaxa + cumulant1.nTaxa,
                        cumulant0.sharedLength + cumulant1.sharedLength +
                                cumulant0.nTaxa * length0 +
                                cumulant1.nTaxa * length1
                );
            }
        };

        abstract BranchCumulant accumulate(BranchCumulant cumulant0, double length0,
                                           BranchCumulant cumulant1, double length1);

        private class BranchCumulant {

            final int nTaxa;
            final double sharedLength;

            BranchCumulant(int nTaxa, double sharedLength) {
                this.nTaxa = nTaxa;
                this.sharedLength = sharedLength;
            }
        }


        private BranchCumulant postOrderAccumulation(Tree tree,
                                                     NodeRef node,
                                                     BranchRates branchRates) {
            if (tree.isExternal(node)) {

                return new BranchCumulant(1, 0);

            } else {

                NodeRef child0 = tree.getChild(node, 0);
                NodeRef child1 = tree.getChild(node, 1);

                BranchCumulant cumulant0 = postOrderAccumulation(tree, child0, branchRates);
                BranchCumulant cumulant1 = postOrderAccumulation(tree, child1, branchRates);

                double length0 = tree.getBranchLength(child0);
                double length1 = tree.getBranchLength(child1);

                if (branchRates != null) {
                    length0 *= branchRates.getBranchRate(tree, child0);
                    length1 *= branchRates.getBranchRate(tree, child1);
                }

                return accumulate(
                        cumulant0, length0,
                        cumulant1, length1
                );
            }
        }
    }

    public static double getVarianceOffDiagonalSum(Tree tree,
                                                   BranchRates branchRates,
                                                   double normalization) {

        Accumulator.BranchCumulant cumulant = Accumulator.OFF_DIAGONAL.postOrderAccumulation(
                tree, tree.getRoot(), branchRates);
        return cumulant.sharedLength * normalization;
    }

    public static double getVarianceDiagonalSum(Tree tree,
                                                BranchRates branchRates,
                                                double normalization) {

        Accumulator.BranchCumulant cumulant = Accumulator.DIAGONAL.postOrderAccumulation(
                tree, tree.getRoot(), branchRates);
        return cumulant.sharedLength * normalization;
    }

    public static double[] getTreeDepths(final Tree tree, final BranchRates branchRates, final double normalization) {

        final int tipCount = tree.getExternalNodeCount();

        double[] depths = new double[tipCount];

        NodeRef rootNode = tree.getRoot();

        recursiveTreeDepth(depths, rootNode, tree, branchRates, normalization);

        return depths;

    }


    public static PostOrderBranchStats recursiveTreeDepth(double[] depths, NodeRef node, final Tree tree, final BranchRates branchRates, final double normalization) {
        if (tree.isExternal(node)) {
            return new PostOrderBranchStats(new int[]{node.getNumber()}, getScaledBranchLength(node, tree, branchRates, normalization));
        }

        PostOrderBranchStats pobsL = recursiveTreeDepth(depths, tree.getChild(node, 0), tree, branchRates, normalization);
        PostOrderBranchStats pobsR = recursiveTreeDepth(depths, tree.getChild(node, 1), tree, branchRates, normalization);

        accumulateBranchLengths(depths, pobsL);
        accumulateBranchLengths(depths, pobsR);


        int nL = pobsL.dims.length;
        int nR = pobsR.dims.length;
        int n = nL + nR;

        int[] dims = new int[n];
        System.arraycopy(pobsL.dims, 0, dims, 0, nL);
        System.arraycopy(pobsR.dims, 0, dims, nL, nR);

        double scaledBranchLength = getScaledBranchLength(node, tree, branchRates, normalization);

        return new PostOrderBranchStats(dims, scaledBranchLength);

    }

    private static void accumulateBranchLengths(double[] depths, PostOrderBranchStats pobs) {
        for (int i : pobs.dims) {
            depths[i] += pobs.branchLength;
        }
    }
}
