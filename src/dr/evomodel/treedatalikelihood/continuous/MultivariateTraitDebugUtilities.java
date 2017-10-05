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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;

import java.util.HashSet;
import java.util.Set;

/**
 * @author Marc A. Suchard
 */
public class MultivariateTraitDebugUtilities {

    public static double getLengthToRoot(final Tree tree, final NodeRef nodeRef) {

        double length = 0;

        if (!tree.isRoot(nodeRef)) {
            NodeRef parent = tree.getParent(nodeRef);
            length += tree.getBranchLength(nodeRef) + getLengthToRoot(tree, parent);
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
                           matrix[child.getNumber()][parent.getNumber()] =-1.0 / (tree.getBranchLength(child) * normalization);
        recurseGraph(tree, child, matrix, normalization);

    }

    public static void recurseGraph(Tree tree, NodeRef node, double[][] matrix, double normalization) {
        if (!tree.isExternal(node)) {
            insertPrecision(tree, node, tree.getChild(node, 0), matrix, normalization);
            insertPrecision(tree, node, tree.getChild(node, 1), matrix, normalization);
        }
    }

    public static double[][] getGraphPrecision(final Tree tree, final double normalization) {

        final int nodeCount = tree.getNodeCount();
        final double[][] precision = new double[nodeCount][nodeCount];

        recurseGraph(tree, tree.getRoot(), precision, normalization);

//        for (int i = 0; i < nodeCount; i++) {
//
//            // Fill in diagonal
//            double marginalTime = getLengthToRoot(tree, tree.getNode(i)) * normalization;
//            variance[i][i] = marginalTime;
//
//            // Fill in upper right triangle,
//            for (int j = i + 1; j < nodeCount; j++) {
//
//                NodeRef mrca = TreeUtils.getCommonAncestor(tree, tree.getNode(i), tree.getNode(j));
//                variance[i][j] = getLengthToRoot(tree, mrca);
//            }
//        }

//        makeSymmetric(variance);
//        addPrior(variance, priorSampleSize);

        return precision;
    }

    public static double[][] getGraphVariance(final Tree tree, final double normalization, final double priorSampleSize) {

        final int nodeCount = tree.getNodeCount();
//        final double[][] variance = new double[nodeCount + 1][nodeCount + 1];
        final double[][] variance = new double[nodeCount][nodeCount];

        for (int i = 0; i < nodeCount; i++) {

            // Fill in diagonal
            double marginalTime = getLengthToRoot(tree, tree.getNode(i)) * normalization;
            variance[i][i] = marginalTime;

            // Fill in upper right triangle,
            for (int j = i + 1; j < nodeCount; j++) {
                NodeRef mrca = TreeUtils.getCommonAncestor(tree, tree.getNode(i), tree.getNode(j));
                variance[i][j] = getLengthToRoot(tree, mrca) * normalization;
            }
        }

        makeSymmetric(variance);

        addPrior(variance, priorSampleSize);
//        variance[nodeCount][nodeCount] = 0.0; // Does not matter, since prior mean is "observed"
//        
//        for (int i = 0; i < nodeCount + 1; ++i) {
//            variance[i][i] += 10.0;
//        }

        return variance;
    }

    public static double[][] getTreeVariance(final Tree tree, final double normalization, final double priorSampleSize) {

        final int tipCount = tree.getExternalNodeCount();
        int length = tipCount;

        double[][] variance = new double[length][length];

        for (int i = 0; i < tipCount; i++) {

            // Fill in diagonal
            double marginalTime = getLengthToRoot(tree, tree.getExternalNode(i)) * normalization;
            variance[i][i] = marginalTime;

            // Fill in upper right triangle,

            for (int j = i + 1; j < tipCount; j++) {
                NodeRef mrca = findMRCA(tree, i, j);
                variance[i][j] = getLengthToRoot(tree, mrca) * normalization;
            }
        }

        // Make symmetric
//        for (int i = 0; i < length; i++) {
//            for (int j = i + 1; j < length; j++) {
//                variance[j][i] = variance[i][j];
//            }
//        }
        makeSymmetric(variance);

//        if (!Double.isInfinite(priorSampleSize)) {
//            for (int i = 0; i < variance.length; ++i) {
//                for (int j = 0; j < variance[i].length; ++j) {
//                    variance[i][j] += 1.0 / priorSampleSize;
//                }
//            }
//        }
        addPrior(variance, priorSampleSize);

        return variance;
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

//    public String getDebugInformation(final FullyConjugateMultivariateTraitLikelihood traitLikelihood) {
//
//        StringBuilder sb = new StringBuilder();
////        sb.append(this.g)
////        System.err.println("Hello");
//        sb.append("Tree:\n");
//        sb.append(traitLikelihood.getId()).append("\t");
//
//        final MutableTreeModel treeModel = traitLikelihood.getTreeModel();
//        sb.append(treeModel.toString());
//        sb.append("\n\n");
//
//        double[][] treeVariance = computeTreeVariance(true);
//        double[][] traitPrecision = traitLikelihood.getDiffusionModel().getPrecisionmatrix();
//        Matrix traitVariance = new Matrix(traitPrecision).inverse();
//
//        double[][] jointVariance = KroneckerOperation.product(treeVariance, traitVariance.toComponents());
//
//        sb.append("Tree variance:\n");
//        sb.append(new Matrix(treeVariance));
//        sb.append(matrixMin(treeVariance)).append("\t").append(matrixMax(treeVariance)).append("\t").append(matrixSum(treeVariance));
//        sb.append("\n\n");
//        sb.append("Trait variance:\n");
//        sb.append(traitVariance);
//        sb.append("\n\n");
////        sb.append("Joint variance:\n");
////        sb.append(new Matrix(jointVariance));
////        sb.append("\n\n");
//
//        sb.append("Tree dim: " + treeVariance.length + "\n");
//        sb.append("data dim: " + jointVariance.length);
//        sb.append("\n\n");
//
//        double[] data = new double[jointVariance.length];
//        System.arraycopy(meanCache, 0, data, 0, jointVariance.length);
//
//        if (nodeToClampMap != null) {
//            int offset = treeModel.getExternalNodeCount() * getDimTrait();
//            for(Map.Entry<NodeRef, RestrictedPartials> clamps : nodeToClampMap.entrySet()) {
//                double[] partials = clamps.getValue().getPartials();
//                for (int i = 0; i < partials.length; ++i) {
//                    data[offset] = partials[i];
//                    ++offset;
//                }
//            }
//        }
//
//        sb.append("Data:\n");
//        sb.append(new Vector(data)).append("\n");
//        sb.append(data.length).append("\t").append(vectorMin(data)).append("\t").append(vectorMax(data)).append("\t").append(vectorSum(data));
//        sb.append(treeModel.getNodeTaxon(treeModel.getExternalNode(0)).getId());
//        sb.append("\n\n");
//
//        MultivariateNormalDistribution mvn = new MultivariateNormalDistribution(new double[data.length], new Matrix(jointVariance).inverse().toComponents());
//        double logDensity = mvn.logPdf(data);
//        sb.append("logLikelihood: " + getLogLikelihood() + " == " + logDensity + "\n\n");
//
//        final WishartSufficientStatistics sufficientStatistics = getWishartStatistics();
//        final double[] outerProducts = sufficientStatistics.getScaleMatrix();
//
//        sb.append("Outer-products (DP):\n");
//        sb.append(new Vector(outerProducts));
//        sb.append(sufficientStatistics.getDf() + "\n");
//
//        Matrix treePrecision = new Matrix(treeVariance).inverse();
//        final int n = data.length / traitPrecision.length;
//        final int p = traitPrecision.length;
//        double[][] tmp = new double[n][p];
//
//        for (int i = 0; i < n; ++i) {
//            for (int j = 0; j < p; ++j) {
//                tmp[i][j] = data[i * p + j];
//            }
//        }
//        Matrix y = new Matrix(tmp);
//
//        Matrix S = null;
//        try {
//            S = y.transpose().product(treePrecision).product(y); // Using Matrix-Normal form
//        } catch (IllegalDimension illegalDimension) {
//            illegalDimension.printStackTrace();
//        }
//        sb.append("Outer-products (from tree variance:\n");
//        sb.append(S);
//        sb.append("\n\n");
//
//        return sb.toString();
//    }

    public static double[][] getTreeDrift(Tree tree, DiffusionProcessDelegate diffusion) {

        final int dim = diffusion.getDiffusionModel(0).getPrecisionParameter().getColumnDimension();
        final double[][] drift = new double[tree.getExternalNodeCount()][dim];

        if (diffusion instanceof DriftDiffusionModelDelegate) {
            for (int tip = 0; tip < tree.getExternalNodeCount(); ++tip) {
                drift[tip] = ((DriftDiffusionModelDelegate) diffusion).getAccumulativeDrift(
                        tree.getExternalNode(tip));
            }
        }

        return drift;
    }

    public static double[][] getGraphDrift(Tree tree, DiffusionProcessDelegate diffusion) {

        final int dim = diffusion.getDiffusionModel(0).getPrecisionParameter().getColumnDimension();
        final double[][] drift = new double[tree.getNodeCount()][dim];

        if (diffusion instanceof DriftDiffusionModelDelegate) {
            for (int node = 0; node < tree.getNodeCount(); ++node) {
                drift[node] = ((DriftDiffusionModelDelegate) diffusion).getAccumulativeDrift(
                        tree.getNode(node));
            }
        }

        return drift;
    }
}
