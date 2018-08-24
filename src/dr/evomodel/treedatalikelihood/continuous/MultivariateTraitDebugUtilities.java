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
                NodeRef mrca = TreeUtils.getCommonAncestorSafely(tree, tree.getNode(i), tree.getNode(j));
                variance[i][j] = getLengthToRoot(tree, mrca) * normalization;
            }
        }

        makeSymmetric(variance);

        addPrior(variance, priorSampleSize);

        return variance;
    }

    public static double[][] getTreeVariance(final Tree tree, final double normalization, final double priorSampleSize) {

        final int tipCount = tree.getExternalNodeCount();

        double[][] variance = new double[tipCount][tipCount];

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

        makeSymmetric(variance);

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
