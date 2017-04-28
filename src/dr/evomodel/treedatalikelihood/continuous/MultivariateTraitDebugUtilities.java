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

import dr.evolution.tree.MultivariateTraitTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.TaxonList;
import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.evomodel.continuous.RestrictedPartials;
import dr.math.KroneckerOperation;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;

import java.util.HashSet;
import java.util.Map;
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
                variance[i][j] = getLengthToRoot(tree, mrca);
            }
        }

//        if (DO_CLAMP && nodeToClampMap != null) {
//            List<RestrictedPartials> partialsList = new ArrayList<RestrictedPartials>();
//            for (Map.Entry<NodeRef, RestrictedPartials> keySet : nodeToClampMap.entrySet()) {
//                partialsList.add(keySet.getValue());
//            }
//
//            for (int i = 0; i < partialsList.size(); ++i) {
//                RestrictedPartials partials = partialsList.get(i);
//                NodeRef node = partials.getNode();
//
//                variance[tipCount + i][tipCount + i] = getRescaledLengthToRoot(node) +
//                        1.0 / partials.getPriorSampleSize();
//
//                for (int j = 0; j < tipCount; ++j) {
//                    NodeRef friend = treeModel.getExternalNode(j);
//                    NodeRef mrca = Tree.Utils.getCommonAncestor(treeModel, node, friend);
//                    variance[j][tipCount + i] = getRescaledLengthToRoot(mrca);
//
//                }
//
//                for (int j = 0; j < i; ++j) {
//                    NodeRef friend = partialsList.get(j).getNode();
//                    NodeRef mrca = Tree.Utils.getCommonAncestor(treeModel, node, friend);
//                    variance[tipCount + j][tipCount + i] = getRescaledLengthToRoot(mrca);
//                }
//            }
//        }

        // Make symmetric
        for (int i = 0; i < length; i++) {
            for (int j = i + 1; j < length; j++) {
                variance[j][i] = variance[i][j];
            }
        }

        if (!Double.isInfinite(priorSampleSize)) {
            for (int i = 0; i < variance.length; ++i) {
                for (int j = 0; j < variance[i].length; ++j) {
                    variance[i][j] += 1.0 / priorSampleSize;
                }
            }
        }

        return variance;
    }

//    public String getDebugInformation(final FullyConjugateMultivariateTraitLikelihood traitLikelihood) {
//
//        StringBuilder sb = new StringBuilder();
////        sb.append(this.g)
////        System.err.println("Hello");
//        sb.append("Tree:\n");
//        sb.append(traitLikelihood.getId()).append("\t");
//
//        final MultivariateTraitTree treeModel = traitLikelihood.getTreeModel();
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
}
