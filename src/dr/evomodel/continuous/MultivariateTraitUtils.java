/*
 * MultivariateTraitUtils.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.continuous;

import dr.evolution.tree.MultivariateTraitTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.math.KroneckerOperation;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.SymmetricMatrix;

import java.util.HashSet;
import java.util.Set;

/**
 * @author Marc A. Suchard
 */
public class MultivariateTraitUtils {

    public static NodeRef findMRCA(FullyConjugateMultivariateTraitLikelihood trait, int iTip, int jTip) {
        MultivariateTraitTree treeModel = trait.getTreeModel();
        Set<String> leafNames = new HashSet<String>();
        leafNames.add(treeModel.getTaxonId(iTip));
        leafNames.add(treeModel.getTaxonId(jTip));
        return Tree.Utils.getCommonAncestorNode(treeModel, leafNames);
    }

    public static double[][] computeTreePrecision(FullyConjugateMultivariateTraitLikelihood trait, boolean conditionOnRoot) {
        return new SymmetricMatrix(computeTreeVariance(trait, conditionOnRoot)).inverse().toComponents();
    }

    public static double[][] computeTreeTraitPrecision(FullyConjugateMultivariateTraitLikelihood trait, boolean conditionOnRoot) {
        double[][] treePrecision = computeTreePrecision(trait, conditionOnRoot);
        double[][] traitPrecision = trait.getDiffusionModel().getPrecisionmatrix();
        return productKronecker(treePrecision, traitPrecision);
    }

    private static double[][] productKronecker(double[][] A, double[][] B) {
        if (B.length > 1) {
            A = KroneckerOperation.product(A, B);
        } else {
            final double b = B[0][0];
            for (int i = 0; i < A.length; ++i) {
                for (int j = 0; j < A[i].length; ++j) {
                    A[i][j] *= b;
                }
            }
        }
        return A;
    }

    private static double[] getShiftContributionToMean(NodeRef node, FullyConjugateMultivariateTraitLikelihood trait) {

        MultivariateTraitTree treeModel = trait.getTreeModel();
        double shiftContribution[] = new double[trait.dimTrait];

        if (!treeModel.isRoot(node)) {
            NodeRef parent = treeModel.getParent(node);
            double shiftContributionParent[] = getShiftContributionToMean(parent, trait);
            for (int i = 0; i < shiftContribution.length; ++i) {
                shiftContribution[i] = trait.getShiftForBranchLength(node)[i] + shiftContributionParent[i];
            }
        }
        return shiftContribution;
    }


    public static double[] computeTreeTraitMean(FullyConjugateMultivariateTraitLikelihood trait, boolean conditionOnRoot) {
        double[] root = trait.getPriorMean();
        if (conditionOnRoot) {
            System.err.println("WARNING: Not yet fully implemented (conditioning on root in simulator)");
            root = new double[root.length];
        }
        final int nTaxa = trait.getTreeModel().getExternalNodeCount();
        double[] mean = new double[root.length * nTaxa];
        for (int i = 0; i < nTaxa; ++i) {
            System.arraycopy(root, 0, mean, i * root.length, root.length);
        }

        if (trait.driftModels != null) {
            MultivariateTraitTree myTreeModel = trait.getTreeModel();
            for (int i = 0; i < nTaxa; ++i) {
                double[] shiftContribution = getShiftContributionToMean(myTreeModel.getExternalNode(i), trait);
                for (int j = 0; j < trait.dimTrait; ++j) {
                    mean[i * trait.dimTrait + j] = mean[i * trait.dimTrait + j] + shiftContribution[j];
                }
            }
        }

        return mean;
    }

    public static double[][] computeTreeTraitVariance(FullyConjugateMultivariateTraitLikelihood trait, boolean conditionOnRoot) {
        double[][] treeVariance = computeTreeVariance(trait, conditionOnRoot);
        double[][] traitVariance =
                new SymmetricMatrix(trait.getDiffusionModel().getPrecisionmatrix()).inverse().toComponents();
        return productKronecker(treeVariance, traitVariance);
    }

    public static double[][] computeTreeVariance(FullyConjugateMultivariateTraitLikelihood trait, boolean conditionOnRoot) {
        MultivariateTraitTree treeModel = trait.getTreeModel();
        final int tipCount = treeModel.getExternalNodeCount();
        double[][] variance = new double[tipCount][tipCount];

        for (int i = 0; i < tipCount; i++) {

            // Fill in diagonal
            double marginalTime = trait.getRescaledLengthToRoot(treeModel.getExternalNode(i));
            variance[i][i] = marginalTime;

            // Fill in upper right triangle,

            for (int j = i + 1; j < tipCount; j++) {
                NodeRef mrca = findMRCA(trait, i, j);
                if (DEBUG) {
                    System.err.println(trait.getTreeModel().getRoot().getNumber());
                    System.err.print("Taxa pair: " + i + " : " + j + " (" + mrca.getNumber() + ") = ");
                }
                double length = trait.getRescaledLengthToRoot(mrca);
                if (DEBUG) {
                    System.err.println(length);
                }
                variance[i][j] = length;
            }
        }

        // Make symmetric
        for (int i = 0; i < tipCount; i++) {
            for (int j = i + 1; j < tipCount; j++) {
                variance[j][i] = variance[i][j];
            }
        }

        if (DEBUG) {
            System.err.println("");
            System.err.println("New tree conditional variance:\n" + new Matrix(variance));
        }

        // TODO Handle missing values
//        variance = removeMissingTipsInTreeVariance(variance); // Automatically prune missing tips
//
//        if (DEBUG) {
//            System.err.println("");
//            System.err.println("New tree (trimmed) conditional variance:\n" + new Matrix(variance));
//        }

        if (!conditionOnRoot) {
            double priorSampleSize = trait.getPriorSampleSize();
            for (int i = 0; i < tipCount; ++i) {
                for (int j = 0; j < tipCount; ++j) {
                    variance[i][j] += 1.0 / priorSampleSize;
                }
            }
            if (DEBUG) {
                System.err.println("");
                System.err.println("New tree unconditional variance:\n" + new Matrix(variance));
            }
        }
        return variance;
    }

    private static final boolean DEBUG = false;

}
