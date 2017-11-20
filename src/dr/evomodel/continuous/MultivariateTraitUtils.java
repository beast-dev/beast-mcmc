/*
 * MultivariateTraitUtils.java
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

package dr.evomodel.continuous;

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeUtils;
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
        MutableTreeModel treeModel = trait.getTreeModel();
        Set<String> leafNames = new HashSet<String>();
        leafNames.add(treeModel.getTaxonId(iTip));
        leafNames.add(treeModel.getTaxonId(jTip));
        return TreeUtils.getCommonAncestorNode(treeModel, leafNames);
    }

    public static double[][] computeTreePrecision(FullyConjugateMultivariateTraitLikelihood trait, boolean conditionOnRoot) {
        if (trait.strengthOfSelection != null) {
            return new SymmetricMatrix(computeTreeVarianceOU(trait, conditionOnRoot)).inverse().toComponents();
        } else {
            return new SymmetricMatrix(computeTreeVariance(trait, conditionOnRoot)).inverse().toComponents();
        }
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

    private static double[][] productMatrices(double[][] A, double[][] B) {
        double[][] C = new double[A.length][B[0].length];


        for (int i = 0; i < A.length; i++) {
            for (int k = 0; k < B[0].length; k++) {
                for (int j = 0; j < A[0].length; j++) {
                    C[i][k] = C[i][k] + A[i][j] * B[j][k];
                }
            }
        }

        return C;
    }

    private static double[][] transposeMatrix(double[][] A) {
        double[][] B = new double[A[0].length][A.length];

        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[0].length; j++) {
                B[j][i] = A[i][j];
            }
        }

        return B;
    }


    private static double[][] computeLinCombMatrix(FullyConjugateMultivariateTraitLikelihood trait) {

        MutableTreeModel treeModel = trait.getTreeModel();
        final int tipCount = treeModel.getExternalNodeCount();
        final int branchCount = 2 * tipCount - 2;
        double[][] linCombMatrix = new double[tipCount][branchCount];

        double tempScalar;
        NodeRef tempNode;
        int tempNodeIndex;

        for (int k = 0; k < tipCount; k++) {
            tempNode = treeModel.getExternalNode(k);
            //check if treeModel.getExternalNode(k).getNumber() == k
            tempScalar = 1;
            tempNodeIndex = k;

            for (int r = 0; r < branchCount; r++) {
                if (r == tempNodeIndex) {
                    linCombMatrix[k][r] = tempScalar;
                    // tempScalar = tempScalar * (1 - treeModel.getBranchLength(tempNode) * trait.getTimeScaledSelection(tempNode));
                    tempScalar = tempScalar * Math.exp(-trait.getTimeScaledSelection(tempNode));
                    tempNode = treeModel.getParent(tempNode);
                    tempNodeIndex = tempNode.getNumber();
                } else {
                    linCombMatrix[k][r] = 0;
                }

            }
        }

        return linCombMatrix;

    }


    private static double[] computeRootMultipliers(FullyConjugateMultivariateTraitLikelihood trait) {
        MutableTreeModel myTreeModel = trait.getTreeModel();
        final int tipCount = myTreeModel.getExternalNodeCount();
        double[] multiplierVect = new double[tipCount];
        NodeRef tempNode;

        for (int i = 0; i < tipCount; i++) {
            tempNode = myTreeModel.getExternalNode(i);
            multiplierVect[i] = Math.exp(-trait.getTimeScaledSelection(tempNode));
            tempNode = myTreeModel.getParent(tempNode);
            while (!myTreeModel.isRoot(tempNode)) {
                multiplierVect[i] = multiplierVect[i] * Math.exp(-trait.getTimeScaledSelection(tempNode));
                tempNode = myTreeModel.getParent(tempNode);
            }
        }


        return multiplierVect;

    }


    private static double[] getShiftContributionToMean(NodeRef node, FullyConjugateMultivariateTraitLikelihood trait) {

        MutableTreeModel treeModel = trait.getTreeModel();
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


    public static double[] computeTreeTraitMean(FullyConjugateMultivariateTraitLikelihood trait, double[] rootValue, boolean conditionOnRoot) {
        double[] root = trait.getPriorMean();
        if (conditionOnRoot) {
            System.err.println("WARNING: Not yet fully implemented (conditioning on root in simulator)");
            //root = new double[root.length];
            root = rootValue;
        }
        final int nTaxa = trait.getTreeModel().getExternalNodeCount();
        double[] mean = new double[root.length * nTaxa];
        for (int i = 0; i < nTaxa; ++i) {
            System.arraycopy(root, 0, mean, i * root.length, root.length);
        }

        if (trait.driftModels != null) {
            MutableTreeModel myTreeModel = trait.getTreeModel();
            for (int i = 0; i < nTaxa; ++i) {
                double[] shiftContribution = getShiftContributionToMean(myTreeModel.getExternalNode(i), trait);
                for (int j = 0; j < trait.dimTrait; ++j) {
                    mean[i * trait.dimTrait + j] = mean[i * trait.dimTrait + j] + shiftContribution[j];
                }
            }
        }

        return mean;
    }

    public static double[] computeTreeTraitMeanOU(FullyConjugateMultivariateTraitLikelihood trait, double[] rootValue, boolean conditionOnRoot) {

        double[] root = trait.getPriorMean();
        MutableTreeModel myTreeModel = trait.getTreeModel();
        double[][] linCombMatrix = computeLinCombMatrix(trait);
        double[] rootMultiplierVect = computeRootMultipliers(trait);

        if (conditionOnRoot) {
            root = rootValue;
        }

        final int nTaxa = myTreeModel.getExternalNodeCount();
        final int branchCount = 2 * nTaxa - 2;
        final int traitDim = trait.dimTrait;
        double[] mean = new double[root.length * nTaxa];
        double[] displacementMeans = new double[branchCount * traitDim];
        double[] tempVect = new double[nTaxa * traitDim];
        NodeRef tempNode;

        for (int k = 0; k < branchCount; k++) {
            tempNode = myTreeModel.getNode(k);
            for (int t = 0; t < traitDim; t++) {
                displacementMeans[k * traitDim + t] = (1 - Math.exp(-trait.getTimeScaledSelection(tempNode))) * trait.getOptimalValue(tempNode)[t];
            }
        }

        //check this
        //multiply linCombMatrix with displacement means
        for (int i = 0; i < nTaxa; i++) {
            for (int j = 0; j < branchCount; j++) {
                for (int k = 0; k < traitDim; k++) {
                    tempVect[i * traitDim + k] = tempVect[i * traitDim + k] + linCombMatrix[i][j] * displacementMeans[j * traitDim + k];
                }
            }

        }

        for (int i = 0; i < nTaxa; ++i) {
            System.arraycopy(root, 0, mean, i * root.length, root.length);

            for (int j = 0; j < traitDim; ++j) {
                mean[i * traitDim + j] = mean[i * traitDim + j] * rootMultiplierVect[i] + tempVect[i * traitDim + j];
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
        MutableTreeModel treeModel = trait.getTreeModel();
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

    public static double[][] computeTreeVarianceOU(FullyConjugateMultivariateTraitLikelihood trait, boolean conditionOnRoot) {
        MutableTreeModel treeModel = trait.getTreeModel();
        final int tipCount = treeModel.getExternalNodeCount();
        final int branchCount = 2 * tipCount - 2;
        double[][] variance = new double[tipCount][tipCount];
        double[][] tempMatrix = new double[tipCount][branchCount];
        double[][] diagMatrix = new double[branchCount][branchCount];


        /*
        double tempScalar;
        NodeRef tempNode;
        int tempNodeIndex;

        for(int k = 0; k < tipCount; k++){
           tempNode = treeModel.getExternalNode(k);
            //check if treeModel.getExternalNode(k).getNumber() == k
           tempScalar = 1;
           tempNodeIndex = k;

           for(int r = 0; r < branchCount; r++){
               if(r == tempNodeIndex){
                   tempMatrix[k][r] = tempScalar;
                   tempScalar = tempScalar*(1-treeModel.getBranchLength(tempNode)*trait.getTimeScaledSelection(tempNode));
                   tempNode = treeModel.getParent(tempNode);
                   tempNodeIndex = tempNode.getNumber();
               }else{
                   tempMatrix[k][r]= 0;
               }

           }
        }
        */

        tempMatrix = computeLinCombMatrix(trait);

        for (int i = 0; i < branchCount; i++) {
            //  diagMatrix[i][i] = treeModel.getBranchLength(treeModel.getNode(i));
            // diagMatrix[i][i] = (2*strengthOfSelection.getBranchRate(treeModel, node) ) / (1 - Math.exp(-2*getTimeScaledSelection(node)));
            diagMatrix[i][i] = (1 - Math.exp(-2 * trait.getTimeScaledSelection(treeModel.getNode(i))))
                    / (2 * trait.strengthOfSelection.getBranchRate(treeModel, treeModel.getNode(i)));
        }


        variance = productMatrices(productMatrices(tempMatrix, diagMatrix), transposeMatrix(tempMatrix));

        return variance;
    }


    private static final boolean DEBUG = false;

}
