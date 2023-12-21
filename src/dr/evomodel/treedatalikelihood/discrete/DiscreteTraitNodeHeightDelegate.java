/*
 * DiscreteTraitNodeHeightDelegate.java
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


package dr.evomodel.treedatalikelihood.discrete;

import beagle.Beagle;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.DifferentiableBranchRates;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;

import java.util.Arrays;

/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class DiscreteTraitNodeHeightDelegate extends DiscreteTraitBranchRateDelegate {

    static final String GRADIENT_TRAIT_NAME = "NodeHeightGradient";

    static final String HESSIAN_TRAIT_NAME = "NodeHeightHessian";

    private final DifferentiableBranchRates branchRates;

    DiscreteTraitNodeHeightDelegate(String name,
                                    Tree tree,
                                    BeagleDataLikelihoodDelegate likelihoodDelegate,
                                    DifferentiableBranchRates branchRates) {
        super(name, tree, likelihoodDelegate);
        this.branchRates = branchRates;
    }

    protected int getGradientLength() {
        return tree.getInternalNodeCount();
    }

    protected void getNodeDerivatives(Tree tree, double[] first, double[] second) {
        double[] branchGradient = new double[tree.getNodeCount() - 1];
        double[] branchDiagonalHessian = second == null ? null : new double[tree.getNodeCount() - 1];
//        double[] branchDiagonalHessian = new double[first.length];
        super.getNodeDerivatives(tree, branchGradient, branchDiagonalHessian);

        Arrays.fill(first, 0.0);

        for (int i = 0; i < tree.getInternalNodeCount(); ++i) {

            final  NodeRef node = tree.getNode(i + tree.getExternalNodeCount());

            for (int j = 0; j < tree.getChildCount(node); j++) {
                NodeRef childNode = tree.getChild(node, j);
                final int childNodeIndex = getParameterIndex(childNode, tree);
                first[i] += branchGradient[childNodeIndex] * branchRates.getBranchRate(tree, childNode);
            }
            if (!tree.isRoot(node)) {
                first[i] -= branchGradient[getParameterIndex(node, tree)] * branchRates.getBranchRate(tree, node);
            }
        }

//        double[] testBranchGradient = new double[tree.getNodeCount()];

        if (second != null) {
            final int internalNodeCount = tree.getInternalNodeCount();
            double[] sisterBranchesSecondDerivatives = new double[internalNodeCount];
            double[] currentAndParentBranchesSecondDerivatives = new double[tree.getNodeCount() - 1];

            double[][] prePartials = new double[tree.getNodeCount()][patternCount * stateCount * categoryCount];
            double[][] postPartials = new double[tree.getNodeCount()][patternCount * stateCount * categoryCount];
            double[][] transitionMatrices = new double[tree.getNodeCount()][stateCount * stateCount * categoryCount];

            for (int i = 0; i < tree.getNodeCount(); i++) {
                beagle.getPartials(getPostOrderPartialIndex(i), Beagle.NONE, postPartials[i]);
                beagle.getTransitionMatrix(evolutionaryProcessDelegate.getMatrixIndex(i), transitionMatrices[i]);
                beagle.getPartials(getPreOrderPartialIndex(i), Beagle.NONE, prePartials[i]);
            }

            double[] Qi = new double[stateCount * stateCount * categoryCount];
            double[] Qj = new double[stateCount * stateCount * categoryCount];
            double[] Qk = new double[stateCount * stateCount * categoryCount];
            double[] tmpLeftPartail = new double[patternCount * stateCount * categoryCount];
            double[] tmpRightPartial = new double[patternCount * stateCount * categoryCount];
            double[] tmpQLeftPartial = new double[patternCount * stateCount * categoryCount];
            double[] tmpQRightPartial = new double[patternCount * stateCount * categoryCount];
            double[] tmpIPartial = new double[patternCount * stateCount * categoryCount];

            double[] testQ = new double[stateCount * stateCount * categoryCount];
            double[][] denominator = new double[tree.getNodeCount()][];
            evolutionaryProcessDelegate.getSubstitutionModel(0).getInfinitesimalMatrix(testQ);

            for (int i = 0; i < internalNodeCount; i++) {
                NodeRef nodeI = tree.getNode(i + tree.getExternalNodeCount());
                NodeRef nodeJ = tree.getChild(nodeI, 0);
                NodeRef nodeK = tree.getChild(nodeI, 1);
                denominator[i] = getDoubleVectorReduction(postPartials[nodeI.getNumber()], prePartials[nodeI.getNumber()], true);

                beagle.getTransitionMatrix(evolutionaryProcessDelegate.getInfinitesimalMatrixBufferIndex(nodeI.getNumber()), Qi);
                beagle.getTransitionMatrix(evolutionaryProcessDelegate.getInfinitesimalMatrixBufferIndex(nodeJ.getNumber()), Qj);
                beagle.getTransitionMatrix(evolutionaryProcessDelegate.getInfinitesimalMatrixBufferIndex(nodeK.getNumber()), Qk);

                getMatrixVectorProduct(transitionMatrices[nodeJ.getNumber()], postPartials[nodeJ.getNumber()], tmpLeftPartail);  // Pj p_j
                getMatrixVectorProduct(Qj, tmpLeftPartail, tmpQLeftPartial); // QjPj p_j

                getMatrixVectorProduct(transitionMatrices[nodeK.getNumber()], postPartials[nodeK.getNumber()], tmpRightPartial); // Pk p_k
                getMatrixVectorProduct(Qk, tmpRightPartial, tmpQRightPartial); // Qk Pk p_k

                sisterBranchesSecondDerivatives[i] = getNormalizedPatternReduction(getTripleVectorReduction(tmpQLeftPartial, tmpQRightPartial, prePartials[nodeI.getNumber()], true), denominator[i])
                        - branchGradient[getParameterIndex(nodeJ, tree)] * branchGradient[getParameterIndex(nodeK, tree)];

                second[i] = branchDiagonalHessian[getParameterIndex(nodeJ, tree)] * branchRates.getBranchRate(tree, nodeJ) * branchRates.getBranchRate(tree, nodeJ)
                        + branchDiagonalHessian[getParameterIndex(nodeK, tree)] * branchRates.getBranchRate(tree, nodeK) * branchRates.getBranchRate(tree, nodeK)
                        + 2 * sisterBranchesSecondDerivatives[i] * branchRates.getBranchRate(tree, nodeJ) * branchRates.getBranchRate(tree, nodeK);


                if (!tree.isRoot(nodeI)) {
                    getMatrixVectorProduct(Qi, prePartials[nodeI.getNumber()], tmpIPartial);  //Qi q_i
                    currentAndParentBranchesSecondDerivatives[getParameterIndex(nodeJ, tree)] = getNormalizedPatternReduction(getTripleVectorReduction(tmpQLeftPartial, tmpRightPartial, tmpIPartial, true), denominator[i])
                            - branchGradient[getParameterIndex(nodeJ, tree)] * branchGradient[getParameterIndex(nodeI, tree)];
                    currentAndParentBranchesSecondDerivatives[getParameterIndex(nodeK, tree)] = getNormalizedPatternReduction(getTripleVectorReduction(tmpQRightPartial, tmpLeftPartail, tmpIPartial, true), denominator[i])
                            - branchGradient[getParameterIndex(nodeK, tree)] * branchGradient[getParameterIndex(nodeI, tree)];


                    second[i] += branchDiagonalHessian[getParameterIndex(nodeI, tree)] * branchRates.getBranchRate(tree, nodeI) * branchRates.getBranchRate(tree, nodeI)
                            - 2 * currentAndParentBranchesSecondDerivatives[getParameterIndex(nodeJ, tree)] * branchRates.getBranchRate(tree, nodeI) * branchRates.getBranchRate(tree, nodeJ)
                            - 2 * currentAndParentBranchesSecondDerivatives[getParameterIndex(nodeK, tree)] * branchRates.getBranchRate(tree, nodeI) * branchRates.getBranchRate(tree, nodeK);
                }

            }
//            for (int i = 0; i < tree.getNodeCount() - 1; i++) {
//                getMatrixVectorProduct(Qi, postPartials[i], tmpIPartial);
//                denominator[i] = getDoubleVectorReduction(postPartials[i], prePartials[i], true);
//                testBranchGradient[i] = getNormalizedPatternReduction(getDoubleVectorReduction(prePartials[i], tmpIPartial, true), denominator[i]);
//            }

//            int test = 0;
        }



    }

    private int getParameterIndex(NodeRef node, Tree tree) {
        return node.getNumber() < tree.getRoot().getNumber() ? node.getNumber() : node.getNumber() - 1;
    }

    private double[] getVectorVectorProduct(double[] first, double[] second) {
        double[] product = new double[first.length];
        for (int i = 0; i < first.length; i++) {
            product[i] = first[i] * second[i];
        }
        return product;
    }

    private int getPrePartialIndex (int nodeIndex, Tree tree) {
        return nodeIndex + tree.getExternalNodeCount();
    }

    private void getMatrixVectorProduct(double[] matrix, double[] vector, double[] result) {
        assert(vector.length == result.length);
        for (int category = 0; category < categoryCount; category++) {
            for (int pattern = 0; pattern < patternCount; pattern++) {
                for (int i = 0; i < stateCount; i++) {
                    double sum = 0;
                    for (int j = 0; j < stateCount; j++) {
                        sum += matrix[category * stateCount * stateCount + i * stateCount + j] * vector[category * patternCount * stateCount + pattern * stateCount + j];
                    }
                    result[category * patternCount * stateCount + pattern * stateCount + i] = sum;
                }
            }
        }
    }

    private double[] getTripleVectorReduction(double[] first, double[] second, double[] third, boolean multipliedCategoryWeight) {
        assert(first.length == second.length);
        assert(second.length == third.length);

        double[] patternSums = new double[patternCount];
        double[] categoryProportions = siteRateModel.getCategoryProportions();

        for (int category = 0; category < categoryCount; category++) {
            final double categoryWeight = categoryProportions[category] * (multipliedCategoryWeight ? 1 : siteRateModel.getRateForCategory(category));

            for (int pattern = 0; pattern < patternCount; pattern++) {

                double sum = 0;
                for (int state = 0; state < stateCount; state++) {
                    final int currentIndex = category * stateCount * patternCount + pattern * stateCount + state;
                    sum += first[currentIndex] * second[currentIndex] * third[currentIndex];
                }

                patternSums[pattern] += categoryWeight * sum;

            }
        }
        return patternSums;

    }

    private double getNormalizedPatternReduction(double[] numerator, double[] denominator) {
        double sum = 0;
        for (int pattern = 0; pattern < patternCount; pattern++) {
            sum += numerator[pattern] / denominator[pattern] * patternList.getPatternWeight(pattern);
        }
        return sum;
    }

    private double[] getDoubleVectorReduction(double[] first, double[] second, boolean multipliedCategoryWeight) {
        assert(first.length == second.length);

        double[] patternSums = new double[patternCount];
        double[] categoryProportions = siteRateModel.getCategoryProportions();
        for (int category = 0; category < categoryCount; category++) {
            final double categoryWeight = categoryProportions[category] * (multipliedCategoryWeight ? 1 : siteRateModel.getRateForCategory(category));

            for (int pattern = 0; pattern < patternCount; pattern++) {

                double sum = 0;
                for (int state = 0; state < stateCount; state++) {
                    final int currentIndex = category * stateCount * patternCount + pattern * stateCount + state;
                    sum += first[currentIndex] * second[currentIndex];
                }
                patternSums[pattern] += sum * categoryWeight;

            }

        }
        return patternSums;
    }

    protected String getGradientTraitName() {
        return GRADIENT_TRAIT_NAME;
    }

    protected String getHessianTraitName() {
        return HESSIAN_TRAIT_NAME;
    }

}
