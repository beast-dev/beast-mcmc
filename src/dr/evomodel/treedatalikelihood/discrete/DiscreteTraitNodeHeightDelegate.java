/*
 * DiscreteTraitNodeHeightDelegate.java
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

        if (first != null) {
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
        }



        if (second != null) {
            final int internalNodeCount = tree.getInternalNodeCount();
            double[][] sisterBranchesSecondDerivatives = new double[internalNodeCount][];
            double[][] currentAndParentBranchesSecondDerivatives = new double[tree.getNodeCount() - 1][];

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

            double[][] branchPatternGradient = new double[tree.getNodeCount()][];
            double[][] branchPatternDiagonalHessian = new double[tree.getNodeCount()][];

            double[] tmpPatternNumerator = new double[patternCount];
            double[] tmpPatternDenominator = new double[patternCount];

            double[] tmpPatternCacheOne = new double[patternCount];
            double[] tmpPatternCacheTwo = new double[patternCount];

            for (int i = 0; i < tree.getNodeCount() - 1; i++) {
                beagle.getTransitionMatrix(evolutionaryProcessDelegate.getInfinitesimalMatrixBufferIndex(i), Qi);
                getMatrixVectorProduct(Qi, postPartials[i], tmpIPartial);
                getMatrixVectorProduct(Qi, tmpIPartial, tmpLeftPartail);
                denominator[i] = getVectorStateReduction(getVectorVectorProduct(postPartials[i], prePartials[i]));
                branchPatternGradient[i] = getVectorVectorDivision(getVectorStateReduction(getVectorVectorProduct(prePartials[i], tmpIPartial)), denominator[i]);
                branchPatternDiagonalHessian[i] = getVectorMinusVector(getVectorVectorDivision(getVectorStateReduction(getVectorVectorProduct(prePartials[i], tmpLeftPartail)), denominator[i]),
                        getVectorVectorProduct(branchPatternGradient[i], branchPatternGradient[i]));
            }
            denominator[tree.getRoot().getNumber()] = getDoubleVectorReduction(postPartials[tree.getRoot().getNumber()], prePartials[tree.getRoot().getNumber()], true);

            for (int i = 0; i < internalNodeCount; i++) {
                NodeRef nodeI = tree.getNode(i + tree.getExternalNodeCount());
                NodeRef nodeJ = tree.getChild(nodeI, 0);
                NodeRef nodeK = tree.getChild(nodeI, 1);

                beagle.getTransitionMatrix(evolutionaryProcessDelegate.getInfinitesimalMatrixBufferIndex(nodeI.getNumber()), Qi);
                beagle.getTransitionMatrix(evolutionaryProcessDelegate.getInfinitesimalMatrixBufferIndex(nodeJ.getNumber()), Qj);
                beagle.getTransitionMatrix(evolutionaryProcessDelegate.getInfinitesimalMatrixBufferIndex(nodeK.getNumber()), Qk);

                getMatrixVectorProduct(transitionMatrices[nodeJ.getNumber()], postPartials[nodeJ.getNumber()], tmpLeftPartail);  // Pj p_j
                getMatrixVectorProduct(Qj, tmpLeftPartail, tmpQLeftPartial); // QjPj p_j

                getMatrixVectorProduct(transitionMatrices[nodeK.getNumber()], postPartials[nodeK.getNumber()], tmpRightPartial); // Pk p_k
                getMatrixVectorProduct(Qk, tmpRightPartial, tmpQRightPartial); // Qk Pk p_k

                getTripleVectorReduction(tmpQLeftPartial, tmpQRightPartial, prePartials[nodeI.getNumber()], tmpPatternNumerator);
                getTripleVectorReduction(tmpLeftPartail, tmpRightPartial, prePartials[nodeI.getNumber()], tmpPatternDenominator);
                getVectorVectorDivision(tmpPatternNumerator, tmpPatternDenominator, tmpPatternCacheOne);
                getVectorVectorProduct(branchPatternGradient[getParameterIndex(nodeJ, tree)], branchPatternGradient[getParameterIndex(nodeK, tree)], tmpPatternCacheTwo);

                sisterBranchesSecondDerivatives[i] = getVectorMinusVector(tmpPatternCacheOne, tmpPatternCacheTwo);

                getVectorPlusScaledVector(branchPatternDiagonalHessian[getParameterIndex(nodeJ, tree)], branchPatternDiagonalHessian[getParameterIndex(nodeK, tree)],
                        branchRates.getBranchRate(tree, nodeJ) * branchRates.getBranchRate(tree, nodeJ),
                        branchRates.getBranchRate(tree, nodeK) * branchRates.getBranchRate(tree, nodeK), tmpPatternCacheOne);
                getVectorPlusScaledVector(tmpPatternCacheOne,
                        sisterBranchesSecondDerivatives[i], 1, 2 * branchRates.getBranchRate(tree, nodeJ) * branchRates.getBranchRate(tree, nodeK),
                        tmpPatternCacheTwo);

                second[i] = getVectorPatternReduction(tmpPatternCacheTwo);


                if (!tree.isRoot(nodeI)) {
                    getMatrixTransformVectorProduct(Qi, prePartials[nodeI.getNumber()], tmpIPartial);  //Qi q_i

                    getTripleVectorReduction(tmpQLeftPartial, tmpRightPartial, tmpIPartial, tmpPatternNumerator);
                    getTripleVectorReduction(tmpLeftPartail, tmpRightPartial, prePartials[nodeI.getNumber()], tmpPatternDenominator);
                    getVectorVectorDivision(tmpPatternNumerator, tmpPatternDenominator, tmpPatternCacheOne);
                    getVectorVectorProduct(branchPatternGradient[getParameterIndex(nodeJ, tree)], branchPatternGradient[getParameterIndex(nodeI, tree)], tmpPatternCacheTwo);

                    currentAndParentBranchesSecondDerivatives[getParameterIndex(nodeJ, tree)] = getVectorMinusVector(tmpPatternCacheOne, tmpPatternCacheTwo);

                    getTripleVectorReduction(tmpQRightPartial, tmpLeftPartail, tmpIPartial, tmpPatternNumerator);
                    getTripleVectorReduction(tmpLeftPartail, tmpRightPartial, prePartials[nodeI.getNumber()], tmpPatternDenominator);
                    getVectorVectorDivision(tmpPatternNumerator, tmpPatternDenominator, tmpPatternCacheOne);
                    getVectorVectorProduct(branchPatternGradient[getParameterIndex(nodeK, tree)], branchPatternGradient[getParameterIndex(nodeI, tree)], tmpPatternCacheTwo);

                    currentAndParentBranchesSecondDerivatives[getParameterIndex(nodeK, tree)] = getVectorMinusVector(tmpPatternCacheOne, tmpPatternCacheTwo);

                    getVectorPlusScaledVector(branchPatternDiagonalHessian[getParameterIndex(nodeI, tree)], currentAndParentBranchesSecondDerivatives[getParameterIndex(nodeJ, tree)],
                            branchRates.getBranchRate(tree, nodeI) * branchRates.getBranchRate(tree, nodeI),
                            -2 * branchRates.getBranchRate(tree, nodeI) * branchRates.getBranchRate(tree, nodeJ), tmpPatternCacheOne);
                    getVectorPlusScaledVector(tmpPatternCacheOne, currentAndParentBranchesSecondDerivatives[getParameterIndex(nodeK, tree)],
                            1, -2 * branchRates.getBranchRate(tree, nodeI) * branchRates.getBranchRate(tree, nodeK),
                            tmpPatternCacheTwo);

                    second[i] +=  getVectorPatternReduction(tmpPatternCacheTwo);

                }

            }

        }
    }

    private double[] getVectorMinusVector(double[] first, double[] second) {
        double[] out = new double[first.length];
        for (int i = 0; i < first.length; i++) {
            out[i] = first[i] - second[i];
        }
        return out;
    }


    private void getVectorPlusScaledVector(double[] first, double[] second, double firstScale, double secondScale, double[] out) {

        if (second != null) {
            for (int i = 0; i < first.length; i++) {
                out[i] = first[i] * firstScale + second[i] * secondScale;
            }
        } else {
            for (int i = 0; i < first.length; i++) {
                out[i] = first[i] * firstScale;
            }
        }

    }


    private int getParameterIndex(NodeRef node, Tree tree) {
        return node.getNumber() < tree.getRoot().getNumber() ? node.getNumber() : node.getNumber() - 1;
    }

    private double[] getVectorStateReduction(double[] vector) {
        double[] out = new double[patternCount];
        for (int category = 0; category < categoryCount; category++) {
            final double categoryWeight = siteRateModel.getProportionForCategory(category);
            for (int pattern = 0; pattern < out.length; pattern++) {
                double sum = 0;
                for (int state = 0; state < stateCount; state++) {
                    sum += vector[category * patternCount * stateCount + pattern * stateCount + state];
                }
                out[pattern] += sum * categoryWeight;
            }
        }
        return out;
    }

    private double[] getVectorVectorProduct(double[] first, double[] second) {
        double[] product = new double[first.length];
        for (int i = 0; i < first.length; i++) {
            product[i] = first[i] * second[i];
        }
        return product;
    }

    private void getVectorVectorProduct(double[] first, double[] second, double[] out) {
        for (int i = 0; i < first.length; i++) {
            out[i] = first[i] * second[i];
        }
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

    private void getMatrixTransformVectorProduct(double[] matrix, double[] vector, double[] result) {
        assert(vector.length == result.length);
        for (int category = 0; category < categoryCount; category++) {
            for (int pattern = 0; pattern < patternCount; pattern++) {
                for (int i = 0; i < stateCount; i++) {
                    double sum = 0;
                    for (int j = 0; j < stateCount; j++) {
                        sum += matrix[category * stateCount * stateCount + j * stateCount + i] * vector[category * patternCount * stateCount + pattern * stateCount + j];
                    }
                    result[category * patternCount * stateCount + pattern * stateCount + i] = sum;
                }
            }
        }
    }

    private double getVectorPatternReduction(double[] vector) {
        double sum = 0;

        for (int pattern = 0; pattern < patternCount; pattern++) {
            sum += vector[pattern] * patternList.getPatternWeight(pattern);
        }
        return sum;
    }


    private void getTripleVectorReduction(double[] first, double[] second, double[] third, double[] out) {
        assert(first.length == second.length);
        assert(second.length == third.length);

        Arrays.fill(out, 0);

        for (int category = 0; category < categoryCount; category++) {
            final double categoryWeight = siteRateModel.getProportionForCategory(category);

            for (int pattern = 0; pattern < patternCount; pattern++) {

                double sum = 0;
                for (int state = 0; state < stateCount; state++) {
                    final int currentIndex = category * stateCount * patternCount + pattern * stateCount + state;
                    sum += first[currentIndex] * second[currentIndex] * third[currentIndex];
                }

                out[pattern] += categoryWeight * sum;

            }
        }
    }


    private double[] getVectorVectorDivision(double[] numerator, double[] denominator) {
        double[] out = new double[numerator.length];
        for (int i = 0; i < numerator.length; i++) {
            out[i] = numerator[i] / denominator[i];
        }
        return out;
    }

    private void getVectorVectorDivision(double[] numerator, double[] denominator, double[] out) {
        for (int i = 0; i < numerator.length; i++) {
            out[i] = numerator[i] / denominator[i];
        }
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