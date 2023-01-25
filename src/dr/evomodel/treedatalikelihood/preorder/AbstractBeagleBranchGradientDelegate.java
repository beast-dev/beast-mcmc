/*
 * DataSimulationDelegate.java
 *
 * Copyright (c) 2002-2020 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.preorder;

import beagle.Beagle;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;

/**
 * AbstractBeagleGradientDelegate - interface for a plugin delegate for data simulation on a tree.
 *
 * @author Xiang Ji
 * @author Marc Suchard
 */
public abstract class AbstractBeagleBranchGradientDelegate extends AbstractBeagleGradientDelegate {

    protected AbstractBeagleBranchGradientDelegate(String name,
                                                   Tree tree,
                                                   BeagleDataLikelihoodDelegate likelihoodDelegate) {
        super(name, tree, likelihoodDelegate);
    }

    abstract protected void cacheDifferentialMassMatrix(Tree tree, boolean cacheSquaredMatrix);

    @Override
    protected int getGradientLength() {
        return tree.getNodeCount() - 1;
    }
    @Override
    protected void getNodeDerivatives(Tree tree, double[] first, double[] second) {

        final int[] postBufferIndices = new int[tree.getNodeCount() - 1];
        final int[] preBufferIndices = new int[tree.getNodeCount() - 1];
        final int[] firstDervIndices = new int[tree.getNodeCount() - 1];
        final int[] secondDeriveIndices = new int[tree.getNodeCount() - 1];

        boolean needsUpdate = !substitutionProcessKnown || second != null;
        if (needsUpdate) {
            cacheDifferentialMassMatrix(tree, second != null);
            substitutionProcessKnown = true;
        }

        int u = 0;
        for (int nodeNum = 0; nodeNum < tree.getNodeCount(); nodeNum++) {
            if (!tree.isRoot(tree.getNode(nodeNum))) {
                postBufferIndices[u] = getPostOrderPartialIndex(nodeNum);
                preBufferIndices[u]  = getPreOrderPartialIndex(nodeNum);
                firstDervIndices[u]  = getFirstDerivativeMatrixBufferIndex(nodeNum);
                secondDeriveIndices[u] = getSecondDerivativeMatrixBufferIndex(nodeNum);
                u++;
            }
        }

        double[] firstSquared = (second != null) ? new double[second.length] : null;

        beagle.calculateEdgeDifferentials(postBufferIndices, preBufferIndices,
                firstDervIndices, new int[] { 0 }, tree.getNodeCount() - 1,
                null, first, firstSquared);

        if (second != null) {
            beagle.calculateEdgeDifferentials(postBufferIndices, preBufferIndices,
                    secondDeriveIndices, new int[] { 0 }, tree.getNodeCount() - 1,
                    null, second, null);

            for (int i = 0; i < second.length; ++i) {
                second[i] -= firstSquared[i];
            }
        }

        if (DEBUG) {
            checkReduction(first);
        }
    }

    private void checkReduction(double[] array) {
        for (double v : array) {
            if (Double.isNaN(v)) {
                double[] postPartial = new double[patternCount * stateCount * categoryCount];
                double[] rootPostPartial = new double[patternCount * stateCount * categoryCount];
                double[] prePartial = new double[patternCount * stateCount * categoryCount];
                double[] differentialMatrix = new double[stateCount * stateCount * categoryCount];
                beagle.getPartials(getPostOrderPartialIndex(0), Beagle.NONE, postPartial);
                beagle.getPartials(getPostOrderPartialIndex(tree.getRoot().getNumber()), Beagle.NONE, rootPostPartial);
                beagle.getPartials(getPreOrderPartialIndex(0), Beagle.NONE, prePartial);
                beagle.getTransitionMatrix(getFirstDerivativeMatrixBufferIndex(0), differentialMatrix);

                double[] grandNumerator = new double[patternCount];
                double[] grandDenominator = new double[patternCount];

                for (int category = 0; category < categoryCount; category++) {
                    final double weight = siteRateModel.getProportionForCategory(category);
                    for (int pattern = 0; pattern < patternCount; pattern++) {
                        double numerator = 0.0;
                        double denominator = 0.0;
                        for (int j = 0; j < stateCount; j++) {
                            double sumOverState = 0.0;
                            for (int k = 0; k < stateCount; k++) {
                                sumOverState += differentialMatrix[stateCount * stateCount * category + stateCount * j + k]
                                        * postPartial[stateCount * patternCount * category + stateCount * pattern + k];

                            }
                            numerator += sumOverState * prePartial[stateCount * patternCount * category + stateCount * pattern + j];
                            denominator += postPartial[stateCount * patternCount * category + stateCount * pattern + j] * prePartial[stateCount * patternCount * category + stateCount * pattern + j];
                        }

                        grandNumerator[pattern] += weight * numerator;
                        grandDenominator[pattern] += weight * denominator;
                    }
                }

//                double sumDeriv = 0.0;
//                for (int j = 0; j < patternCount; j++) {
//                    sumDeriv += grandNumerator[j] / grandDenominator[j] * patternList.getPatternWeight(j);
//                }
                double[] rootFrequencies = evolutionaryProcessDelegate.getRootStateFrequencies();
                double[] patternProb = new double[patternCount];
                for (int category = 0; category < categoryCount; category++) {
                    final double weight = siteRateModel.getProportionForCategory(category);
                    for (int pattern = 0; pattern < patternCount; pattern++) {
                        double sumOverState = 0.0;
                        for (int j = 0; j < stateCount; j++) {
                            sumOverState += rootPostPartial[stateCount * patternCount * category + stateCount * pattern + j] * rootFrequencies[j];
                        }
                        patternProb[pattern] += sumOverState * weight;
                    }
                }
            }
        }
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
                return Intent.BRANCH;
            }

            @Override
            public double[] getTrait(Tree tree, NodeRef node) {
                return getGradient(node);
            }
        });

        treeTraitHelper.addTrait(new TreeTrait.DA() {
            @Override
            public String getTraitName() {
                return getHessianTraitName();
            }

            @Override
            public Intent getIntent() {
                return Intent.BRANCH;
            }

            @Override
            public double[] getTrait(Tree tree, NodeRef node) {
                return getHessian(tree, node);
            }
        });
    }
}
