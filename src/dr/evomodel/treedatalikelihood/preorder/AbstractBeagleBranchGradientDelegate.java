/*
 * AbstractBeagleBranchGradientDelegate.java
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

package dr.evomodel.treedatalikelihood.preorder;

import beagle.Beagle;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.AugmentedNodeRegistry;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;

import java.util.Arrays;

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

        if (epochProcessDelegate != null) {
            getSegmentDerivatives(tree, first, second);
            return;
        }

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

    /**
     * With degree-2 nodes at the epoch transition times, a branch is a chain of segments, each in one epoch. The
     * derivative with respect to the length of a segment uses the infinitesimal matrix of the epoch of the segment.
     * A branch has one rate, so the derivative with respect to its length is the average of the derivatives of its
     * segments, weighted by the fraction of the branch that each segment covers. The chain rule of a caller
     * multiplies by the length of the branch, giving the derivative with respect to the rate of the branch.
     */
    private void getSegmentDerivatives(Tree tree, double[] first, double[] second) {

        if (second != null) {
            throw new UnsupportedOperationException("Second derivatives are not yet supported with augmented " +
                    "epoch nodes, since the segments of a branch share its rate");
        }

        if (!substitutionProcessKnown) {
            cacheDifferentialMassMatrix(tree, false);
            substitutionProcessKnown = true;
        }

        final AugmentedNodeRegistry registry = epochProcessDelegate.getAugmentedNodeRegistry();

        int segmentCount = 0;
        for (int nodeNum = 0; nodeNum < tree.getNodeCount(); nodeNum++) {
            if (!tree.isRoot(tree.getNode(nodeNum))) {
                segmentCount += 1 + registry.getChainLength(nodeNum);
            }
        }

        final int[] postBufferIndices = new int[segmentCount];
        final int[] preBufferIndices = new int[segmentCount];
        final int[] firstDervIndices = new int[segmentCount];
        final int[] branch = new int[segmentCount]; // the position of the branch of each segment in first
        final double[] weight = new double[segmentCount];

        final int[] chain = new int[registry.getBoundaryCount()];

        int s = 0;
        int u = 0;
        for (int nodeNum = 0; nodeNum < tree.getNodeCount(); nodeNum++) {
            final NodeRef node = tree.getNode(nodeNum);
            if (!tree.isRoot(node)) {

                final int count = registry.copyChain(nodeNum, chain);
                final int firstBoundary = registry.getFirstBoundary(nodeNum);

                final double lowestHeight = tree.getNodeHeight(node);
                final double highestHeight = tree.getNodeHeight(tree.getParent(node));

                for (int j = 0; j <= count; ++j) {
                    // the segment above the node (j = 0) and above each of its augmented nodes
                    final int id = (j == 0) ? nodeNum : chain[j - 1];
                    final double lower = (j == 0) ? lowestHeight : registry.getBoundary(firstBoundary + j - 1);
                    final double upper = (j == count) ? highestHeight : registry.getBoundary(firstBoundary + j);

                    postBufferIndices[s] = getPostOrderPartialIndex(id);
                    preBufferIndices[s] = getPreOrderPartialIndex(id);
                    firstDervIndices[s] = getFirstDerivativeMatrixBufferIndex(id);
                    branch[s] = u;
                    weight[s] = (count == 0) ? 1.0 : (upper - lower) / (highestHeight - lowestHeight);
                    ++s;
                }
                u++;
            }
        }

        final double[] segmentDerivatives = new double[segmentCount];

        beagle.calculateEdgeDifferentials(postBufferIndices, preBufferIndices,
                firstDervIndices, new int[] { 0 }, segmentCount,
                null, segmentDerivatives, null);

        Arrays.fill(first, 0.0);
        for (s = 0; s < segmentCount; ++s) {
            first[branch[s]] += weight[s] * segmentDerivatives[s];
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
