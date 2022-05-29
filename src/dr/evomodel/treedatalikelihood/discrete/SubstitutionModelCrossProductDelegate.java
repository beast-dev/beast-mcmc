/*
 * DiscreteTraitBranchSubstitutionParameterDelegate.java
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.preorder.AbstractBeagleGradientDelegate;

import java.util.Arrays;

/**
 * @author Marc A. Suchard
 */
public class SubstitutionModelCrossProductDelegate extends AbstractBeagleGradientDelegate {

    private final String name;
    private final Tree tree;
    private final BranchRateModel branchRateModel;
    private final BranchModel branchModel;
    private final int stateCount;
    private final int substitutionModelCount;

    private static final String GRADIENT_TRAIT_NAME = "substitutionModelCrossProductGradient";

    public SubstitutionModelCrossProductDelegate(String name,
                                                 Tree tree,
                                                 BeagleDataLikelihoodDelegate likelihoodDelegate,
                                                 BranchRateModel branchRateModel,
                                                 int stateCount) {

        super(name, tree, likelihoodDelegate);
        this.name = name;
        this.tree = tree;
        this.stateCount = stateCount;
        this.branchRateModel = branchRateModel;
        this.branchModel = likelihoodDelegate.getBranchModel();
        this.substitutionModelCount = branchModel.getSubstitutionModels().size();
    }

    private double getBranchLength(NodeRef node) {

        final double branchRate;
        synchronized (branchRateModel) {
            branchRate = branchRateModel.getBranchRate(tree, node);
        }

        double parentHeight = tree.getNodeHeight(tree.getParent(node));
        double nodeHeight = tree.getNodeHeight(node);

        return branchRate * (parentHeight - nodeHeight);
    }

    @Override
    protected int getGradientLength() {
        return stateCount * stateCount * substitutionModelCount;
    }

    private int coverWholeTree(int[] postBufferIndices,
                                int[] preBufferIndices,
                                double[] branchLengths) {
        int u = 0;
        for (int nodeNum = 0; nodeNum < tree.getNodeCount(); nodeNum++) {
            NodeRef node = tree.getNode(nodeNum);
            if (!tree.isRoot(tree.getNode(nodeNum))) {
                postBufferIndices[u] = getPostOrderPartialIndex(nodeNum);
                preBufferIndices[u]  = getPreOrderPartialIndex(nodeNum);
                branchLengths[u] = getBranchLength(node);
                u++;
            }
        }
        return tree.getNodeCount() - 1;
    }

    private int coverPartialTree(int modelNumber,
                                 int[] postBufferIndices,
                                 int[] preBufferIndices,
                                 double[] branchLengths) {
        int u = 0;
        for (int nodeNum = 0; nodeNum < tree.getNodeCount(); nodeNum++) {
            NodeRef node = tree.getNode(nodeNum);
            if (!tree.isRoot(tree.getNode(nodeNum))) {

                BranchModel.Mapping mapping = branchModel.getBranchModelMapping(node);
                int[] order = mapping.getOrder();
                double[] weights = mapping.getWeights();
                double sum = 0.0;
                for (double w : weights) {
                    sum += w;
                }

                for (int k = 0; k < order.length; ++k) {
                    if (order[k] == modelNumber) {
                        postBufferIndices[u] = getPostOrderPartialIndex(nodeNum);
                        preBufferIndices[u]  = getPreOrderPartialIndex(nodeNum);
                        branchLengths[u] = getBranchLength(node) * weights[k] / sum;
                        u++;
                    }
                }
            }
        }
        return u;
    }

    @Override
    protected void getNodeDerivatives(Tree tree, double[] first, double[] second) {

        final int length = stateCount * stateCount;

        assert (first.length >= length * substitutionModelCount);
        assert (second == null || second.length >= stateCount * stateCount * substitutionModelCount);

        if (second != null) {
            throw new RuntimeException("Not yet implemented");
        }

        final int[] postBufferIndices = new int[tree.getNodeCount() - 1];
        final int[] preBufferIndices = new int[tree.getNodeCount() - 1];
        final double[] branchWeights = new double[tree.getNodeCount() - 1];

        if (substitutionModelCount == 1) {

            Arrays.fill(first, 0, first.length, 0.0);

            int count = coverWholeTree(postBufferIndices, preBufferIndices, branchWeights);
            beagle.calculateCrossProductDifferentials(postBufferIndices, preBufferIndices,
                    new int[] { 0 }, new int[] { 0 },
                    branchWeights,
                    count,
                    first, null);
        } else {

            double[] buffer = new double[length];

            for (int i = 0; i < substitutionModelCount; ++i) {

                Arrays.fill(buffer, 0, buffer.length, 0.0);
                int count = coverPartialTree(i, postBufferIndices, preBufferIndices,
                        branchWeights);
                beagle.calculateCrossProductDifferentials(postBufferIndices, preBufferIndices,
                        new int[] { 0 }, new int[] { 0 },
                        branchWeights,
                        count,
                        buffer, null);

                System.arraycopy(buffer, 0, first, i * length, length);
            }
        }
        // TOOD handle `firstSquared` and `second`
    }

    @Override
    protected String getGradientTraitName() {
        return GRADIENT_TRAIT_NAME + "." + name;
    }

    public static String getName(String name) {
        return GRADIENT_TRAIT_NAME + "." + name;
    }

    @Override
    protected void constructTraits(Helper treeTraitHelper) {

        treeTraitHelper.addTrait(new TreeTrait.DA() {
            @Override
            public String getTraitName() {
                return getGradientTraitName();
            }

            @Override
            public TreeTrait.Intent getIntent() {
                return Intent.WHOLE_TREE;
            }

            @Override
            public double[] getTrait(Tree tree, NodeRef node) {
                return getGradient(node);
            }
        });
    }
}
