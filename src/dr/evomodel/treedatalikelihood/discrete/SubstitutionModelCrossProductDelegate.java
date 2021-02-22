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
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.substmodel.GLMSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.preorder.AbstractBeagleGradientDelegate;

import java.util.Arrays;

import static dr.evomodel.treedatalikelihood.discrete.DiscreteTraitBranchRateDelegate.scaleInfinitesimalMatrixByRates;

/**
 * @author Marc A. Suchard
 */
public class SubstitutionModelCrossProductDelegate extends AbstractBeagleGradientDelegate {

    private final String name;
    private final Tree tree;
    private final BranchRateModel branchRateModel;
    private final int stateCount;

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
        return stateCount * stateCount;
    }

    @Override
    protected void getNodeDerivatives(Tree tree, double[] first, double[] second) {

        assert (first.length >= stateCount * stateCount);
        assert (second == null || second.length >= stateCount * stateCount);

        final int[] postBufferIndices = new int[tree.getNodeCount() - 1];
        final int[] preBufferIndices = new int[tree.getNodeCount() - 1];
        final double[] branchLengths = new double[tree.getNodeCount() - 1];

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

        Arrays.fill(first, 0, first.length, 0.0);
        double[] firstSquared = (second != null) ? new double[second.length] : null;

        beagle.calculateCrossProductDifferentials(postBufferIndices, preBufferIndices,
                new int[] { 0 }, new int[] { 0 },
                branchLengths,
                tree.getNodeCount() - 1,
                first, firstSquared);

        if (second != null) {
//            beagle.calculateEdgeDifferentials(postBufferIndices, preBufferIndices,
//                    secondDeriveIndices, new int[] { 0 }, tree.getNodeCount() - 1,
//                    null, second, null);
//
//            for (int i = 0; i < second.length; ++i) {
//                second[i] -= firstSquared[i];
//            }
            throw new RuntimeException("Not yet implemented");
        }
    }

//    private void cacheDifferentialMassMatrix(Tree tree, boolean b) {
//    }


    protected String getGradientTraitName() {
        return GRADIENT_TRAIT_NAME + "." + name;
    }

//    protected String getHessianTraitName() {
//        return HESSIAN_TRAIT_NAME + ":" + name;
//    }

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

//    private static final boolean DEBUG = false;
}
