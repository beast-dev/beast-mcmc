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


    protected void getNodeDerivatives(Tree tree, double[] first, double[] second) {
        double[] branchGradient = new double[first.length];
        double[] branchDiagonalHessian = second == null ? null : new double[first.length];
        super.getNodeDerivatives(tree, branchGradient, branchDiagonalHessian);

        for (int i = 0; i < tree.getInternalNodeCount(); ++i) {

            final  NodeRef node = tree.getNode(i + tree.getExternalNodeCount());

            for (int j = 0; j < tree.getChildCount(node); j++) {
                NodeRef childNode = tree.getChild(node, j);
                final int childNodeIndex = childNode.getNumber();
                first[i] += branchGradient[childNodeIndex] * branchRates.getBranchRate(tree, childNode);
            }
            if (!tree.isRoot(node)) {
                first[i] -= branchGradient[i + tree.getExternalNodeCount()] * branchRates.getBranchRate(tree, node);
            }
        }

        if (second != null) {
            final int internalNodeCount = tree.getInternalNodeCount();
            double[] sisterBranchesSecondDerivatives = new double[internalNodeCount];
            double[] currentAndParentBranchesSecondDerivatives = new double[tree.getNodeCount() - 2];

            double[][] prePartials = new double[internalNodeCount][patternCount * stateCount * categoryCount];
            double[][] postPartials = new double[internalNodeCount][patternCount * stateCount * categoryCount];
            double[][] transitionMatrices = new double[internalNodeCount][stateCount * stateCount * categoryCount];

            for (int i = 0; i < internalNodeCount; i++) {
                beagle.getPartials(getPostOrderPartialIndex(i + tree.getExternalNodeCount()), Beagle.NONE, postPartials[i]);
                beagle.getPartials(getPreOrderPartialIndex(i + tree.getExternalNodeCount()), Beagle.NONE, prePartials[i]);
                beagle.getTransitionMatrix(evolutionaryProcessDelegate.getMatrixIndex(i + tree.getExternalNodeCount()), transitionMatrices[i]);
            }

            double[] Qi = new double[stateCount * stateCount * categoryCount];
            double[] Qj = new double[stateCount * stateCount * categoryCount];
            double[] Qk = new double[stateCount * stateCount * categoryCount];
            double[] tmpLeftPartail = new double[patternCount * stateCount * categoryCount];
            double[] tmpRightPartial = new double[patternCount * stateCount * categoryCount];
            double[] tmpQLeftPartial = new double[patternCount * stateCount * categoryCount];
            double[] tmpQRightPartial = new double[patternCount * stateCount * categoryCount];
            double[] tmpIPartial = new double[patternCount * stateCount * categoryCount];

            for (int i = 0; i < internalNodeCount; i++) {
                NodeRef nodeI = tree.getNode(i + tree.getExternalNodeCount());
                NodeRef nodeJ = tree.getChild(nodeI, 0);
                NodeRef nodeK = tree.getChild(nodeI, 1);

                beagle.getTransitionMatrix(evolutionaryProcessDelegate.getInfinitesimalMatrixBufferIndex(i + tree.getExternalNodeCount()), Qi);
                beagle.getTransitionMatrix(evolutionaryProcessDelegate.getInfinitesimalMatrixBufferIndex(nodeJ.getNumber()), Qj);
                beagle.getTransitionMatrix(evolutionaryProcessDelegate.getInfinitesimalMatrixBufferIndex(nodeK.getNumber()), Qk);

                getMatrixVectorProduct(transitionMatrices[nodeJ.getNumber()], postPartials[nodeJ.getNumber()], tmpLeftPartail);
                getMatrixVectorProduct(Qj, tmpLeftPartail, tmpQLeftPartial);

                getMatrixVectorProduct(transitionMatrices[nodeK.getNumber()], postPartials[nodeK.getNumber()], tmpRightPartial);
                getMatrixVectorProduct(Qk, tmpRightPartial, tmpQRightPartial);

                getMatrixVectorProduct(Qi, prePartials[i + tree.getExternalNodeCount()], tmpIPartial);

                sisterBranchesSecondDerivatives[i] = getTripleVectorReduction(tmpQLeftPartial, tmpQRightPartial, prePartials[i + tree.getExternalNodeCount()]);
                currentAndParentBranchesSecondDerivatives[nodeJ.getNumber() - tree.getExternalNodeCount()] = getTripleVectorReduction(tmpQLeftPartial, tmpRightPartial, tmpIPartial);
                currentAndParentBranchesSecondDerivatives[nodeK.getNumber() - tree.getExternalNodeCount()] = getTripleVectorReduction(tmpQRightPartial, tmpLeftPartail, tmpIPartial);
            }

            for (int i = 0; i < internalNodeCount; i++) {
                NodeRef nodeI = tree.getNode(i + tree.getExternalNodeCount());
                NodeRef nodeJ = tree.getChild(nodeI, 0);
                NodeRef nodeK = tree.getChild(nodeI, 1);


                if (tree.isRoot(nodeI)) {

                }
            }
        }

        int test = 0;


    }

    private void getMatrixVectorProduct(double[] matrix, double[] vector, double[] result) {

    }

    private double getTripleVectorReduction(double[] first, double[] second, double[] third) {
        assert(first.length == second.length);
        assert(second.length == third.length);
        double tmp = 0;
        for (int i = 0; i < first.length; i++) {
            tmp += first[i] * second[i] * third[i];
        }
        return tmp;

    }

    protected String getGradientTraitName() {
        return GRADIENT_TRAIT_NAME;
    }

    protected String getHessianTraitName() {
        return HESSIAN_TRAIT_NAME;
    }

}
