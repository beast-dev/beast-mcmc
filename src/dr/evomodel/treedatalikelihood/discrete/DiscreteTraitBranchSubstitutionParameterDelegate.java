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
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.preorder.AbstractDiscreteTraitDelegate;

/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class DiscreteTraitBranchSubstitutionParameterDelegate extends AbstractDiscreteTraitDelegate {

    private final BranchRateModel branchRateModel;
    private final BranchDifferentialMassProvider branchDifferentialMassProvider;
    private final String name;

    private static final String GRADIENT_TRAIT_NAME = "BranchSubstitutionGradient";
    private static final String HESSIAN_TRAIT_NAME = "BranchSubstitutionHessian";

    DiscreteTraitBranchSubstitutionParameterDelegate(String name,
                                                     Tree tree,
                                                     BeagleDataLikelihoodDelegate likelihoodDelegate,
                                                     BranchRateModel branchRateModel,
                                                     BranchDifferentialMassProvider branchDifferentialMassProvider) {
        super(name, tree, likelihoodDelegate);
        this.name = name;
        this.branchRateModel = branchRateModel;
        this.branchDifferentialMassProvider = branchDifferentialMassProvider;
    }

    @Override
    protected void cacheDifferentialMassMatrix(Tree tree, boolean cacheSquaredMatrix) {
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {

                final double time = tree.getBranchLength(node) * branchRateModel.getBranchRate(tree, node);
                double[] differentialMassMatrix = branchDifferentialMassProvider.getDifferentialMassMatrixForBranch(node, time);
                evolutionaryProcessDelegate.cacheFirstOrderDifferentialMatrix(beagle, i, differentialMassMatrix);
            }
        }
        if (cacheSquaredMatrix) {
            throw new RuntimeException("Not yet implemented!");
        }
    }

    @Override
    protected int getFirstDerivativeMatrixBufferIndex(int nodeNum) {
        return evolutionaryProcessDelegate.getFirstOrderDifferentialMatrixBufferIndex(nodeNum);
    }

    @Override
    protected int getSecondDerivativeMatrixBufferIndex(int nodeNum) {
        return evolutionaryProcessDelegate.getSecondOrderDifferentialMatrixBufferIndex(nodeNum);
    }

    protected String getGradientTraitName() {
        return GRADIENT_TRAIT_NAME + ":" + name;
    }

    protected String getHessianTraitName() {
        return HESSIAN_TRAIT_NAME + ":" + name;
    }

    public static String getName(String name) {
        return GRADIENT_TRAIT_NAME + ":" + name;
    }

}
