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
import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;

/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class DiscreteTraitNodeHeightDelegate extends DiscreteTraitBranchRateDelegate {

    static final String GRADIENT_TRAIT_NAME = "NodeHeightGradient";

    static final String HESSIAN_TRAIT_NAME = "NodeHeightHessian";

    DiscreteTraitNodeHeightDelegate(String name, Tree tree, BeagleDataLikelihoodDelegate likelihoodDelegate) {
        super(name, tree, likelihoodDelegate);
    }


    protected void getNodeDerivatives(Tree tree, double[] first, double[] second) {
        super.getNodeDerivatives(tree, first, second);
        final int internalNodeCount = tree.getInternalNodeCount();

        double[][] prePartials = new double[internalNodeCount][patternCount * stateCount * categoryCount];
        double[][] postPartials = new double[internalNodeCount][patternCount * stateCount * categoryCount];
        double[][] transitionMatrices = new double[internalNodeCount][stateCount * stateCount * categoryCount];

        for (int i = 0; i < internalNodeCount; i++) {
            beagle.getPartials(getPostOrderPartialIndex(i + tree.getExternalNodeCount()), Beagle.NONE, postPartials[i]);
            beagle.getPartials(getPreOrderPartialIndex(i + tree.getExternalNodeCount()), Beagle.NONE, prePartials[i]);
            beagle.getTransitionMatrix(evolutionaryProcessDelegate.getMatrixIndex(i + tree.getExternalNodeCount()), transitionMatrices[i]);
        }

    }

    protected String getGradientTraitName() {
        return GRADIENT_TRAIT_NAME;
    }

    protected String getHessianTraitName() {
        return HESSIAN_TRAIT_NAME;
    }

}
