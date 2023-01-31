/*
 * BranchSubstitutionParameterLocationGradient.java
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
import dr.evomodel.branchratemodel.BranchSpecificFixedEffects;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.BranchParameter;


/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class BranchSubstitutionParameterLocationGradient extends HyperParameterBranchSubstitutionParameterGradient {

    private final BranchSpecificFixedEffects fixedEffects;

    public BranchSubstitutionParameterLocationGradient(String traitName,
                                                       TreeDataLikelihood treeDataLikelihood,
                                                       BeagleDataLikelihoodDelegate likelihoodDelegate,
//                                                       Parameter branchSubstitutionParameter,
                                                       BranchParameter branchParameter,
                                                       Double tolerance,
                                                       boolean useHessian,
                                                       BranchSpecificFixedEffects fixedEffects) {
        super(traitName, treeDataLikelihood, likelihoodDelegate, branchParameter,
                fixedEffects.getFixedEffectsParameter(), tolerance, useHessian);

        this.fixedEffects = fixedEffects;
    }

    @Override
    double[] getDifferential(Tree tree, NodeRef node) {
        double rate = branchRateModel.getBranchRate(tree, node);
        double[] results = fixedEffects.getDifferential(rate, tree, node);

        return results;
    }
}
