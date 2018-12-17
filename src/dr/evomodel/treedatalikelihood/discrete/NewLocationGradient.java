/*
 * NewLocationGradient.java
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
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

//TODO: use it to replace code duplication
public class NewLocationGradient extends HyperParameterGradient {

    private final BranchSpecificFixedEffects fixedEffects;

    public NewLocationGradient(TreeDataLikelihood treeDataLikelihood,
                               GradientWrtParameterProvider gradientWrtParameterProvider,
                               BranchSpecificFixedEffects fixedEffects,
                               TreeParameterModel branchParameter,
                               boolean useHessian) {

        super(treeDataLikelihood, gradientWrtParameterProvider, fixedEffects.getFixedEffectsParameter(), branchParameter, useHessian);

        this.fixedEffects = fixedEffects;

    }

    @Override
    double[] getDifferential(Tree tree, NodeRef node) {
        double rate = branchParameter.getNodeValue(tree, node);
        double[] results = fixedEffects.getDifferential(rate, tree, node);

        return results;
    }

    @Override
    public double[][] getHessianLogDensity() {
        throw new RuntimeException("Not yet implemented!");
    }
}
