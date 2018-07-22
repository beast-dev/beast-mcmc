/*
 * BranchRateGradientForDiscreteTrait.java
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

import dr.evolution.tree.*;
import dr.evomodel.treedatalikelihood.*;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.loggers.Loggable;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

import static dr.math.MachineAccuracy.SQRT_EPSILON;

/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class BranchRateGradientForDiscreteTrait extends DiscreteTraitBranchRateGradient
        implements GradientWrtParameterProvider, HessianWrtParameterProvider, Reportable, Loggable {

    // TODO Refactor / remove code duplication with BranchRateGradient
    // TODO Maybe use:  AbstractBranchRateGradient, DiscreteTraitBranchRateGradient, ContinuousTraitBranchRateGradient
    public BranchRateGradientForDiscreteTrait(String traitName,
                                              TreeDataLikelihood treeDataLikelihood,
                                              BeagleDataLikelihoodDelegate likelihoodDelegate,
                                              Parameter rateParameter,
                                              boolean useHessian) {
        super(traitName, treeDataLikelihood, likelihoodDelegate, rateParameter, useHessian);
    }

    @Override
    public double[] getDiagonalHessianLogDensity() {

        double[] result = new double[rateParameter.getDimension()];

        //Do single call to traitProvider with node == null (get full tree)
        double[] gradient =  (double[]) treeTraitProvider.getTrait(tree, null);
        double[] diagonalHessian =  (double[]) treeDataLikelihood.getTreeTrait("Hessian").getTrait(tree, null);

        int v =0;
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                final int destinationIndex = getParameterIndexFromNode(node);
                final double differential = branchRateModel.getBranchRateDifferential(tree, node);
                final double secondDifferential = branchRateModel.getBranchRateSecondDifferential(tree, node);
                final double branchLength = tree.getBranchLength(node);
                result[destinationIndex] =  branchLength *
                        (gradient[v] * secondDifferential + diagonalHessian[v] * differential * differential * branchLength);
                v++;
            }
        }

        return result;
    }

    @Override
    public double[] getGradientLogDensity() {

        double[] result = new double[rateParameter.getDimension()];

        //Do single call to traitProvider with node == null (get full tree)
        double[] gradient =  (double[]) treeTraitProvider.getTrait(tree, null);

        int v =0;
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                final int destinationIndex = getParameterIndexFromNode(node);
                final double differential = branchRateModel.getBranchRateDifferential(tree, node);
                final double nodeResult = gradient[v] * differential * tree.getBranchLength(node);
                if (Double.isNaN(nodeResult) && !Double.isInfinite(treeDataLikelihood.getLogLikelihood())) {
                    System.err.println("Check Gradient calculation please.");
                }
                result[destinationIndex] = nodeResult;
                v++;
            }
        }

        if (COUNT_TOTAL_OPERATIONS) {
            ++getGradientLogDensityCount;
        }

        return result;
    }

}
