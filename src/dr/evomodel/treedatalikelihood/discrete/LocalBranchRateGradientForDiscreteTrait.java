/*
 * LocalBranchRateGradientForDiscreteTrait.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.loggers.Loggable;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

/**
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class LocalBranchRateGradientForDiscreteTrait extends DiscreteTraitBranchRateGradient
        implements GradientWrtParameterProvider, Reportable, Loggable {

    public LocalBranchRateGradientForDiscreteTrait(String traitName,
                                                   TreeDataLikelihood treeDataLikelihood,
                                                   BeagleDataLikelihoodDelegate likelihoodDelegate,
                                                   Parameter rateParameter,
                                                   boolean useHessian) {
        super(traitName, treeDataLikelihood, likelihoodDelegate, rateParameter, useHessian);
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] nodeGradient = super.getGradientLogDensity();
        double[] result = new double[tree.getNodeCount() - 1];

        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if(!tree.isRoot(node)){
                final int parameterIndex = branchRateModel.getParameterIndexFromNode(node);
                result[parameterIndex] = getSubTreeGradient(tree, node, nodeGradient);
            }
        }
        return result;
    }

    private double getSubTreeGradient(Tree tree, NodeRef node, final double[] nodeGradient) {
        final double chainGradient =  branchRateModel.getBranchRateDifferential(tree, node);
        final double branchGradient = nodeGradient[branchRateModel.getParameterIndexFromNode(node)];
        if (tree.isExternal(node)) {
            return chainGradient * branchGradient;
        } else {
            double sum = chainGradient * branchGradient;
            for (int i = 0; i < tree.getChildCount(node); i++) {
                sum += getSubTreeGradient(tree, tree.getChild(node, i), nodeGradient);
            }
            return sum;
        }
    }
}
