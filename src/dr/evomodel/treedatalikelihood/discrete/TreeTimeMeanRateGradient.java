/*
 * TreeTimeMeanRateGradient.java
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.ScaledByTreeTimeBranchRateModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;

/**
 * Gradient of a tree likelihood with respect to the mean rate of a branch-rate model scaled by tree time (the
 * location of the shrinkage local clock). The mean rate m multiplies every branch rate, so
 * d log L / d m = sum_e (d log L / d rate_e) rate_e / m.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class TreeTimeMeanRateGradient extends DiscreteTraitBranchRateGradient {

    public TreeTimeMeanRateGradient(String traitName,
                                    TreeDataLikelihood treeDataLikelihood,
                                    BeagleDataLikelihoodDelegate likelihoodDelegate) {
        super(traitName, treeDataLikelihood, likelihoodDelegate, getMeanRate(treeDataLikelihood), false);
    }

    private static dr.inference.model.Parameter getMeanRate(TreeDataLikelihood treeDataLikelihood) {
        if (!(treeDataLikelihood.getBranchRateModel() instanceof ScaledByTreeTimeBranchRateModel)) {
            throw new IllegalArgumentException("The branch rate model of " + treeDataLikelihood.getId()
                    + " is not scaled by tree time");
        }
        dr.inference.model.Parameter meanRate =
                ((ScaledByTreeTimeBranchRateModel) treeDataLikelihood.getBranchRateModel()).getMeanRateParameter();
        if (meanRate == null) {
            throw new IllegalArgumentException("The branch rate model of " + treeDataLikelihood.getId()
                    + " has no mean rate parameter");
        }
        return meanRate;
    }

    // result[v] is d log L / d rate_v, indexed by branch (node number with the root skipped)
    @Override
    double[] updateBranchRateGradientLogDensity(double[] result) {
        final double meanRate = rateParameter.getParameterValue(0);
        double gradient = 0.0;
        int v = 0;
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                gradient += result[v] * treeDataLikelihood.getBranchRateModel().getBranchRate(tree, node) / meanRate;
                v++;
            }
        }
        return new double[]{gradient};
    }
}
