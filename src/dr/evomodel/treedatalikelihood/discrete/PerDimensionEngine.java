/*
 * PerDimensionEngine.java
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

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.substmodel.DifferentiableSubstitutionModel;
import dr.evomodel.substmodel.DifferentialMassProvider;
import dr.evomodel.substmodel.DifferentialMassProvider.DifferentialWrapper.WrtParameter;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Branch term one dimension at a time through the exact Psi-form differential mass matrices (one pre-order pass
 * per dimension). Fallback and parity reference for {@link AdjointEngine}.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
final class PerDimensionEngine implements BranchTermEngine {

    private final Tree tree;
    private final List<TreeTrait> traits = new ArrayList<>();

    PerDimensionEngine(TreeDataLikelihood treeDataLikelihood,
                       BeagleDataLikelihoodDelegate likelihoodDelegate,
                       DifferentiableSubstitutionModel substitutionModel,
                       List<WrtParameter> wrts,
                       String namePrefix) {

        String reason = whyIneligible(likelihoodDelegate);
        if (reason != null) {
            throw new IllegalArgumentException(reason);
        }

        this.tree = treeDataLikelihood.getTree();

        for (int d = 0; d < wrts.size(); ++d) {
            final String traitName = namePrefix + ":" + d;
            final String name = BranchSubstitutionParameterDelegate.getName(traitName);

            if (treeDataLikelihood.getTreeTrait(name) == null) {
                DifferentialMassProvider provider = new DifferentialMassProvider.DifferentialWrapper(
                        substitutionModel, wrts.get(d), DifferentialMassProvider.Mode.EXACT);
                BranchDifferentialMassProvider branchProvider = new BranchDifferentialMassProvider(null,
                        Collections.singletonList(provider));
                ProcessSimulationDelegate gradientDelegate = new BranchSubstitutionParameterDelegate(traitName,
                        tree, likelihoodDelegate, treeDataLikelihood.getBranchRateModel(), branchProvider);
                TreeTraitProvider traitProvider = new ProcessSimulation(treeDataLikelihood, gradientDelegate);
                treeDataLikelihood.addTraits(traitProvider.getTreeTraits());
            }
            traits.add(treeDataLikelihood.getTreeTrait(name));
        }
    }

    static String whyIneligible(BeagleDataLikelihoodDelegate likelihoodDelegate) {
        if (!likelihoodDelegate.isUsePreOrder()
                || !likelihoodDelegate.getPreOrderSettings().isBranchInfinitesimalDerivative()) {
            return "the per-dimension path needs usePreOrder=\"true\" and branchInfinitesimalDerivative=\"true\"";
        }
        if (likelihoodDelegate.getBranchModel().getSubstitutionModels().size() != 1) {
            return "the per-dimension path needs a single substitution model";
        }
        return null;
    }

    @Override
    public double[] getBranchGradient() {
        final double[] gradient = new double[traits.size()];
        for (int d = 0; d < traits.size(); ++d) {
            double[] branches = (double[]) traits.get(d).getTrait(tree, null);
            for (double value : branches) {
                gradient[d] += value; // sum over all branches
            }
        }
        return gradient;
    }

    @Override
    public String getName() {
        return "perDimension";
    }
}
