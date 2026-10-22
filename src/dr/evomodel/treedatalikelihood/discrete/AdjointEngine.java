/*
 * AdjointEngine.java
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
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.substmodel.DifferentiableSubstitutionModel;
import dr.evomodel.substmodel.DifferentialMassProvider.DifferentialWrapper.WrtParameter;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.math.matrixAlgebra.WrappedMatrix;

import java.util.List;

/**
 * Branch term from one adjoint pass: M = d log L / d Q~ (summed over branches and rate categories), then
 * g_d = sum_ij M_ij (dQ~ / d theta_d)_ij for each dimension.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
final class AdjointEngine implements BranchTermEngine {

    // Reserved: the trait factory keys by name only, so this must not collide with a user's first-order gradient.
    // Every provider on one likelihood shares this single adjoint pass.
    static final String ADJOINT_TRAIT = "substitutionModelAdjoint";

    private final SubstitutionModelGradientTraitFactory.Context context;
    private final Tree tree;
    private final DifferentiableSubstitutionModel substitutionModel;
    private final List<WrtParameter> wrts;
    private final int stateCount;

    AdjointEngine(TreeDataLikelihood treeDataLikelihood,
                  BeagleDataLikelihoodDelegate likelihoodDelegate,
                  DifferentiableSubstitutionModel substitutionModel,
                  List<WrtParameter> wrts,
                  int stateCount) {
        this.tree = treeDataLikelihood.getTree();
        this.substitutionModel = substitutionModel;
        this.wrts = wrts;
        this.stateCount = stateCount;
        this.context = SubstitutionModelGradientTraitFactory.create(ADJOINT_TRAIT, treeDataLikelihood,
                likelihoodDelegate, stateCount, true, false, "an exact adjoint");
    }

    /**
     * @return null if the adjoint applies to this likelihood, otherwise the reason it does not
     */
    static String whyIneligible(BeagleDataLikelihoodDelegate likelihoodDelegate) {
        if (!likelihoodDelegate.isUsePreOrder()) {
            return "the adjoint needs usePreOrder=\"true\"";
        }
        if (likelihoodDelegate.getRescalingScheme() == PartialsRescalingScheme.AUTO) {
            return "the adjoint does not support BEAGLE auto-scaling";
        }
        if (likelihoodDelegate.getPreOrderSettings().useHierarchicalRepresentation()) {
            return "the adjoint does not support the hierarchical representation";
        }
        BranchModel branchModel = likelihoodDelegate.getBranchModel();
        if (branchModel.getSubstitutionModels().size() != 1 || branchModel.requiresMatrixConvolution()) {
            return "the adjoint needs a single substitution model without matrix convolution";
        }
        return null;
    }

    @Override
    public double[] getBranchGradient() {

        context.prepareEvaluation();
        final double[] adjoint = (double[]) context.getTreeTrait().getTrait(tree, null);

        final int length = stateCount * stateCount;
        for (int ij = 0; ij < length; ++ij) {
            if (!Double.isFinite(adjoint[ij])) {
                throw new NonFiniteAdjointException("non-finite adjoint entry " + ij + ": " + adjoint[ij]);
            }
        }

        final double[] gradient = new double[wrts.size()];
        for (int d = 0; d < wrts.size(); ++d) {
            WrappedMatrix differential = substitutionModel.getInfinitesimalDifferentialMatrix(wrts.get(d));
            double sum = 0.0;
            for (int i = 0; i < stateCount; ++i) {
                for (int j = 0; j < stateCount; ++j) {
                    sum += adjoint[i * stateCount + j] * differential.get(i, j);
                }
            }
            gradient[d] = sum;
        }
        return gradient;
    }

    @Override
    public String getName() {
        return "adjoint";
    }

    static final class NonFiniteAdjointException extends RuntimeException {
        NonFiniteAdjointException(String message) {
            super(message);
        }
    }
}
