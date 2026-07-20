/*
 * CanonicalOUBranchRateGradient.java
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

package dr.evomodel.treedatalikelihood.continuousNew;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.BranchRateGradient;

import dr.inference.model.Parameter;

/**
 * Branch-rate gradient for the canonical OU tree likelihood path.
 *
 * <p>The canonical message passer gives derivatives with respect to effective
 * branch time. This adapter converts those into the existing BEAST branch-rate
 * parameter convention while preserving the old {@code BranchRateGradient}
 * XML-facing type.</p>
 */
public final class CanonicalOUBranchRateGradient extends BranchRateGradient {

    private final CtsCompatibility likelihoodDelegate;

    public CanonicalOUBranchRateGradient(final TreeDataLikelihood treeDataLikelihood,
                                         final CtsCompatibility likelihoodDelegate,
                                         final Parameter rateParameter) {
        super(treeDataLikelihood, rateParameter);
        if (!likelihoodDelegate.usesCanonicalOULikelihood()) {
            throw new IllegalArgumentException(
                    "Canonical OU branch-rate gradient requires canonical OU likelihood mode.");
        }
        if (branchRateModel == null) {
            throw new IllegalArgumentException(
                    "Canonical OU branch-rate gradient requires a differentiable branch-rate model.");
        }
        if (rateParameter == null) {
            throw new IllegalArgumentException(
                    "Canonical OU branch-rate gradient requires a branch-rate model with an explicit rate parameter.");
        }
        this.likelihoodDelegate = likelihoodDelegate;
    }

    @Override
    public double[] getGradientLogDensity() {
        treeDataLikelihood.makeDirty();
        final double[] effectiveTimeGradient =
                likelihoodDelegate.getCanonicalOUGradientAdapter()
                        .getEffectiveBranchTimeGradient(tree.getNodeCount());
        return CanonicalOUBranchTimeGradient.branchRateParameterGradient(
                tree,
                branchRateModel,
                likelihoodDelegate.getRateTransformation(),
                effectiveTimeGradient);
    }

    @Override
    public double[] getDiagonalHessianLogDensity() {
        throw new UnsupportedOperationException(
                "Canonical OU branch-rate Hessians are not implemented.");
    }

    @Override
    public String getReport() {
        double[] testGradient = getNumericalGradient();

        StringBuilder sb = new StringBuilder();
        sb.append("Canonical OU branch-rate gradient: ")
                .append(new dr.math.matrixAlgebra.Vector(getGradientLogDensity()));
        sb.append("\n");
        sb.append("numeric: ").append(new dr.math.matrixAlgebra.Vector(testGradient));
        sb.append("\n");

        return sb.toString();
    }
}
