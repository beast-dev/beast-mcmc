/*
 * CanonicalConjugateRootPriorAdapter.java
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

package dr.evomodel.treedatalikelihood.continuous.adapter;

import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalRootPrior;
import dr.inference.timeseries.engine.gaussian.CanonicalGaussianMessageOps;
import dr.inference.timeseries.representation.CanonicalGaussianState;
import dr.inference.timeseries.representation.CanonicalGaussianUtils;

/**
 * Canonical-form adapter over {@link ConjugateRootTraitPrior}.
 *
 * <p>This class keeps the canonical OU pathway rooted in precision/information
 * form without crossing back through a moment-form root-likelihood bridge. Root
 * seeding for fixed-root and diffuse-root cases is intentionally left unsupported
 * in this first cut because they require exact delta/flat canonical states rather
 * than finite normalized Gaussian states.
 */
public final class CanonicalConjugateRootPriorAdapter implements CanonicalRootPrior {

    private final ConjugateRootTraitPrior prior;
    private final int dim;

    private final double[] scratchMean;
    private final double[][] scratchCovariance;
    private final CanonicalGaussianState scratchPriorState;
    private final CanonicalGaussianState scratchCombinedState;
    private final CanonicalGaussianMessageOps.Workspace workspace;

    public CanonicalConjugateRootPriorAdapter(final ConjugateRootTraitPrior prior, final int dim) {
        this.prior = prior;
        this.dim = dim;
        this.scratchMean = new double[dim];
        this.scratchCovariance = new double[dim][dim];
        this.scratchPriorState = new CanonicalGaussianState(dim);
        this.scratchCombinedState = new CanonicalGaussianState(dim);
        this.workspace = new CanonicalGaussianMessageOps.Workspace(dim);
    }

    @Override
    public int getDimension() {
        return dim;
    }

    @Override
    public void fillRootPriorState(final double[][] traitCovariance, final CanonicalGaussianState out) {
        final double kappa0 = prior.getPseudoObservations();
        if (!(kappa0 > 0.0) || Double.isInfinite(kappa0)) {
            throw new UnsupportedOperationException(
                    "Canonical root seeding currently supports only finite positive pseudo-observations.");
        }

        final double[] mean = prior.getMean();
        final double[][] scaledCovariance = scratchCovariance;
        final double invKappa0 = 1.0 / kappa0;
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                scaledCovariance[i][j] = traitCovariance[i][j] * invKappa0;
            }
            scratchMean[i] = mean[i];
        }
        CanonicalGaussianUtils.fillStateFromMoments(scratchMean, scaledCovariance, out);
    }

    @Override
    public double computeLogMarginalLikelihood(final CanonicalGaussianState rootMessage,
                                               final double[][] traitCovariance) {
        if (isFixedRoot()) {
            final double[] mean = prior.getMean();
            double quadratic = 0.0;
            double linear = 0.0;
            for (int i = 0; i < dim; i++) {
                linear += rootMessage.information[i] * mean[i];
                for (int j = 0; j < dim; j++) {
                    quadratic += mean[i] * rootMessage.precision[i][j] * mean[j];
                }
            }
            return -0.5 * quadratic + linear - rootMessage.logNormalizer;
        }

        fillRootPriorState(traitCovariance, scratchPriorState);
        CanonicalGaussianMessageOps.combineStates(rootMessage, scratchPriorState, scratchCombinedState);
        return CanonicalGaussianMessageOps.normalizationShift(scratchCombinedState, workspace);
    }

    @Override
    public boolean isFixedRoot() {
        return Double.isInfinite(prior.getPseudoObservations());
    }

    @Override
    public void fillFixedRootValue(final double[] out) {
        final double[] mean = prior.getMean();
        System.arraycopy(mean, 0, out, 0, dim);
    }

    @Override
    public void accumulateRootMeanGradient(final CanonicalGaussianState rootMessage,
                                           final double[][] traitCovariance,
                                           final double[] gradOut) {
        if (isFixedRoot()) {
            final double[] mean = prior.getMean();
            for (int i = 0; i < dim; i++) {
                double gradient = rootMessage.information[i];
                for (int j = 0; j < dim; j++) {
                    gradient -= rootMessage.precision[i][j] * mean[j];
                }
                gradOut[i] += gradient;
            }
            return;
        }

        fillRootPriorState(traitCovariance, scratchPriorState);
        CanonicalGaussianMessageOps.combineStates(rootMessage, scratchPriorState, scratchCombinedState);
        CanonicalGaussianUtils.fillMomentsFromCanonical(scratchCombinedState, scratchMean, scratchCovariance);

        final double[] mean = prior.getMean();
        for (int i = 0; i < dim; i++) {
            double delta = scratchMean[i] - mean[i];
            double gradient = 0.0;
            for (int j = 0; j < dim; j++) {
                gradient += scratchPriorState.precision[i][j] * (scratchMean[j] - mean[j]);
            }
            gradOut[i] += gradient;
        }
    }

    @Override
    public double getDiffusionScale() {
        final double kappa0 = prior.getPseudoObservations();
        if (!(kappa0 > 0.0) || Double.isInfinite(kappa0)) {
            return 0.0;
        }
        return 1.0 / kappa0;
    }
}
