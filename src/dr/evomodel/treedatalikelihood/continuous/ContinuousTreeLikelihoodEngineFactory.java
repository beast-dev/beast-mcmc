/*
 * ContinuousTreeLikelihoodEngineFactory.java
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.canonical.adapter.CanonicalOUGradientAdapter;
import dr.evomodel.treedatalikelihood.continuous.canonical.adapter.CanonicalOUIntegrator;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;

/**
 * Centralizes the light engine-selection policy for continuous trait tree likelihoods.
 */
final class ContinuousTreeLikelihoodEngineFactory {

    private ContinuousTreeLikelihoodEngineFactory() {
    }

    static Selection legacy(final ContinuousDataLikelihoodDelegate delegate) {
        return new Selection(
                new LegacyCdiContinuousTreeLikelihoodEngine(delegate),
                null,
                null,
                false,
                false,
                false);
    }

    static Selection canonical(final ContinuousDataLikelihoodDelegate delegate) {
        final DiffusionProcessDelegate diffusionProcessDelegate = delegate.getDiffusionProcessDelegate();
        if (diffusionProcessDelegate instanceof OUDiffusionModelDelegate) {
            return canonicalOU(delegate, (OUDiffusionModelDelegate) diffusionProcessDelegate);
        }
        if (diffusionProcessDelegate instanceof HomogeneousDiffusionModelDelegate) {
            return cdiBM(delegate);
        }
        throw new UnsupportedOperationException(
                "Canonical continuous likelihood currently supports only homogeneous BM and OU models.");
    }

    static Selection cdiBM(final ContinuousDataLikelihoodDelegate delegate) {
        validateSingleTrait(delegate, "CDI/moment-form BM likelihood");
        if (!(delegate.getDiffusionProcessDelegate() instanceof HomogeneousDiffusionModelDelegate)) {
            throw new UnsupportedOperationException(
                    "CDI/moment-form BM likelihood is only available for homogeneous BM diffusion delegates.");
        }

        return new Selection(
                new CdiBMContinuousTreeLikelihoodEngine(delegate),
                null,
                null,
                true,
                false,
                false);
    }

    static Selection momentBM(final ContinuousDataLikelihoodDelegate delegate) {
        validateSingleTrait(delegate, "Native moment-form BM likelihood");
        if (!(delegate.getDiffusionProcessDelegate() instanceof HomogeneousDiffusionModelDelegate)) {
            throw new UnsupportedOperationException(
                    "Native moment-form BM likelihood is the native moment-form BM backend and supports only homogeneous BM.");
        }
        if (delegate.getPrecisionType() != PrecisionType.SCALAR) {
            throw new UnsupportedOperationException(
                    "Native moment-form BM likelihood currently supports only scalar tip precision.");
        }

        return new Selection(
                new MomentBMContinuousTreeLikelihoodEngine(delegate),
                null,
                null,
                false,
                true,
                false);
    }

    static Selection canonicalOU(final ContinuousDataLikelihoodDelegate delegate,
                                 final OUDiffusionModelDelegate ouDelegate) {
        validateSingleTrait(delegate, "Canonical OU likelihood");
        if (ouDelegate instanceof IntegratedOUDiffusionModelDelegate) {
            throw new UnsupportedOperationException(
                    "Canonical OU likelihood does not yet support integrated OU delegates.");
        }

        final CanonicalOUIntegrator integrator = delegate.createCanonicalOUIntegrator(ouDelegate);
        return new Selection(
                new CanonicalOUContinuousTreeLikelihoodEngine(delegate, integrator),
                integrator,
                new CanonicalOUGradientAdapter(integrator, delegate.getTraitDim()),
                false,
                false,
                true);
    }

    private static void validateSingleTrait(final ContinuousDataLikelihoodDelegate delegate,
                                            final String label) {
        if (delegate.getTraitCount() != 1) {
            throw new UnsupportedOperationException(label + " currently supports only a single trait partition.");
        }
    }

    static final class Selection {

        private final ContinuousTreeLikelihoodEngine engine;
        private final CanonicalOUIntegrator canonicalOUIntegrator;
        private final CanonicalOUGradientAdapter canonicalOUGradientAdapter;
        private final boolean useCdiBMLikelihood;
        private final boolean useMomentBMLikelihood;
        private final boolean useCanonicalOULikelihood;

        private Selection(final ContinuousTreeLikelihoodEngine engine,
                          final CanonicalOUIntegrator canonicalOUIntegrator,
                          final CanonicalOUGradientAdapter canonicalOUGradientAdapter,
                          final boolean useCdiBMLikelihood,
                          final boolean useMomentBMLikelihood,
                          final boolean useCanonicalOULikelihood) {
            this.engine = engine;
            this.canonicalOUIntegrator = canonicalOUIntegrator;
            this.canonicalOUGradientAdapter = canonicalOUGradientAdapter;
            this.useCdiBMLikelihood = useCdiBMLikelihood;
            this.useMomentBMLikelihood = useMomentBMLikelihood;
            this.useCanonicalOULikelihood = useCanonicalOULikelihood;
        }

        ContinuousTreeLikelihoodEngine getEngine() {
            return engine;
        }

        CanonicalOUIntegrator getCanonicalOUIntegrator() {
            return canonicalOUIntegrator;
        }

        CanonicalOUGradientAdapter getCanonicalOUGradientAdapter() {
            return canonicalOUGradientAdapter;
        }

        boolean usesCdiBMLikelihood() {
            return useCdiBMLikelihood;
        }

        boolean usesMomentBMLikelihood() {
            return useMomentBMLikelihood;
        }

        boolean usesCanonicalOULikelihood() {
            return useCanonicalOULikelihood;
        }
    }
}
