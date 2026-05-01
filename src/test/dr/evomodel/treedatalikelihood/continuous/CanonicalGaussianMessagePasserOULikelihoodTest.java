/*
 * CanonicalGaussianMessagePasserOULikelihoodTest.java
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

package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTipObservation;

public class CanonicalGaussianMessagePasserOULikelihoodTest extends CanonicalOUMessagePasserTestSupport {

    public CanonicalGaussianMessagePasserOULikelihoodTest(final String name) {
        super(name);
    }

    public void testLogLikelihoodConjugateRootFullyObservedIsFinite() {
        System.out.println("\nTest: canonical OU log-likelihood is finite (conjugate root, fully observed)");
        final OUSetup setup = buildOUSetup("canonLogL_full", buildFullyObservedTips());
        final double logLikelihood = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);
        assertTrue("Fully observed canonical OU log-likelihood must be finite", Double.isFinite(logLikelihood));
    }

    public void testLogLikelihoodConjugateRootPartiallyObservedIsFinite() {
        System.out.println("\nTest: canonical OU log-likelihood is finite (conjugate root, partially observed)");
        final OUSetup setup = buildOUSetup("canonLogL_partial", buildPartiallyObservedTips());
        final double logLikelihood = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);
        assertTrue("Partially observed canonical OU log-likelihood must be finite", Double.isFinite(logLikelihood));
    }

    public void testLogLikelihoodConjugateRootGaussianLinkFullyObservedIsFinite() {
        System.out.println("\nTest: canonical OU log-likelihood is finite (conjugate root, Gaussian-link tips)");
        final CanonicalTipObservation[] tips = buildFullyObservedTips();
        final OUSetup setup = buildOUSetup("canonLogL_gaussian_link", tips);
        installIdentityGaussianLinkTips(setup, tips, 0.05);

        final double logLikelihood = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);

        assertTrue("Gaussian-link canonical OU log-likelihood must be finite", Double.isFinite(logLikelihood));
    }

    public void testLogLikelihoodConjugateRootGaussianLinkDiffersFromCollapsedIdentity() {
        System.out.println("\nTest: canonical OU Gaussian-link log-likelihood differs from collapsed identity tips");
        final CanonicalTipObservation[] tips = buildFullyObservedTips();
        final OUSetup identitySetup = buildOUSetup("canonLogL_gaussian_identity", tips);
        final OUSetup gaussianSetup = buildOUSetup("canonLogL_gaussian_noise", tips);
        installIdentityGaussianLinkTips(gaussianSetup, tips, 0.25);

        final double identityLogLikelihood =
                identitySetup.passer.computePostOrderLogLikelihood(identitySetup.provider, identitySetup.rootPrior);
        final double gaussianLogLikelihood =
                gaussianSetup.passer.computePostOrderLogLikelihood(gaussianSetup.provider, gaussianSetup.rootPrior);

        assertTrue("Collapsed identity log-likelihood must be finite", Double.isFinite(identityLogLikelihood));
        assertTrue("Gaussian-link log-likelihood must be finite", Double.isFinite(gaussianLogLikelihood));
        assertTrue("Finite observation noise should change the integrated likelihood",
                Math.abs(identityLogLikelihood - gaussianLogLikelihood) > 1.0e-6);
    }

    public void testOrthogonalBlockConjugateRootFullyObservedLikelihoodMatchesDense() {
        System.out.println("\nTest: canonical OU orthogonal-block log-likelihood matches dense (fully observed)");
        checkOrthogonalBlockLikelihoodMatchesDense(
                buildOrthogonalBlockOUSetup("canonOrthLike_full", buildFullyObservedTips()),
                "fully observed");
    }

    public void testOrthogonalBlockConjugateRootPartiallyObservedLikelihoodMatchesDense() {
        System.out.println("\nTest: canonical OU orthogonal-block log-likelihood matches dense (partially observed)");
        checkOrthogonalBlockLikelihoodMatchesDense(
                buildOrthogonalBlockOUSetup("canonOrthLike_partial", buildPartiallyObservedTips()),
                "partially observed");
    }
}
