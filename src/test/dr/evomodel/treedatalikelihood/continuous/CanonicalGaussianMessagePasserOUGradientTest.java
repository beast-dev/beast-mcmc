/*
 * CanonicalGaussianMessagePasserOUGradientTest.java
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

import dr.evomodel.treedatalikelihood.continuous.canonical.adapter.CanonicalOUTreeLikelihoodIntegrator;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalPreparedBranchSnapshot;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalRootPrior;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipObservation;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTransitionCachePhases;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTreeMessagePasser;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockPreparedBranchBasis;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockPreparedBranchHandle;

/**
 * Validates canonical OU tree-passer gradients against numerical finite differences.
 */
public class CanonicalGaussianMessagePasserOUGradientTest extends CanonicalOUMessagePasserTestSupport {

    private static final String TRANSITION_CACHE_DIAGNOSTICS_PROPERTY =
            "beast.debug.canonicalTransitionCache";

    public CanonicalGaussianMessagePasserOUGradientTest(final String name) {
        super(name);
    }

    public void testGradientQConjugateRootFullyObservedMatchesNumerical() {
        System.out.println("\nTest: canonical OU ∂logL/∂Q analytic vs numerical (conjugate root, fully observed)");
        checkGradientQ(buildOUSetup("canonGradQ_full", buildFullyObservedTips()), "fully observed");
    }

    public void testGradientAConjugateRootFullyObservedMatchesNumerical() {
        System.out.println("\nTest: canonical OU ∂logL/∂A analytic vs numerical (conjugate root, fully observed)");
        checkGradientA(buildOUSetup("canonGradA_full", buildFullyObservedTips()), "fully observed");
    }

    public void testGradientMuConjugateRootFullyObservedMatchesNumerical() {
        System.out.println("\nTest: canonical OU ∂logL/∂mu analytic vs numerical (conjugate root, fully observed)");
        checkGradientMu(buildOUSetup("canonGradMu_full", buildFullyObservedTips()), "fully observed");
    }

    public void testGradientBranchLengthsConjugateRootFullyObservedMatchesNumerical() {
        System.out.println("\nTest: canonical OU ∂logL/∂t analytic vs numerical (conjugate root, fully observed)");
        checkGradientBranchLengths(buildOUSetup("canonGradT_full", buildFullyObservedTips()), "fully observed");
    }

    public void testGradientBranchLengthsConjugateRootGaussianLinkMatchesNumerical() {
        System.out.println("\nTest: canonical OU ∂logL/∂t analytic vs numerical (conjugate root, Gaussian-link tips)");
        final CanonicalTipObservation[] tips = buildFullyObservedTips();
        final OUSetup setup = buildOUSetup("canonGradT_gaussian_link", tips);
        installIdentityGaussianLinkTips(setup, tips, 0.05);
        checkGradientBranchLengths(setup, "Gaussian-link tips");
    }

    public void testGradientQConjugateRootPartiallyObservedMatchesNumerical() {
        System.out.println("\nTest: canonical OU ∂logL/∂Q analytic vs numerical (conjugate root, partially observed)");
        checkGradientQ(buildOUSetup("canonGradQ_partial", buildPartiallyObservedTips()), "partially observed");
    }

    public void testGradientQFixedRootPartiallyObservedMatchesNumerical() {
        System.out.println("\nTest: canonical OU ∂logL/∂Q analytic vs numerical (fixed root, partially observed)");
        checkGradientQ(buildOUSetup("canonGradQ_fixed_partial", buildPartiallyObservedTips(), rootPriorInf),
                "fixed root, partially observed");
    }

    public void testGradientAConjugateRootPartiallyObservedMatchesNumerical() {
        System.out.println("\nTest: canonical OU ∂logL/∂A analytic vs numerical (conjugate root, partially observed)");
        checkGradientA(buildOUSetup("canonGradA_partial", buildPartiallyObservedTips()), "partially observed");
    }

    public void testGradientMuConjugateRootPartiallyObservedMatchesNumerical() {
        System.out.println("\nTest: canonical OU ∂logL/∂mu analytic vs numerical (conjugate root, partially observed)");
        checkGradientMu(buildOUSetup("canonGradMu_partial", buildPartiallyObservedTips()), "partially observed");
    }

    public void testGradientRootMeanConjugateRootPartiallyObservedMatchesNumerical() {
        System.out.println("\nTest: canonical OU ∂logL/∂rootMean analytic vs numerical (conjugate root, partially observed)");

        final OUSetup setup = buildOUSetup("canonGradRootMean_partial", buildPartiallyObservedTips());
        setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);
        final double[] traitCovariance = new double[dimTrait * dimTrait];
        setup.provider.fillTraitCovarianceFlat(traitCovariance);
        final double[] gradAnalytic = new double[dimTrait];
        setup.rootPrior.accumulateRootMeanGradient(
                setup.passer.getPostOrderState(treeModel.getRoot().getNumber()),
                traitCovariance,
                gradAnalytic);

        final double[] gradNumeric = new double[dimTrait];
        final dr.inference.model.Parameter rootMean = rootPrior.getMeanParameter();
        for (int i = 0; i < dimTrait; i++) {
            final double original = rootMean.getParameterValue(i);

            rootMean.setParameterValue(i, original + FD_DELTA);
            final double logLPlus = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);

            rootMean.setParameterValue(i, original - FD_DELTA);
            final double logLMinus = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);

            rootMean.setParameterValue(i, original);
            gradNumeric[i] = (logLPlus - logLMinus) / (2.0 * FD_DELTA);
        }

        System.out.printf("  Analytic ∂logL/∂rootMean: %s%n", java.util.Arrays.toString(gradAnalytic));
        System.out.printf("  Numeric  ∂logL/∂rootMean: %s%n", java.util.Arrays.toString(gradNumeric));
        for (int i = 0; i < dimTrait; i++) {
            assertEquals("∂logL/∂rootMean[" + i + "] (partially observed)",
                    gradNumeric[i], gradAnalytic[i], TOL);
        }
    }

    public void testGradientBranchLengthsConjugateRootPartiallyObservedMatchesNumerical() {
        System.out.println("\nTest: canonical OU ∂logL/∂t analytic vs numerical (conjugate root, partially observed)");
        checkGradientBranchLengths(buildOUSetup("canonGradT_partial", buildPartiallyObservedTips()), "partially observed");
    }

    public void testOrthogonalBlockConjugateRootFullyObservedNativeSelectionGradientMatchesNumerical() {
        System.out.println("\nTest: canonical OU orthogonal-block native ∂logL/∂A matches numerical (fully observed)");
        checkOrthogonalBlockNativeSelectionGradient(
                buildOrthogonalBlockOUSetup("canonOrthGrad_full", buildFullyObservedTips()),
                "fully observed");
    }

    public void testOrthogonalBlockConjugateRootPartiallyObservedNativeSelectionGradientMatchesNumerical() {
        System.out.println("\nTest: canonical OU orthogonal-block native ∂logL/∂A matches numerical (partially observed)");
        checkOrthogonalBlockNativeSelectionGradient(
                buildOrthogonalBlockOUSetup("canonOrthGrad_partial", buildPartiallyObservedTips()),
                "partially observed");
    }

    public void testOrthogonalBlockConjugateRootPartiallyObservedGradientQMatchesDense() {
        System.out.println("\nTest: canonical OU orthogonal-block ∂logL/∂Q matches dense (partially observed)");
        checkOrthogonalBlockGradientQMatchesDense(
                buildOrthogonalBlockOUSetup("canonOrthGradQ_partial", buildPartiallyObservedTips()),
                "partially observed");
    }

    public void testStandaloneIntegratorReusesOneJointGradientPassAcrossAQMu() {
        System.out.println("\nTest: canonical OU integrator reuses one joint A/Q/mu gradient pass");

        final OUSetup setup = buildOUSetup("canonJointGradReuse_partial", buildPartiallyObservedTips());
        final CountingPasser countingPasser = new CountingPasser(setup.passer);
        final CanonicalOUTreeLikelihoodIntegrator integrator =
                new CanonicalOUTreeLikelihoodIntegrator(countingPasser, setup.provider, setup.rootPrior);

        final double[] jointA = new double[dimTrait * dimTrait];
        final double[] jointQ = new double[dimTrait * dimTrait];
        final double[] jointMu = new double[dimTrait];
        integrator.computeSelectionGradient(jointA);
        integrator.computeDiffusionGradient(jointQ);
        integrator.computeStationaryMeanGradient(jointMu);

        assertEquals("joint post-order reuse", 1, countingPasser.postOrderCalls);
        assertEquals("joint pre-order reuse", 1, countingPasser.preOrderCalls);
        assertEquals("joint branch-adjoint reuse", 1, countingPasser.jointGradientCalls);

        refreshMessages(setup);
        final double[] directA = new double[dimTrait * dimTrait];
        final double[] directQ = new double[dimTrait * dimTrait];
        final double[] directMu = new double[dimTrait];
        setup.passer.computeJointGradients(setup.provider, directA, directQ, directMu);

        assertArrayEquals("cached vs direct joint A", directA, jointA, 1.0e-11);
        assertArrayEquals("cached vs direct joint Q", directQ, jointQ, 1.0e-11);
        assertArrayEquals("cached vs direct joint mu", directMu, jointMu, 1.0e-11);
    }

    public void testStandaloneIntegratorModelDirtyInvalidatesJointGradientCache() {
        System.out.println("\nTest: canonical OU integrator model dirty invalidates joint gradient cache");

        final OUSetup setup = buildOUSetup("canonJointGradModelDirty_partial", buildPartiallyObservedTips());
        final CountingPasser countingPasser = new CountingPasser(setup.passer);
        final CanonicalOUTreeLikelihoodIntegrator integrator =
                new CanonicalOUTreeLikelihoodIntegrator(countingPasser, setup.provider, setup.rootPrior);

        integrator.computeSelectionGradient(new double[dimTrait * dimTrait]);
        integrator.computeDiffusionGradient(new double[dimTrait * dimTrait]);
        assertEquals("initial joint gradient computation", 1, countingPasser.jointGradientCalls);

        integrator.markModelDirty();
        integrator.computeStationaryMeanGradient(new double[dimTrait]);
        assertEquals("joint gradient recomputed after model dirty", 2, countingPasser.jointGradientCalls);
    }

    public void testStandaloneIntegratorReloadedTipsInvalidateJointGradientCache() {
        System.out.println("\nTest: canonical OU integrator tip reload invalidates joint gradient cache");

        final OUSetup setup = buildOUSetup("canonJointGradTipReload_partial", buildPartiallyObservedTips());
        final CountingPasser countingPasser = new CountingPasser(setup.passer);
        final CanonicalOUTreeLikelihoodIntegrator integrator =
                new CanonicalOUTreeLikelihoodIntegrator(countingPasser, setup.provider, setup.rootPrior);

        integrator.computeSelectionGradient(new double[dimTrait * dimTrait]);
        assertEquals("initial joint gradient computation", 1, countingPasser.jointGradientCalls);

        integrator.reloadTips(buildFullyObservedTips());
        integrator.computeDiffusionGradient(new double[dimTrait * dimTrait]);
        assertEquals("joint gradient recomputed after tip reload", 2, countingPasser.jointGradientCalls);
    }

    public void testStandaloneIntegratorStoreRestorePreservesCachedJointGradients() {
        System.out.println("\nTest: canonical OU integrator store/restore preserves cached joint gradients");

        final OUSetup setup = buildOUSetup("canonJointGradStoreRestore_partial", buildPartiallyObservedTips());
        final CountingPasser countingPasser = new CountingPasser(setup.passer);
        final CanonicalOUTreeLikelihoodIntegrator integrator =
                new CanonicalOUTreeLikelihoodIntegrator(countingPasser, setup.provider, setup.rootPrior);

        integrator.computeSelectionGradient(new double[dimTrait * dimTrait]);
        assertEquals("initial joint gradient computation", 1, countingPasser.jointGradientCalls);

        integrator.storeState();
        integrator.markModelDirty();
        integrator.computeDiffusionGradient(new double[dimTrait * dimTrait]);
        assertEquals("joint gradient recomputed after model dirty", 2, countingPasser.jointGradientCalls);

        integrator.restoreState();
        integrator.computeStationaryMeanGradient(new double[dimTrait]);
        assertEquals("restored cached joint gradient reused", 2, countingPasser.jointGradientCalls);
    }

    public void testOrthogonalBlockJointGradientsMatchAcrossBranchGradientThreadCounts() {
        System.out.println(
                "\nTest: canonical OU orthogonal-block joint gradients match across branch-gradient thread counts");

        final CanonicalTipObservation[] tips = buildPartiallyObservedTips();
        final OrthogonalBlockOUSetup baseSetup =
                buildOrthogonalBlockOUSetup("canonOrthJointGradThreads_partial", tips);
        final OUSetup sequentialSetup = withPasserParallelism(baseSetup.orthogonal, tips, 1);
        final OUSetup parallelSetup = withPasserParallelism(baseSetup.orthogonal, tips, 3);

        refreshMessages(sequentialSetup);
        refreshMessages(parallelSetup);

        final int nativeSelectionDim = baseSetup.nativeParameter.getDimension();
        final int diffusionDim = dimTrait * dimTrait;
        final double[] sequentialA = new double[nativeSelectionDim];
        final double[] sequentialQ = new double[diffusionDim];
        final double[] sequentialMu = new double[dimTrait];
        final double[] parallelA = new double[nativeSelectionDim];
        final double[] parallelQ = new double[diffusionDim];
        final double[] parallelMu = new double[dimTrait];

        sequentialSetup.passer.computeJointGradients(
                sequentialSetup.provider,
                sequentialA,
                sequentialQ,
                sequentialMu);
        parallelSetup.passer.computeJointGradients(
                parallelSetup.provider,
                parallelA,
                parallelQ,
                parallelMu);

        assertArrayEquals("orthogonal joint A threads", sequentialA, parallelA, 1.0e-11);
        assertArrayEquals("orthogonal joint Q threads", sequentialQ, parallelQ, 1.0e-11);
        assertArrayEquals("orthogonal joint mu threads", sequentialMu, parallelMu, 1.0e-11);
    }

    public void testOrthogonalBlockPreparedBranchSnapshotReusedUntilBranchLengthChanges() {
        System.out.println("\nTest: orthogonal-block prepared branch snapshot cache reuse and branch-length invalidation");

        final OrthogonalBlockOUSetup setup =
                buildOrthogonalBlockOUSetup("canonOrthSnapshotCache_partial", buildPartiallyObservedTips());
        final int childNodeIndex = firstNonRootNodeIndex();

        final CanonicalPreparedBranchSnapshot first =
                setup.orthogonal.provider.getPreparedBranchSnapshot(childNodeIndex);
        final CanonicalPreparedBranchSnapshot second =
                setup.orthogonal.provider.getPreparedBranchSnapshot(childNodeIndex);

        assertSame("prepared branch snapshot reused on cache hit", first, second);
        assertSame("orthogonal prepared basis reused on cache hit",
                orthogonalPreparedBasis(first), orthogonalPreparedBasis(second));

        setup.orthogonal.branchScale[childNodeIndex] += 0.125;
        final CanonicalPreparedBranchSnapshot afterLengthChange =
                setup.orthogonal.provider.getPreparedBranchSnapshot(childNodeIndex);

        assertNotSame("branch length change refreshes snapshot", first, afterLengthChange);
        assertSame("prepared basis buffer is retained across snapshot refresh",
                orthogonalPreparedBasis(first), orthogonalPreparedBasis(afterLengthChange));
    }

    public void testJointGradientsReuseTransitionCacheWithoutRebuilds() {
        System.out.println("\nTest: canonical OU joint gradients reuse prepared transition cache");

        final String previous = System.getProperty(TRANSITION_CACHE_DIAGNOSTICS_PROPERTY);
        System.setProperty(TRANSITION_CACHE_DIAGNOSTICS_PROPERTY, "true");
        try {
            final OrthogonalBlockOUSetup setup =
                    buildOrthogonalBlockOUSetup("canonOrthGradientCacheCounters", buildPartiallyObservedTips());
            refreshMessages(setup.orthogonal);

            final double[] gradA = new double[setup.nativeParameter.getDimension()];
            final double[] gradQ = new double[dimTrait * dimTrait];
            final double[] gradMu = new double[dimTrait];

            setup.orthogonal.passer.computeJointGradients(
                    setup.orthogonal.provider,
                    gradA,
                    gradQ,
                    gradMu);
            final long requestsAfterFirst =
                    setup.orthogonal.provider.getTransitionCacheRequestCount(
                            CanonicalTransitionCachePhases.GRADIENT_PREP);
            final long transitionRebuildsAfterFirst =
                    setup.orthogonal.provider.getTransitionCacheRebuildCount();
            final long preparedRebuildsAfterFirst =
                    setup.orthogonal.provider.getPreparedBranchRebuildCount();

            assertTrue("gradient prep records cache requests", requestsAfterFirst > 0L);
            assertTrue("initial traversal prepared branch bases", preparedRebuildsAfterFirst > 0L);

            setup.orthogonal.passer.computeJointGradients(
                    setup.orthogonal.provider,
                    gradA,
                    gradQ,
                    gradMu);

            assertTrue("second gradient pass records additional cache requests",
                    setup.orthogonal.provider.getTransitionCacheRequestCount(
                            CanonicalTransitionCachePhases.GRADIENT_PREP)
                            > requestsAfterFirst);
            assertEquals("cached gradient prep does not rebuild transitions",
                    transitionRebuildsAfterFirst,
                    setup.orthogonal.provider.getTransitionCacheRebuildCount());
            assertEquals("cached gradient prep does not rebuild prepared branch bases",
                    preparedRebuildsAfterFirst,
                    setup.orthogonal.provider.getPreparedBranchRebuildCount());
        } finally {
            if (previous == null) {
                System.clearProperty(TRANSITION_CACHE_DIAGNOSTICS_PROPERTY);
            } else {
                System.setProperty(TRANSITION_CACHE_DIAGNOSTICS_PROPERTY, previous);
            }
        }
    }

    private OrthogonalBlockPreparedBranchBasis orthogonalPreparedBasis(
            final CanonicalPreparedBranchSnapshot snapshot) {
        assertTrue("expected an orthogonal-block prepared branch handle",
                snapshot.getPreparedBranchHandle() instanceof OrthogonalBlockPreparedBranchHandle);
        return ((OrthogonalBlockPreparedBranchHandle) snapshot.getPreparedBranchHandle()).getBasis();
    }

    private int firstNonRootNodeIndex() {
        final int rootIndex = treeModel.getRoot().getNumber();
        for (int nodeIndex = 0; nodeIndex < treeModel.getNodeCount(); ++nodeIndex) {
            if (nodeIndex != rootIndex) {
                return nodeIndex;
            }
        }
        throw new IllegalStateException("Tree has no non-root node");
    }

    private void assertArrayEquals(final String label,
                                   final double[] expected,
                                   final double[] actual,
                                   final double tolerance) {
        assertEquals(label + " length", expected.length, actual.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals(label + "[" + i + "]", expected[i], actual[i], tolerance);
        }
    }

    private static final class CountingPasser implements CanonicalTreeMessagePasser {
        private final CanonicalTreeMessagePasser delegate;
        private int postOrderCalls = 0;
        private int preOrderCalls = 0;
        private int jointGradientCalls = 0;

        private CountingPasser(final CanonicalTreeMessagePasser delegate) {
            this.delegate = delegate;
        }

        @Override
        public int getDimension() {
            return delegate.getDimension();
        }

        @Override
        public int getTipCount() {
            return delegate.getTipCount();
        }

        @Override
        public void setTipObservation(final int tipIndex, final CanonicalTipObservation observation) {
            delegate.setTipObservation(tipIndex, observation);
        }

        @Override
        public double computePostOrderLogLikelihood(final CanonicalBranchTransitionProvider transitionProvider,
                                                    final CanonicalRootPrior rootPrior) {
            postOrderCalls++;
            return delegate.computePostOrderLogLikelihood(transitionProvider, rootPrior);
        }

        @Override
        public CanonicalGaussianState getPostOrderState(final int nodeIndex) {
            return delegate.getPostOrderState(nodeIndex);
        }

        @Override
        public void computePreOrder(final CanonicalBranchTransitionProvider transitionProvider,
                                    final CanonicalRootPrior rootPrior) {
            preOrderCalls++;
            delegate.computePreOrder(transitionProvider, rootPrior);
        }

        @Override
        public void computeGradientBranchLengths(final CanonicalBranchTransitionProvider transitionProvider,
                                                 final double[] gradT) {
            delegate.computeGradientBranchLengths(transitionProvider, gradT);
        }

        @Override
        public void computeJointGradients(final CanonicalBranchTransitionProvider transitionProvider,
                                          final double[] gradA,
                                          final double[] gradQ,
                                          final double[] gradMu) {
            jointGradientCalls++;
            delegate.computeJointGradients(transitionProvider, gradA, gradQ, gradMu);
        }

        @Override
        public void storeState() {
            delegate.storeState();
        }

        @Override
        public void restoreState() {
            delegate.restoreState();
        }

        @Override
        public void acceptState() {
            delegate.acceptState();
        }
    }
}
