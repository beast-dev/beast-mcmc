/*
 * CanonicalOUFocusedParityAndCacheTest.java
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

import dr.evomodel.treedatalikelihood.continuous.canonical.gradient.BranchGradientInputs;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipObservation;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTransitionCachePhases;

import java.util.Arrays;

public class CanonicalOUFocusedParityAndCacheTest extends CanonicalOUMessagePasserTestSupport {

    private static final String TRANSITION_CACHE_DIAGNOSTICS_PROPERTY =
            "beast.debug.canonicalTransitionCache";
    private static final int BRANCH_GRADIENT_BENCHMARK_REPEATS = 32;

    public CanonicalOUFocusedParityAndCacheTest(final String name) {
        super(name);
    }

    public void testDenseAndOrthogonalBlockLikelihoodParityFocused() {
        System.out.println("\nTest: focused dense vs orthogonal-block canonical OU likelihood parity");

        checkOrthogonalBlockLikelihoodMatchesDense(
                buildOrthogonalBlockOUSetup("focusedParityLike_full", buildFullyObservedTips()),
                "fully observed focused suite");
        checkOrthogonalBlockLikelihoodMatchesDense(
                buildOrthogonalBlockOUSetup("focusedParityLike_partial", buildPartiallyObservedTips()),
                "partially observed focused suite");
    }

    public void testDenseAndOrthogonalBlockGradientParityFocused() {
        System.out.println("\nTest: focused dense vs orthogonal-block canonical OU gradient parity");

        final OrthogonalBlockOUSetup setup =
                buildOrthogonalBlockOUSetup("focusedParityGrad_partial", buildPartiallyObservedTips());

        checkOrthogonalBlockGradientQMatchesDense(setup, "partially observed focused suite");
        checkMeanGradientMatchesDense(setup, "partially observed focused suite");
        checkNativeSelectionGradientMatchesDenseFiniteDifference(setup, "partially observed focused suite");
    }

    public void testTransitionCacheInvalidationCountersAreReasonSpecific() {
        System.out.println("\nTest: canonical OU transition-cache invalidation counters are reason-specific");

        withTransitionCacheDiagnostics(new CheckedRunnable() {
            @Override
            public void run() {
                final OUSetup setup = buildOUSetup("focusedCacheInvalidation", buildPartiallyObservedTips());
                refreshMessages(setup);

                final int childNodeIndex = firstNonRootNodeIndex();
                setup.provider.getPreparedBranchSnapshot(childNodeIndex);

                final long selectionClears = clearCount(setup, "SELECTION_CHANGED");
                setup.aMatrix.setParameterValue(0, 0, setup.aMatrix.getParameterValue(0, 0) + 0.03125);
                setup.provider.getPreparedBranchSnapshot(childNodeIndex);
                assertEquals("selection change invalidates transition cache",
                        selectionClears + 1L, clearCount(setup, "SELECTION_CHANGED"));

                final long diffusionClears = clearCount(setup, "DIFFUSION_CHANGED");
                setup.precisionMatrix.setParameterValue(0, 0,
                        setup.precisionMatrix.getParameterValue(0, 0) * 1.015);
                setup.provider.getPreparedBranchSnapshot(childNodeIndex);
                assertEquals("diffusion change invalidates transition cache",
                        diffusionClears + 1L, clearCount(setup, "DIFFUSION_CHANGED"));

                final long meanClears = clearCount(setup, "STATIONARY_MEAN_CHANGED");
                setup.muParam.setParameterValue(0, setup.muParam.getParameterValue(0) + 0.0625);
                setup.provider.getPreparedBranchSnapshot(childNodeIndex);
                assertEquals("stationary mean change invalidates transition cache",
                        meanClears + 1L, clearCount(setup, "STATIONARY_MEAN_CHANGED"));

                final long branchLengthClears = clearCount(setup, "BRANCH_LENGTH_CHANGED");
                setup.branchScale[childNodeIndex] += 0.125;
                setup.provider.getPreparedBranchSnapshot(childNodeIndex);
                assertEquals("branch length change invalidates branch transition cache entry",
                        branchLengthClears + 1L, clearCount(setup, "BRANCH_LENGTH_CHANGED"));

                refreshMessages(setup);
                final long rebuildsBeforeTipReload = setup.provider.getTransitionCacheRebuildCount();
                setup.passer.setTipObservation(0, observedTip(-1.25, 2.25, 3.25));
                setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);
                setup.passer.computePreOrder(setup.provider, setup.rootPrior);
                assertEquals("tip changes do not rebuild transition cache",
                        rebuildsBeforeTipReload, setup.provider.getTransitionCacheRebuildCount());
            }
        });
    }

    public void testNoTransitionCacheMissesDuringGradientPrepFocused() {
        System.out.println("\nTest: canonical OU gradient-prep phase has no transition-cache misses");

        withTransitionCacheDiagnostics(new CheckedRunnable() {
            @Override
            public void run() {
                final OrthogonalBlockOUSetup setup =
                        buildOrthogonalBlockOUSetup("focusedGradientPrepCache", buildPartiallyObservedTips());
                refreshMessages(setup.orthogonal);

                final double[] gradA = new double[setup.nativeParameter.getDimension()];
                final double[] gradQ = new double[dimTrait * dimTrait];
                final double[] gradMu = new double[dimTrait];
                setup.orthogonal.passer.computeJointGradients(
                        setup.orthogonal.provider,
                        gradA,
                        gradQ,
                        gradMu);

                assertTrue("gradient prep records transition-cache requests",
                        setup.orthogonal.provider.getTransitionCacheRequestCount(
                                CanonicalTransitionCachePhases.GRADIENT_PREP) > 0L);
                assertEquals("gradient prep should reuse prepared traversal transitions",
                        0L,
                        setup.orthogonal.provider.getTransitionCacheMissCount(
                                CanonicalTransitionCachePhases.GRADIENT_PREP));
            }
        });
    }

    public void testDenseGradientRestoreMatchesColdWorkspaceAfterSelectionPerturbation() {
        System.out.println("\nTest: dense canonical OU gradient cache is fresh after selection restore");

        final OUSetup setup = buildOUSetup("focusedDenseRestoreCache", buildPartiallyObservedTips());
        refreshMessages(setup);

        final double[] baselineA = new double[dimTrait * dimTrait];
        final double[] baselineQ = new double[dimTrait * dimTrait];
        final double[] baselineMu = new double[dimTrait];
        setup.passer.computeJointGradients(setup.provider, baselineA, baselineQ, baselineMu);

        setup.aMatrix.storeParameterValues();
        try {
            setup.aMatrix.setParameterValue(0, 1, setup.aMatrix.getParameterValue(0, 1) + 0.0625);
            refreshMessages(setup);
            setup.passer.computeJointGradients(
                    setup.provider,
                    new double[dimTrait * dimTrait],
                    new double[dimTrait * dimTrait],
                    new double[dimTrait]);

            setup.aMatrix.restoreParameterValues();
            setup.aMatrix.fireParameterChangedEvent();
            refreshMessages(setup);

            final double[] restoredA = new double[dimTrait * dimTrait];
            final double[] restoredQ = new double[dimTrait * dimTrait];
            final double[] restoredMu = new double[dimTrait];
            setup.passer.computeJointGradients(setup.provider, restoredA, restoredQ, restoredMu);

            final OUSetup cold = buildOUSetup("focusedDenseRestoreCacheCold", buildPartiallyObservedTips());
            refreshMessages(cold);
            final double[] coldA = new double[dimTrait * dimTrait];
            final double[] coldQ = new double[dimTrait * dimTrait];
            final double[] coldMu = new double[dimTrait];
            cold.passer.computeJointGradients(cold.provider, coldA, coldQ, coldMu);

            assertVectorEquals("baseline vs restored A", baselineA, restoredA, 1.0e-11);
            assertVectorEquals("baseline vs restored Q", baselineQ, restoredQ, 1.0e-11);
            assertVectorEquals("baseline vs restored mu", baselineMu, restoredMu, 1.0e-11);
            assertVectorEquals("cold vs restored A", coldA, restoredA, 1.0e-11);
            assertVectorEquals("cold vs restored Q", coldQ, restoredQ, 1.0e-11);
            assertVectorEquals("cold vs restored mu", coldMu, restoredMu, 1.0e-11);
        } finally {
            setup.aMatrix.acceptParameterValues();
        }
    }

    public void testDenseGradientSharedWorkspaceHandlesDistinctBranchLengths() {
        System.out.println("\nTest: dense canonical OU shared gradient workspace handles distinct branch lengths");

        final CanonicalTipObservation[] tips = buildPartiallyObservedTips();
        final OUSetup base = buildOUSetup("focusedDenseSharedBranchLengths", tips);
        final int rootIndex = treeModel.getRoot().getNumber();
        for (int nodeIndex = 0; nodeIndex < treeModel.getNodeCount(); ++nodeIndex) {
            if (nodeIndex != rootIndex) {
                base.branchScale[nodeIndex] = 0.75 + 0.125 * nodeIndex;
            }
        }

        final OUSetup sharedWorkspace = withPasserParallelism(base, tips, 1);
        final OUSetup freshWorkspace = withPasserParallelism(base, tips, 3);
        refreshMessages(sharedWorkspace);
        refreshMessages(freshWorkspace);

        final double[] sharedA = new double[dimTrait * dimTrait];
        final double[] sharedQ = new double[dimTrait * dimTrait];
        final double[] sharedMu = new double[dimTrait];
        final double[] freshA = new double[dimTrait * dimTrait];
        final double[] freshQ = new double[dimTrait * dimTrait];
        final double[] freshMu = new double[dimTrait];

        sharedWorkspace.passer.computeJointGradients(sharedWorkspace.provider, sharedA, sharedQ, sharedMu);
        freshWorkspace.passer.computeJointGradients(freshWorkspace.provider, freshA, freshQ, freshMu);

        assertVectorEquals("shared vs fresh A with distinct branch lengths", freshA, sharedA, 1.0e-11);
        assertVectorEquals("shared vs fresh Q with distinct branch lengths", freshQ, sharedQ, 1.0e-11);
        assertVectorEquals("shared vs fresh mu with distinct branch lengths", freshMu, sharedMu, 1.0e-11);
    }

    public void testBranchGradientCopyRegressionBenchmarkReusesPreparedInputs() {
        System.out.println("\nBenchmark: branch gradient prepared-input reuse");

        withTransitionCacheDiagnostics(new CheckedRunnable() {
            @Override
            public void run() {
                final OrthogonalBlockOUSetup setup =
                        buildOrthogonalBlockOUSetup("focusedBranchGradientCopyBench", buildPartiallyObservedTips());
                final OUSetup parallel =
                        withPasserParallelism(setup.orthogonal, buildPartiallyObservedTips(), 3);
                refreshMessages(parallel);

                final BranchGradientInputs inputs =
                        parallel.passer.prepareBranchGradientInputs(parallel.provider);
                final long prepRequests =
                        parallel.provider.getTransitionCacheRequestCount(CanonicalTransitionCachePhases.GRADIENT_PREP);
                final long prepMisses =
                        parallel.provider.getTransitionCacheMissCount(CanonicalTransitionCachePhases.GRADIENT_PREP);
                final long transitionRebuilds = parallel.provider.getTransitionCacheRebuildCount();
                final long preparedBranchRebuilds = parallel.provider.getPreparedBranchRebuildCount();

                final double[] gradA = new double[setup.nativeParameter.getDimension()];
                final double[] gradQ = new double[dimTrait * dimTrait];
                final double[] gradMu = new double[dimTrait];

                final long startNanos = System.nanoTime();
                for (int repeat = 0; repeat < BRANCH_GRADIENT_BENCHMARK_REPEATS; ++repeat) {
                    parallel.passer.computeJointGradients(parallel.provider, inputs, gradA, gradQ, gradMu);
                }
                final long elapsedNanos = System.nanoTime() - startNanos;

                System.out.printf("  Prepared-input branch gradients: repeats=%d activeBranches=%d nanosPerRun=%d%n",
                        BRANCH_GRADIENT_BENCHMARK_REPEATS,
                        inputs.getActiveBranchCount(),
                        elapsedNanos / BRANCH_GRADIENT_BENCHMARK_REPEATS);
                System.out.printf("  Final gradients A=%s Q=%s mu=%s%n",
                        Arrays.toString(gradA), Arrays.toString(gradQ), Arrays.toString(gradMu));

                assertEquals("prepared-input benchmark does not redo gradient prep requests",
                        prepRequests,
                        parallel.provider.getTransitionCacheRequestCount(CanonicalTransitionCachePhases.GRADIENT_PREP));
                assertEquals("prepared-input benchmark has no gradient prep misses",
                        prepMisses,
                        parallel.provider.getTransitionCacheMissCount(CanonicalTransitionCachePhases.GRADIENT_PREP));
                assertEquals("prepared-input benchmark does not rebuild transitions",
                        transitionRebuilds,
                        parallel.provider.getTransitionCacheRebuildCount());
                assertEquals("prepared-input benchmark does not rebuild prepared branch bases",
                        preparedBranchRebuilds,
                        parallel.provider.getPreparedBranchRebuildCount());
            }
        });
    }

    private void checkMeanGradientMatchesDense(final OrthogonalBlockOUSetup setup,
                                               final String label) {
        refreshMessages(setup.orthogonal);
        refreshMessages(setup.dense);

        final double[] orthogonalGradient = new double[dimTrait];
        final double[] denseGradient = new double[dimTrait];
        computeJointGradientMu(setup.orthogonal, orthogonalGradient, setup.nativeParameter.getDimension());
        computeJointGradientMu(setup.dense, denseGradient, dimTrait * dimTrait);

        System.out.printf("  Orthogonal dlogL/dmu: %s%n", Arrays.toString(orthogonalGradient));
        System.out.printf("  Dense      dlogL/dmu: %s%n", Arrays.toString(denseGradient));
        assertVectorEquals("orthogonal-block vs dense dlogL/dmu (" + label + ")",
                denseGradient, orthogonalGradient, TOL);
    }

    private void checkNativeSelectionGradientMatchesDenseFiniteDifference(final OrthogonalBlockOUSetup setup,
                                                                          final String label) {
        refreshMessages(setup.orthogonal);
        final double[] orthogonalGradient = new double[setup.nativeParameter.getDimension()];
        computeJointGradientA(setup.orthogonal, orthogonalGradient);

        final double[] denseNativeGradient = new double[setup.nativeParameter.getDimension()];
        for (int i = 0; i < setup.nativeParameter.getDimension(); ++i) {
            final double original = setup.nativeParameter.getParameterValue(i);

            setup.nativeParameter.setParameterValue(i, original + FD_DELTA);
            copyBlockSelectionToDenseSelection(setup);
            final double logLPlus =
                    setup.dense.passer.computePostOrderLogLikelihood(setup.dense.provider, setup.dense.rootPrior);

            setup.nativeParameter.setParameterValue(i, original - FD_DELTA);
            copyBlockSelectionToDenseSelection(setup);
            final double logLMinus =
                    setup.dense.passer.computePostOrderLogLikelihood(setup.dense.provider, setup.dense.rootPrior);

            setup.nativeParameter.setParameterValue(i, original);
            copyBlockSelectionToDenseSelection(setup);

            denseNativeGradient[i] = (logLPlus - logLMinus) / (2.0 * FD_DELTA);
        }

        System.out.printf("  Orthogonal native dlogL/dA: %s%n", Arrays.toString(orthogonalGradient));
        System.out.printf("  Dense native numeric dlogL/dA: %s%n", Arrays.toString(denseNativeGradient));
        assertVectorEquals("orthogonal-block native selection gradient vs dense finite difference (" + label + ")",
                denseNativeGradient, orthogonalGradient, TOL);
    }

    private void copyBlockSelectionToDenseSelection(final OrthogonalBlockOUSetup setup) {
        fillMatrixParameter(setup.dense.aMatrix, setup.blockSelection.getParameterAsMatrix());
    }

    private CanonicalTipObservation observedTip(final double x0, final double x1, final double x2) {
        final CanonicalTipObservation observation = new CanonicalTipObservation(dimTrait);
        observation.setObserved(new double[]{x0, x1, x2});
        return observation;
    }

    private long clearCount(final OUSetup setup, final String reason) {
        return setup.provider.getTransitionCacheClearCount(reason);
    }

    private int firstNonRootNodeIndex() {
        final int rootIndex = treeModel.getRoot().getNumber();
        for (int nodeIndex = 0; nodeIndex < treeModel.getNodeCount(); ++nodeIndex) {
            if (nodeIndex != rootIndex) {
                return nodeIndex;
            }
        }
        throw new IllegalStateException("test tree has no non-root nodes");
    }

    private void assertVectorEquals(final String label,
                                    final double[] expected,
                                    final double[] actual,
                                    final double tolerance) {
        assertEquals(label + " length", expected.length, actual.length);
        for (int i = 0; i < expected.length; ++i) {
            assertEquals(label + "[" + i + "]", expected[i], actual[i], tolerance);
        }
    }

    private void withTransitionCacheDiagnostics(final CheckedRunnable body) {
        final String previous = System.getProperty(TRANSITION_CACHE_DIAGNOSTICS_PROPERTY);
        System.setProperty(TRANSITION_CACHE_DIAGNOSTICS_PROPERTY, "true");
        try {
            body.run();
        } finally {
            if (previous == null) {
                System.clearProperty(TRANSITION_CACHE_DIAGNOSTICS_PROPERTY);
            } else {
                System.setProperty(TRANSITION_CACHE_DIAGNOSTICS_PROPERTY, previous);
            }
        }
    }

    private interface CheckedRunnable {
        void run();
    }
}
