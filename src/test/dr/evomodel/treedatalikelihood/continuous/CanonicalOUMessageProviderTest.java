/*
 * CanonicalOUMessageProviderTest.java
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

import dr.evolution.tree.NodeRef;
import dr.evomodel.treedatalikelihood.continuous.canonical.adapter.CanonicalOUTreeLikelihoodIntegrator;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.message.ContinuousCanonicalPostOrderMessageProvider;
import dr.evomodel.treedatalikelihood.continuous.message.ContinuousCanonicalPreOrderMessageProvider;

public final class CanonicalOUMessageProviderTest extends CanonicalOUMessagePasserTestSupport {

    private static final String NATIVE_DIAGONAL_MESSAGES_PROPERTY =
            "beast.experimental.canonicalUseNativeDiagonalOUMessages";

    public CanonicalOUMessageProviderTest(final String name) {
        super(name);
    }

    public void testStandaloneIntegratorExposesCanonicalMessageProviders() {
        System.out.println("\nTest: canonical OU integrator exposes pre/post message providers");

        final OUSetup setup = buildOUSetup("messageProviderDense", buildPartiallyObservedTips());
        final CanonicalOUTreeLikelihoodIntegrator integrator =
                new CanonicalOUTreeLikelihoodIntegrator(setup.passer, setup.provider, setup.rootPrior);

        assertTrue("integrator is a post-order provider",
                integrator instanceof ContinuousCanonicalPostOrderMessageProvider);
        assertTrue("integrator is a pre-order provider",
                integrator instanceof ContinuousCanonicalPreOrderMessageProvider);

        final int childNodeIndex = firstNonRootInternalNodeIndex();
        final CanonicalGaussianState providerState = new CanonicalGaussianState(dimTrait);
        final CanonicalGaussianState directState = new CanonicalGaussianState(dimTrait);

        integrator.fillPostOrderNode(childNodeIndex, providerState);
        setup.passer.fillPostOrderNodeState(childNodeIndex, directState);
        assertStateEquals("post-order node provider", directState, providerState, 1.0e-12);

        integrator.fillPostOrderBranchBottom(childNodeIndex, providerState);
        setup.passer.fillPostOrderBranchBottomState(childNodeIndex, directState);
        assertStateEquals("post-order branch-bottom provider", directState, providerState, 1.0e-12);

        integrator.fillPostOrderBranchTop(childNodeIndex, providerState);
        setup.passer.fillPostOrderBranchTopState(childNodeIndex, setup.provider, directState);
        assertStateEquals("post-order branch-top provider", directState, providerState, 1.0e-12);

        integrator.fillPreOrderBranchTop(childNodeIndex, providerState);
        setup.passer.fillPreOrderBranchTopState(childNodeIndex, directState);
        assertStateEquals("pre-order branch-top provider", directState, providerState, 1.0e-12);

        integrator.fillPreOrderBranchBottom(childNodeIndex, providerState);
        setup.passer.fillPreOrderBranchBottomState(childNodeIndex, setup.provider, directState);
        assertStateEquals("pre-order branch-bottom provider", directState, providerState, 1.0e-12);

        integrator.fillRootPreOrder(providerState);
        setup.passer.fillRootPreOrderState(directState);
        assertStateEquals("root pre-order provider", directState, providerState, 1.0e-12);
    }

    public void testDiagonalNativeProviderMessagesMatchFullTransitionFallback() {
        System.out.println("\nTest: native diagonal provider messages match full-transition fallback");

        final String previous = System.getProperty(NATIVE_DIAGONAL_MESSAGES_PROPERTY);
        try {
            System.setProperty(NATIVE_DIAGONAL_MESSAGES_PROPERTY, "true");
            final OUSetup nativeSetup =
                    buildDiagonalOUSetup("messageProviderDiagonalNative", buildFullyObservedTips());
            final CanonicalOUTreeLikelihoodIntegrator nativeIntegrator =
                    new CanonicalOUTreeLikelihoodIntegrator(
                            nativeSetup.passer,
                            nativeSetup.provider,
                            nativeSetup.rootPrior);
            assertTrue("native diagonal messages enabled",
                    nativeSetup.provider.usesNativeBranchMessages());

            final StateSnapshot[] nativePostTop = allocateSnapshots(treeModel.getNodeCount());
            final StateSnapshot[] nativePreTop = allocateSnapshots(treeModel.getNodeCount());
            final StateSnapshot[] nativePreBottom = allocateSnapshots(treeModel.getNodeCount());
            snapshotProviderMessages(nativeIntegrator, nativePostTop, nativePreTop, nativePreBottom);

            System.setProperty(NATIVE_DIAGONAL_MESSAGES_PROPERTY, "false");
            final OUSetup fallbackSetup =
                    buildDiagonalOUSetup("messageProviderDiagonalFallback", buildFullyObservedTips());
            final CanonicalOUTreeLikelihoodIntegrator fallbackIntegrator =
                    new CanonicalOUTreeLikelihoodIntegrator(
                            fallbackSetup.passer,
                            fallbackSetup.provider,
                            fallbackSetup.rootPrior);
            assertFalse("full-transition fallback enabled",
                    fallbackSetup.provider.usesNativeBranchMessages());

            final CanonicalGaussianState state = new CanonicalGaussianState(dimTrait);
            for (int nodeIndex = 0; nodeIndex < treeModel.getNodeCount(); ++nodeIndex) {
                if (nodeIndex == treeModel.getRoot().getNumber()) {
                    continue;
                }

                fallbackIntegrator.fillPostOrderBranchTop(nodeIndex, state);
                nativePostTop[nodeIndex].assertMatches(
                        "diagonal post-order branch top[" + nodeIndex + "]", state, 1.0e-9);

                fallbackIntegrator.fillPreOrderBranchTop(nodeIndex, state);
                nativePreTop[nodeIndex].assertMatches(
                        "diagonal pre-order branch top[" + nodeIndex + "]", state, 1.0e-9);

                fallbackIntegrator.fillPreOrderBranchBottom(nodeIndex, state);
                nativePreBottom[nodeIndex].assertMatches(
                        "diagonal pre-order branch bottom[" + nodeIndex + "]", state, 1.0e-9);
            }
        } finally {
            if (previous == null) {
                System.clearProperty(NATIVE_DIAGONAL_MESSAGES_PROPERTY);
            } else {
                System.setProperty(NATIVE_DIAGONAL_MESSAGES_PROPERTY, previous);
            }
        }
    }

    public void testExactIdentityTipBottomRejectsFiniteCanonicalStateRequest() {
        System.out.println("\nTest: exact identity tip bottom is not exposed as finite canonical state");

        final OUSetup setup = buildOUSetup("messageProviderExactTip", buildFullyObservedTips());
        final CanonicalOUTreeLikelihoodIntegrator integrator =
                new CanonicalOUTreeLikelihoodIntegrator(setup.passer, setup.provider, setup.rootPrior);
        final int tipNodeIndex = treeModel.getExternalNode(0).getNumber();

        try {
            integrator.fillPostOrderBranchBottom(tipNodeIndex, new CanonicalGaussianState(dimTrait));
            fail("exact identity tip branch-bottom message should be degenerate");
        } catch (UnsupportedOperationException expected) {
            assertTrue("exception explains finite-state limitation",
                    expected.getMessage().contains("finite canonical child state"));
        }
    }

    private void snapshotProviderMessages(final CanonicalOUTreeLikelihoodIntegrator integrator,
                                          final StateSnapshot[] postTop,
                                          final StateSnapshot[] preTop,
                                          final StateSnapshot[] preBottom) {
        final CanonicalGaussianState state = new CanonicalGaussianState(dimTrait);
        for (int nodeIndex = 0; nodeIndex < treeModel.getNodeCount(); ++nodeIndex) {
            if (nodeIndex == treeModel.getRoot().getNumber()) {
                continue;
            }

            integrator.fillPostOrderBranchTop(nodeIndex, state);
            postTop[nodeIndex] = new StateSnapshot(state);

            integrator.fillPreOrderBranchTop(nodeIndex, state);
            preTop[nodeIndex] = new StateSnapshot(state);

            integrator.fillPreOrderBranchBottom(nodeIndex, state);
            preBottom[nodeIndex] = new StateSnapshot(state);
        }
    }

    private int firstNonRootInternalNodeIndex() {
        final int rootIndex = treeModel.getRoot().getNumber();
        for (int nodeIndex = 0; nodeIndex < treeModel.getNodeCount(); ++nodeIndex) {
            if (nodeIndex == rootIndex) {
                continue;
            }
            final NodeRef node = treeModel.getNode(nodeIndex);
            if (!treeModel.isExternal(node)) {
                return nodeIndex;
            }
        }
        throw new IllegalStateException("Tree has no non-root internal node.");
    }

    private StateSnapshot[] allocateSnapshots(final int count) {
        return new StateSnapshot[count];
    }

    private void assertStateEquals(final String label,
                                   final CanonicalGaussianState expected,
                                   final CanonicalGaussianState actual,
                                   final double tolerance) {
        assertEquals(label + " log normalizer",
                expected.logNormalizer, actual.logNormalizer, tolerance);
        assertArrayEquals(label + " information",
                expected.information, actual.information, tolerance);
        assertArrayEquals(label + " precision",
                expected.precision, actual.precision, tolerance);
    }

    private void assertArrayEquals(final String label,
                                   final double[] expected,
                                   final double[] actual,
                                   final double tolerance) {
        assertEquals(label + " length", expected.length, actual.length);
        for (int i = 0; i < expected.length; ++i) {
            assertEquals(label + "[" + i + "]", expected[i], actual[i], tolerance);
        }
    }

    private final class StateSnapshot {
        private final double[] precision;
        private final double[] information;
        private final double logNormalizer;

        private StateSnapshot(final CanonicalGaussianState state) {
            this.precision = state.precision.clone();
            this.information = state.information.clone();
            this.logNormalizer = state.logNormalizer;
        }

        private void assertMatches(final String label,
                                   final CanonicalGaussianState actual,
                                   final double tolerance) {
            assertEquals(label + " log normalizer",
                    logNormalizer, actual.logNormalizer, tolerance);
            assertArrayEquals(label + " information",
                    information, actual.information, tolerance);
            assertArrayEquals(label + " precision",
                    precision, actual.precision, tolerance);
        }
    }
}
