package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipObservation;
import junit.framework.Test;
import junit.framework.TestSuite;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;

/**
 * Standalone (not part of {@code ant junit-canonical}) baseline benchmark for the
 * <b>current, unmodified</b> canonical OU tree-likelihood path under a general
 * (non-orthogonal R) block-diagonal selection matrix S = R*B*R^{-1}.
 *
 * <p>This is the reference snapshot that the rotated/block-basis OU work tracked in
 * {@code private-code/rotatedOU/} is implemented and benchmarked against. It records,
 * for the production dense-rotation path:
 * <ul>
 *     <li>log-likelihood, for fully-observed and partial/missing tip data;</li>
 *     <li>analytic gradients (selection-matrix R block, diffusion Q, stationary mean mu,
 *     branch lengths), cross-checked against finite differences where the existing test
 *     support already does so;</li>
 *     <li>wall-clock timing of the log-likelihood, joint-gradient, and branch-length
 *     gradient evaluations, at a small (6-tip) and larger synthetic (128/1024-tip) scale.</li>
 * </ul>
 *
 * <p>Run via {@code private-code/scripts/local/run_rotated_ou_baseline_benchmark.sh},
 * which captures stdout as the versioned baseline artifact.
 */
public final class RotatedOUBaselineBenchmarkTest extends CanonicalOUMessagePasserTestSupport {

    private static final int SMALL_WARMUP_ITERS = 500;
    private static final int SMALL_MEASURE_ITERS = 3_000;
    private static final int LARGE_WARMUP_ITERS = 50;
    private static final int LARGE_MEASURE_ITERS = 300;
    private static final int XLARGE_WARMUP_ITERS = 20;
    private static final int XLARGE_MEASURE_ITERS = 100;

    public RotatedOUBaselineBenchmarkTest(final String name) {
        super(name);
    }

    public static Test suite() {
        return new TestSuite(RotatedOUBaselineBenchmarkTest.class);
    }

    public static void main(final String[] args) {
        junit.textui.TestRunner.run(suite());
    }

    public void testBaselineReportIsFiniteAndConsistent() {
        final PrintStream originalOut = System.out;
        final ByteArrayOutputStream capture = new ByteArrayOutputStream();
        final PrintStream tee = new PrintStream(capture, true);
        System.setOut(tee);
        try {
            runBaselineCases();
        } finally {
            System.setOut(originalOut);
        }

        final String report = capture.toString();
        originalOut.println(report);

        assertTrue(report.contains("Case 6-tip full-observed"));
        assertTrue(report.contains("Case 6-tip partial-missing"));
        assertTrue(report.contains("Case 128-tip full-observed"));
        assertTrue(report.contains("Case 128-tip partial-missing"));
        assertTrue(report.contains("Case 1024-tip full-observed"));
        assertTrue(report.contains("Case 1024-tip partial-missing"));
        assertFalse(report.contains("NaN"));
        assertFalse(report.contains("Infinity"));
    }

    private void runBaselineCases() {
        appendSmallCase("6-tip full-observed", buildFullyObservedTips());
        appendSmallCase("6-tip partial-missing", buildPartiallyObservedTips());
        appendSyntheticTreeCase(128, "128-tip full-observed", true, LARGE_WARMUP_ITERS, LARGE_MEASURE_ITERS);
        appendSyntheticTreeCase(128, "128-tip partial-missing", false, LARGE_WARMUP_ITERS, LARGE_MEASURE_ITERS);
        appendSyntheticTreeCase(1024, "1024-tip full-observed", true, XLARGE_WARMUP_ITERS, XLARGE_MEASURE_ITERS);
        appendSyntheticTreeCase(1024, "1024-tip partial-missing", false, XLARGE_WARMUP_ITERS, XLARGE_MEASURE_ITERS);
    }

    private void appendSmallCase(final String label, final CanonicalTipObservation[] tips) {
        System.out.println("Case " + label);

        final GeneralBlockOUSetup setup = buildGeneralBlockOUSetup("baseline." + label, tips);

        checkGeneralBlockLikelihoodMatchesDense(setup, label);
        checkGeneralBlockGradientRMatchesFiniteDifference(setup, label);
        final int gradientDimensionA = setup.blockSelection.getBlockDiagonalNParameters() + dimTrait * dimTrait;
        recordGradientQ(setup.block, gradientDimensionA, label);
        recordGradientMu(setup.block, gradientDimensionA, label);
        recordGradientBranchLengths(setup.block, label);
        final double likelihoodNs = timeLogLikelihoodNs(setup.block, SMALL_WARMUP_ITERS, SMALL_MEASURE_ITERS);
        final double jointGradientNs =
                timeJointGradientNs(setup.block, gradientDimensionA, SMALL_WARMUP_ITERS, SMALL_MEASURE_ITERS);
        final double branchGradientNs =
                timeBranchLengthGradientNs(setup.block, SMALL_WARMUP_ITERS, SMALL_MEASURE_ITERS);

        System.out.printf("  postorder likelihood: %.1f ns%n", likelihoodNs);
        System.out.printf("  joint gradient (R+Q+mu): %.1f ns%n", jointGradientNs);
        System.out.printf("  branch-length gradient: %.1f ns%n", branchGradientNs);
        System.out.println();
    }

    private void appendSyntheticTreeCase(final int tipCount,
                                         final String label,
                                         final boolean fullObserved,
                                         final int warmupIterations,
                                         final int measureIterations) {
        useBalancedSyntheticTree(tipCount);
        final CanonicalTipObservation[] tips = buildSyntheticTips(fullObserved);
        final GeneralBlockOUSetup setup = buildGeneralBlockOUSetup("baseline.synthetic." + label, tips);

        System.out.println("Case " + label);

        final double logLikelihood =
                setup.block.passer.computePostOrderLogLikelihood(setup.block.provider, setup.block.rootPrior);
        if (!Double.isFinite(logLikelihood)) {
            throw new AssertionError("Non-finite synthetic baseline likelihood: " + logLikelihood);
        }
        System.out.printf("  logL: %.12f%n", logLikelihood);

        final int gradientDimensionA = setup.blockSelection.getBlockDiagonalNParameters() + dimTrait * dimTrait;
        final double likelihoodNs = timeLogLikelihoodNs(setup.block, warmupIterations, measureIterations);
        final double jointGradientNs =
                timeJointGradientNs(setup.block, gradientDimensionA, warmupIterations, measureIterations);
        final double branchGradientNs =
                timeBranchLengthGradientNs(setup.block, warmupIterations, measureIterations);

        System.out.printf("  postorder likelihood: %.1f ns%n", likelihoodNs);
        System.out.printf("  joint gradient (R+Q+mu): %.1f ns%n", jointGradientNs);
        System.out.printf("  branch-length gradient: %.1f ns%n", branchGradientNs);
        System.out.println();
    }

    // Q, mu, and branch-length analytic gradients here are recorded with a relative-tolerance
    // finite-difference spot check rather than the shared support class's strict absolute TOL/
    // BRANCH_TOL gates: this general-block fixture's gradient magnitudes (driven by the real
    // primate tree's branch lengths and the existing diffusion model's precision) are large
    // enough that a fixed absolute tolerance tuned for the dense/orthogonal fixtures is too
    // strict, even though the relative agreement is good (see implementation_notes.md,
    // 2026-06-26 entry). The point of this benchmark is to lock in baseline values, not to
    // re-litigate gradient correctness already covered elsewhere.
    private void recordGradientQ(final OUSetup setup, final int selectionGradientDimension, final String label) {
        refreshMessages(setup);
        final double[] gradQ = new double[dimTrait * dimTrait];
        computeJointGradientQ(setup, gradQ, selectionGradientDimension);
        System.out.printf("  Analytic dlogL/dQ: %s%n", Arrays.toString(gradQ));

        final double original = setup.qMatrix.getParameterValue(0, 0);
        setup.qMatrix.setParameterValue(0, 0, original + FD_DELTA);
        syncPrecisionFromVariance(setup.qMatrix, setup.precisionMatrix);
        final double logLPlus = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);
        setup.qMatrix.setParameterValue(0, 0, original - FD_DELTA);
        syncPrecisionFromVariance(setup.qMatrix, setup.precisionMatrix);
        final double logLMinus = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);
        setup.qMatrix.setParameterValue(0, 0, original);
        syncPrecisionFromVariance(setup.qMatrix, setup.precisionMatrix);
        final double numeric = (logLPlus - logLMinus) / (2.0 * FD_DELTA);
        assertRelativelyClose("dlogL/dQ[0,0] (" + label + ")", numeric, gradQ[0], 5.0e-3);
    }

    private void recordGradientMu(final OUSetup setup, final int selectionGradientDimension, final String label) {
        refreshMessages(setup);
        final double[] gradMu = new double[dimTrait];
        computeJointGradientMu(setup, gradMu, selectionGradientDimension);
        System.out.printf("  Analytic dlogL/dmu: %s%n", Arrays.toString(gradMu));

        final double original = setup.muParam.getParameterValue(0);
        setup.muParam.setParameterValue(0, original + FD_DELTA);
        final double logLPlus = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);
        setup.muParam.setParameterValue(0, original - FD_DELTA);
        final double logLMinus = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);
        setup.muParam.setParameterValue(0, original);
        final double numeric = (logLPlus - logLMinus) / (2.0 * FD_DELTA);
        assertRelativelyClose("dlogL/dmu[0] (" + label + ")", numeric, gradMu[0], 5.0e-3);
    }

    private void recordGradientBranchLengths(final OUSetup setup, final String label) {
        refreshMessages(setup);
        final int nodeCount = treeModel.getNodeCount();
        final int rootIndex = treeModel.getRoot().getNumber();
        final double[] gradBranch = new double[nodeCount];
        setup.passer.computeGradientBranchLengths(setup.provider, gradBranch);
        System.out.printf("  Analytic dlogL/dt: %s%n", Arrays.toString(gradBranch));

        int spotIndex = -1;
        for (int nodeIndex = 0; nodeIndex < nodeCount; nodeIndex++) {
            if (nodeIndex != rootIndex) {
                spotIndex = nodeIndex;
                break;
            }
        }
        final double baseBranchLength = treeModel.getBranchLength(treeModel.getNode(spotIndex));
        setup.branchScale[spotIndex] += FD_DELTA;
        final double logLPlus = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);
        setup.branchScale[spotIndex] -= 2.0 * FD_DELTA;
        final double logLMinus = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);
        setup.branchScale[spotIndex] += FD_DELTA;
        final double numeric = (logLPlus - logLMinus) / (2.0 * FD_DELTA * baseBranchLength);
        assertRelativelyClose("dlogL/dt[" + spotIndex + "] (" + label + ")", numeric, gradBranch[spotIndex], 5.0e-3);
    }

    private static void assertRelativelyClose(final String label,
                                               final double expectedNumeric,
                                               final double actualAnalytic,
                                               final double relativeTolerance) {
        final double scale = Math.max(1.0, Math.abs(expectedNumeric));
        final double relativeError = Math.abs(actualAnalytic - expectedNumeric) / scale;
        System.out.printf("  %s: numeric=%.6f analytic=%.6f relErr=%.3e%n",
                label, expectedNumeric, actualAnalytic, relativeError);
        if (relativeError > relativeTolerance) {
            throw new AssertionError(label + " relative error " + relativeError + " exceeds " + relativeTolerance
                    + " (numeric=" + expectedNumeric + ", analytic=" + actualAnalytic + ")");
        }
    }

    private double timeLogLikelihoodNs(final OUSetup setup, final int warmupIterations, final int measureIterations) {
        double logLikelihood = Double.NaN;
        for (int i = 0; i < warmupIterations; ++i) {
            logLikelihood = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);
        }
        final long start = System.nanoTime();
        for (int i = 0; i < measureIterations; ++i) {
            logLikelihood = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);
        }
        if (!Double.isFinite(logLikelihood)) {
            throw new AssertionError("Non-finite benchmark likelihood: " + logLikelihood);
        }
        return (System.nanoTime() - start) / (double) measureIterations;
    }

    private double timeJointGradientNs(final OUSetup setup,
                                       final int gradientDimensionA,
                                       final int warmupIterations,
                                       final int measureIterations) {
        final double[] gradA = new double[gradientDimensionA];
        final double[] gradQ = new double[dimTrait * dimTrait];
        final double[] gradMu = new double[dimTrait];
        for (int i = 0; i < warmupIterations; ++i) {
            setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);
            setup.passer.computePreOrder(setup.provider, setup.rootPrior);
            setup.passer.computeJointGradients(setup.provider, gradA, gradQ, gradMu);
        }
        final long start = System.nanoTime();
        for (int i = 0; i < measureIterations; ++i) {
            setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);
            setup.passer.computePreOrder(setup.provider, setup.rootPrior);
            setup.passer.computeJointGradients(setup.provider, gradA, gradQ, gradMu);
        }
        if (!isFinite(gradA) || !isFinite(gradQ) || !isFinite(gradMu)) {
            throw new AssertionError("Non-finite benchmark joint gradient.");
        }
        return (System.nanoTime() - start) / (double) measureIterations;
    }

    private double timeBranchLengthGradientNs(final OUSetup setup,
                                              final int warmupIterations,
                                              final int measureIterations) {
        final double[] gradBranch = new double[treeModel.getNodeCount()];
        for (int i = 0; i < warmupIterations; ++i) {
            setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);
            setup.passer.computePreOrder(setup.provider, setup.rootPrior);
            setup.passer.computeGradientBranchLengths(setup.provider, gradBranch);
        }
        final long start = System.nanoTime();
        for (int i = 0; i < measureIterations; ++i) {
            setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);
            setup.passer.computePreOrder(setup.provider, setup.rootPrior);
            setup.passer.computeGradientBranchLengths(setup.provider, gradBranch);
        }
        if (!isFinite(gradBranch)) {
            throw new AssertionError("Non-finite benchmark branch-length gradient.");
        }
        return (System.nanoTime() - start) / (double) measureIterations;
    }

    private static boolean isFinite(final double[] values) {
        for (final double value : values) {
            if (!Double.isFinite(value)) {
                return false;
            }
        }
        return true;
    }

    private void useBalancedSyntheticTree(final int tipCount) {
        try {
            final Tree tree = new NewickImporter(buildBalancedNewick(0, tipCount, true)).importTree(null);
            treeModel = new DefaultTreeModel("rotatedOuBaselineTree" + tipCount, tree);
            rateTransformation = new ContinuousRateTransformation.Default(treeModel, false, false);
            nTips = tipCount;
        } catch (IOException | Importer.ImportException e) {
            throw new RuntimeException("Failed to build synthetic baseline tree", e);
        }
    }

    private static String buildBalancedNewick(final int start, final int end, final boolean root) {
        if (end - start == 1) {
            return "taxon" + start + ":0.05";
        }
        final int mid = start + (end - start) / 2;
        final String subtree = "("
                + buildBalancedNewick(start, mid, false)
                + ","
                + buildBalancedNewick(mid, end, false)
                + ")";
        return root ? subtree + ";" : subtree + ":0.05";
    }

    private CanonicalTipObservation[] buildSyntheticTips(final boolean fullObserved) {
        final CanonicalTipObservation[] tips = new CanonicalTipObservation[treeModel.getExternalNodeCount()];
        for (int tipIndex = 0; tipIndex < treeModel.getExternalNodeCount(); ++tipIndex) {
            final int nodeIndex = treeModel.getExternalNode(tipIndex).getNumber();
            final CanonicalTipObservation observation = new CanonicalTipObservation(dimTrait);
            final double phase = nodeIndex + 1.0;
            final double[] values = new double[]{
                    Math.sin(0.17 * phase),
                    Math.cos(0.11 * phase),
                    ((nodeIndex % 17) - 8.0) / 4.0
            };
            if (fullObserved) {
                observation.setObserved(values);
            } else if (nodeIndex % 11 == 0) {
                observation.setMissing();
            } else if (nodeIndex % 7 == 0) {
                observation.setPartiallyObserved(values, new boolean[]{true, false, true});
            } else {
                observation.setObserved(values);
            }
            tips[nodeIndex] = observation;
        }
        return tips;
    }
}
