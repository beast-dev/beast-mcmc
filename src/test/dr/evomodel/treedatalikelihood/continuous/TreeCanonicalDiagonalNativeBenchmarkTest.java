package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.AbstractBranchRateModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitDataModel;
import dr.evomodel.treedatalikelihood.continuous.DiffusionProcessDelegate;
import dr.evomodel.treedatalikelihood.continuous.OUDiffusionModelDelegate;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTransitionCachePhases;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodel.treedatalikelihood.continuous.cdi.SafeMultivariateDiagonalActualizedWithDriftIntegrator;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipObservation;
import dr.inference.model.CompoundParameter;
import dr.inference.model.DiagonalMatrix;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import junit.framework.Test;
import junit.framework.TestSuite;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class TreeCanonicalDiagonalNativeBenchmarkTest extends CanonicalOUMessagePasserTestSupport {

    private static final String NATIVE_DIAGONAL_MESSAGES_PROPERTY =
            "beast.experimental.canonicalUseNativeDiagonalOUMessages";
    private static final String TRANSITION_CACHE_DIAGNOSTICS_PROPERTY =
            "beast.debug.canonicalTransitionCache";
    private static final int SMALL_WARMUP_ITERS = 1_000;
    private static final int SMALL_MEASURE_ITERS = 10_000;
    private static final int SMALL_CHANGING_BRANCH_MEASURE_ITERS = 2_000;
    private static final int LARGE_WARMUP_ITERS = 300;
    private static final int LARGE_MEASURE_ITERS = 2_000;
    private static final int LARGE_CHANGING_BRANCH_MEASURE_ITERS = 500;
    private static final int LEGACY_WARMUP_ITERS = 100;
    private static final int LEGACY_MEASURE_ITERS = 500;
    private static final int LEGACY_CHANGING_BRANCH_MEASURE_ITERS = 200;

    public TreeCanonicalDiagonalNativeBenchmarkTest(final String name) {
        super(name);
    }

    public void testBenchmarkReportIsFinite() {
        final String report = buildReport();
        System.out.println(report);
        assertTrue(report.contains("Case 6-tip full-observed"));
        assertTrue(report.contains("Case 6-tip partial-missing"));
        assertTrue(report.contains("Case 128-tip full-observed"));
        assertTrue(report.contains("Case 128-tip partial-missing"));
        assertTrue(report.contains("SafeMultivariateDiagonalActualizedWithDriftIntegrator"));
        assertFalse(report.contains("NaN"));
        assertFalse(report.contains("Infinity"));
    }

    public static Test suite() {
        return new TestSuite(TreeCanonicalDiagonalNativeBenchmarkTest.class);
    }

    public static void main(final String[] args) {
        junit.textui.TestRunner.run(suite());
    }

    private String buildReport() {
        final StringBuilder sb = new StringBuilder();
        appendCase(sb, "6-tip full-observed", true,
                buildFullyObservedTips(),
                SMALL_WARMUP_ITERS,
                SMALL_MEASURE_ITERS,
                SMALL_CHANGING_BRANCH_MEASURE_ITERS);
        appendCase(sb, "6-tip partial-missing", false,
                buildPartiallyObservedTips(),
                SMALL_WARMUP_ITERS,
                SMALL_MEASURE_ITERS,
                SMALL_CHANGING_BRANCH_MEASURE_ITERS);
        appendSyntheticTreeCase(sb, 128, "128-tip full-observed", true);
        appendSyntheticTreeCase(sb, 128, "128-tip partial-missing", false);
        return sb.toString();
    }

    private void appendCase(final StringBuilder sb,
                            final String label,
                            final boolean fullObserved,
                            final CanonicalTipObservation[] tips,
                            final int warmupIterations,
                            final int measureIterations,
                            final int changingBranchIterations) {
        final String previous = System.getProperty(NATIVE_DIAGONAL_MESSAGES_PROPERTY);
        final String previousDiagnostics = System.getProperty(TRANSITION_CACHE_DIAGNOSTICS_PROPERTY);
        try {
            System.setProperty(TRANSITION_CACHE_DIAGNOSTICS_PROPERTY, "true");
            System.setProperty(NATIVE_DIAGONAL_MESSAGES_PROPERTY, "true");
            final OUSetup nativeSetup = buildDiagonalOUSetup(
                    "benchmarkNative." + label,
                    tips);
            final double nativeNs = timePostorder(nativeSetup, warmupIterations, measureIterations);
            final double nativeChangingBranchNs = timePostorderWithChangingBranches(
                    nativeSetup,
                    warmupIterations,
                    changingBranchIterations);
            final double nativeLog =
                    nativeSetup.passer.computePostOrderLogLikelihood(nativeSetup.provider, nativeSetup.rootPrior);

            System.setProperty(NATIVE_DIAGONAL_MESSAGES_PROPERTY, "false");
            final OUSetup fallbackSetup = buildDiagonalOUSetup(
                    "benchmarkFallback." + label,
                    copyTips(tips));
            final double fallbackNs = timePostorder(fallbackSetup, warmupIterations, measureIterations);
            final double fallbackChangingBranchNs = timePostorderWithChangingBranches(
                    fallbackSetup,
                    warmupIterations,
                    changingBranchIterations);
            final double fallbackLog =
                    fallbackSetup.passer.computePostOrderLogLikelihood(
                            fallbackSetup.provider,
                            fallbackSetup.rootPrior);

            assertEquals("native vs fallback log-likelihood (" + label + ")",
                    fallbackLog,
                    nativeLog,
                    1.0e-9);

            sb.append("Case ").append(label).append('\n');
            sb.append("  steady postorder likelihood\n");
            sb.append(String.format("    native diagonal messages: %.1f ns%n", nativeNs));
            sb.append(String.format("    full-transition fallback: %.1f ns%n", fallbackNs));
            sb.append(String.format("    speedup: %.2fx%n", fallbackNs / nativeNs));
            sb.append("  changing branch times\n");
            sb.append(String.format("    native diagonal messages: %.1f ns%n", nativeChangingBranchNs));
            sb.append(String.format("    full-transition fallback: %.1f ns%n", fallbackChangingBranchNs));
            sb.append(String.format("    speedup: %.2fx%n", fallbackChangingBranchNs / nativeChangingBranchNs));
            appendCacheLine(sb, "native", nativeSetup);
            appendCacheLine(sb, "fallback", fallbackSetup);
            sb.append('\n');
        } finally {
            if (previous == null) {
                System.clearProperty(NATIVE_DIAGONAL_MESSAGES_PROPERTY);
            } else {
                System.setProperty(NATIVE_DIAGONAL_MESSAGES_PROPERTY, previous);
            }
            if (previousDiagnostics == null) {
                System.clearProperty(TRANSITION_CACHE_DIAGNOSTICS_PROPERTY);
            } else {
                System.setProperty(TRANSITION_CACHE_DIAGNOSTICS_PROPERTY, previousDiagnostics);
            }
        }
    }

    private void appendSyntheticTreeCase(final StringBuilder sb,
                                         final int tipCount,
                                         final String label,
                                         final boolean fullObserved) {
        useBalancedSyntheticTree(tipCount);
        final CanonicalTipObservation[] tips = buildSyntheticTips(fullObserved);
        appendCase(sb, label, fullObserved,
                tips,
                LARGE_WARMUP_ITERS,
                LARGE_MEASURE_ITERS,
                LARGE_CHANGING_BRANCH_MEASURE_ITERS);
        appendSafeMultivariateComparisonCase(sb, label, tips);
    }

    private static double timePostorder(final OUSetup setup,
                                        final int warmupIterations,
                                        final int measureIterations) {
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

    private static double timePostorderWithChangingBranches(final OUSetup setup,
                                                            final int warmupIterations,
                                                            final int measureIterations) {
        double logLikelihood = Double.NaN;
        for (int i = 0; i < warmupIterations; ++i) {
            perturbBranchScales(setup, i);
            logLikelihood = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);
        }
        final long start = System.nanoTime();
        for (int i = 0; i < measureIterations; ++i) {
            perturbBranchScales(setup, warmupIterations + i);
            logLikelihood = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);
        }
        if (!Double.isFinite(logLikelihood)) {
            throw new AssertionError("Non-finite changing-branch benchmark likelihood: " + logLikelihood);
        }
        return (System.nanoTime() - start) / (double) measureIterations;
    }

    private static void perturbBranchScales(final OUSetup setup, final int iteration) {
        perturbBranchScales(setup.branchScale, iteration);
    }

    private static void perturbBranchScales(final double[] branchScale, final int iteration) {
        for (int nodeIndex = 0; nodeIndex < branchScale.length; ++nodeIndex) {
            branchScale[nodeIndex] = 1.0 + 1.0e-4 * ((iteration + nodeIndex) % 17);
        }
    }

    private static void appendCacheLine(final StringBuilder sb,
                                        final String label,
                                        final OUSetup setup) {
        sb.append(String.format(
                "  %s full-transition cache: requests=%d misses=%d rebuilds=%d%n",
                label,
                setup.provider.getTransitionCacheRequestCount(CanonicalTransitionCachePhases.POSTORDER),
                setup.provider.getTransitionCacheMissCount(CanonicalTransitionCachePhases.POSTORDER),
                setup.provider.getTransitionCacheRebuildCount()));
    }

    private void useBalancedSyntheticTree(final int tipCount) {
        try {
            final Tree tree = new NewickImporter(buildBalancedNewick(0, tipCount, true)).importTree(null);
            treeModel = new DefaultTreeModel("diagonalOuBenchmarkTree" + tipCount, tree);
            rateTransformation = new ContinuousRateTransformation.Default(treeModel, false, false);
            nTips = tipCount;
        } catch (IOException | Importer.ImportException e) {
            throw new RuntimeException("Failed to build synthetic benchmark tree", e);
        }
    }

    private static String buildBalancedNewick(final int start,
                                              final int end,
                                              final boolean root) {
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
        final CanonicalTipObservation[] tips =
                new CanonicalTipObservation[treeModel.getExternalNodeCount()];
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

    private static CanonicalTipObservation[] copyTips(final CanonicalTipObservation[] tips) {
        final CanonicalTipObservation[] out = new CanonicalTipObservation[tips.length];
        for (int i = 0; i < tips.length; ++i) {
            out[i] = new CanonicalTipObservation(tips[i].dim);
            out[i].copyFrom(tips[i]);
        }
        return out;
    }

    private void appendSafeMultivariateComparisonCase(final StringBuilder sb,
                                                      final String label,
                                                      final CanonicalTipObservation[] tips) {
        final String previous = System.getProperty(NATIVE_DIAGONAL_MESSAGES_PROPERTY);
        try {
            System.setProperty(NATIVE_DIAGONAL_MESSAGES_PROPERTY, "true");
            final TreeLikelihoodSetup canonicalSetup =
                    buildTreeLikelihoodSetup("treeCanonical." + label, copyTips(tips), true);
            final TreeLikelihoodSetup legacySetup =
                    buildTreeLikelihoodSetup("treeLegacy." + label, copyTips(tips), false);

            assertTrue("legacy diagonal OU benchmark should use SafeMultivariateDiagonalActualizedWithDriftIntegrator",
                    legacySetup.delegate.getIntegrator()
                            instanceof SafeMultivariateDiagonalActualizedWithDriftIntegrator);
            assertTrue("canonical OU benchmark should enable canonical OU",
                    canonicalSetup.delegate.usesCanonicalOULikelihood());
            assertFalse("legacy diagonal OU benchmark should not enable canonical OU",
                    legacySetup.delegate.usesCanonicalOULikelihood());

            final double canonicalLog = forceLogLikelihood(canonicalSetup);
            final double legacyLog = forceLogLikelihood(legacySetup);
            assertEquals("canonical diagonal OU vs legacy diagonal CDI log-likelihood (" + label + ")",
                    legacyLog,
                    canonicalLog,
                    1.0e-7);

            final double canonicalNs = timeTreeLikelihood(
                    canonicalSetup,
                    LEGACY_WARMUP_ITERS,
                    LEGACY_MEASURE_ITERS);
            final double legacyNs = timeTreeLikelihood(
                    legacySetup,
                    LEGACY_WARMUP_ITERS,
                    LEGACY_MEASURE_ITERS);
            final double canonicalChangingBranchNs = timeTreeLikelihoodWithChangingBranches(
                    canonicalSetup,
                    LEGACY_WARMUP_ITERS,
                    LEGACY_CHANGING_BRANCH_MEASURE_ITERS);
            final double legacyChangingBranchNs = timeTreeLikelihoodWithChangingBranches(
                    legacySetup,
                    LEGACY_WARMUP_ITERS,
                    LEGACY_CHANGING_BRANCH_MEASURE_ITERS);

            sb.append("SafeMultivariate comparison ").append(label).append('\n');
            sb.append(String.format("  log likelihood: canonical=%.12f legacy=%.12f%n",
                    canonicalLog,
                    legacyLog));
            sb.append("  steady full TreeDataLikelihood\n");
            sb.append(String.format("    canonical native diagonal OU: %.1f ns%n", canonicalNs));
            sb.append(String.format("    legacy %s: %.1f ns%n",
                    legacySetup.delegate.getIntegrator().getClass().getSimpleName(),
                    legacyNs));
            sb.append(String.format("    legacy/canonical: %.2fx%n", legacyNs / canonicalNs));
            sb.append("  changing branch times\n");
            sb.append(String.format("    canonical native diagonal OU: %.1f ns%n", canonicalChangingBranchNs));
            sb.append(String.format("    legacy %s: %.1f ns%n",
                    legacySetup.delegate.getIntegrator().getClass().getSimpleName(),
                    legacyChangingBranchNs));
            sb.append(String.format("    legacy/canonical: %.2fx%n%n",
                    legacyChangingBranchNs / canonicalChangingBranchNs));
        } finally {
            if (previous == null) {
                System.clearProperty(NATIVE_DIAGONAL_MESSAGES_PROPERTY);
            } else {
                System.setProperty(NATIVE_DIAGONAL_MESSAGES_PROPERTY, previous);
            }
        }
    }

    private TreeLikelihoodSetup buildTreeLikelihoodSetup(final String tag,
                                                         final CanonicalTipObservation[] tips,
                                                         final boolean canonical) {
        final double[] branchScale = new double[treeModel.getNodeCount()];
        Arrays.fill(branchScale, 1.0);
        final BranchRateModel branchRateModel =
                new BenchmarkBranchRateModel("branchRate." + tag, branchScale);
        final ContinuousTraitDataModel benchmarkDataModel =
                buildContinuousTraitDataModel("data." + tag, tips);
        final DiffusionProcessDelegate diffusionDelegate =
                buildDiagonalOUDiffusionDelegate(tag);
        final ContinuousDataLikelihoodDelegate delegate =
                new ContinuousDataLikelihoodDelegate(
                        treeModel,
                        diffusionDelegate,
                        benchmarkDataModel,
                        buildRootPrior(tag),
                        rateTransformation,
                        branchRateModel,
                        false,
                        false);
        if (canonical) {
            delegate.enableCanonicalOULikelihood();
        }
        final TreeDataLikelihood likelihood =
                new TreeDataLikelihood(delegate, treeModel, branchRateModel);
        return new TreeLikelihoodSetup(delegate, likelihood, branchScale);
    }

    private ContinuousTraitDataModel buildContinuousTraitDataModel(final String tag,
                                                                  final CanonicalTipObservation[] tips) {
        final Parameter[] dataTraits = new Parameter[treeModel.getExternalNodeCount()];
        final boolean[] missingIndicators = new boolean[treeModel.getExternalNodeCount() * dimTrait];
        for (int tipIndex = 0; tipIndex < treeModel.getExternalNodeCount(); ++tipIndex) {
            final NodeRef node = treeModel.getExternalNode(tipIndex);
            final int nodeIndex = node.getNumber();
            final CanonicalTipObservation tip = tips[nodeIndex];
            if (tip == null) {
                throw new IllegalArgumentException("Missing benchmark tip observation for node " + nodeIndex);
            }
            dataTraits[nodeIndex] = new Parameter.Default(
                    treeModel.getNodeTaxon(node).getId(),
                    tip.values.clone());
            final int offset = nodeIndex * dimTrait;
            for (int i = 0; i < dimTrait; ++i) {
                missingIndicators[offset + i] = !tip.observed[i];
            }
        }
        return new ContinuousTraitDataModel(
                "dataModel." + tag,
                new CompoundParameter("trait." + tag, dataTraits),
                missingIndicators,
                true,
                dimTrait,
                PrecisionType.FULL);
    }

    private DiffusionProcessDelegate buildDiagonalOUDiffusionDelegate(final String tag) {
        final Parameter attenuation =
                new Parameter.Default("A.diagonal.values." + tag, new double[]{2.0, 3.0, 1.5});
        final DiagonalMatrix diagonalSelection = new DiagonalMatrix(attenuation);
        return new OUDiffusionModelDelegate(
                treeModel,
                buildDiffusionModel(tag),
                buildConstantOptima(tag),
                new MultivariateElasticModel(diagonalSelection));
    }

    private MultivariateDiffusionModel buildDiffusionModel(final String tag) {
        final MatrixParameter precisionMatrix =
                new MatrixParameter("precision." + tag, dimTrait, dimTrait);
        fillMatrixParameter(precisionMatrix, diffusionModel.getPrecisionMatrix());
        return new MultivariateDiffusionModel(precisionMatrix);
    }

    private List<BranchRateModel> buildConstantOptima(final String tag) {
        final List<BranchRateModel> optima = new ArrayList<BranchRateModel>(dimTrait);
        for (int i = 0; i < dimTrait; ++i) {
            optima.add(new StrictClockBranchRates(
                    new Parameter.Default("optimum." + tag + "." + i, new double[]{MU_VALUES[i]})));
        }
        return optima;
    }

    private static ConjugateRootTraitPrior buildRootPrior(final String tag) {
        return new ConjugateRootTraitPrior(
                new Parameter.Default("rootMean." + tag, new double[]{-1.0, -3.0, 2.5}),
                new Parameter.Default("rootSampleSize." + tag, 10.0));
    }

    private static double forceLogLikelihood(final TreeLikelihoodSetup setup) {
        setup.likelihood.makeDirty();
        final double logLikelihood = setup.likelihood.getLogLikelihood();
        if (!Double.isFinite(logLikelihood)) {
            throw new AssertionError("Non-finite full tree benchmark likelihood: " + logLikelihood);
        }
        return logLikelihood;
    }

    private static double timeTreeLikelihood(final TreeLikelihoodSetup setup,
                                             final int warmupIterations,
                                             final int measureIterations) {
        double logLikelihood = Double.NaN;
        for (int i = 0; i < warmupIterations; ++i) {
            logLikelihood = forceLogLikelihood(setup);
        }
        final long start = System.nanoTime();
        for (int i = 0; i < measureIterations; ++i) {
            logLikelihood = forceLogLikelihood(setup);
        }
        if (!Double.isFinite(logLikelihood)) {
            throw new AssertionError("Non-finite full tree benchmark likelihood: " + logLikelihood);
        }
        return (System.nanoTime() - start) / (double) measureIterations;
    }

    private static double timeTreeLikelihoodWithChangingBranches(final TreeLikelihoodSetup setup,
                                                                 final int warmupIterations,
                                                                 final int measureIterations) {
        double logLikelihood = Double.NaN;
        for (int i = 0; i < warmupIterations; ++i) {
            perturbBranchScales(setup.branchScale, i);
            logLikelihood = forceLogLikelihood(setup);
        }
        final long start = System.nanoTime();
        for (int i = 0; i < measureIterations; ++i) {
            perturbBranchScales(setup.branchScale, warmupIterations + i);
            logLikelihood = forceLogLikelihood(setup);
        }
        if (!Double.isFinite(logLikelihood)) {
            throw new AssertionError("Non-finite changing-branch full tree benchmark likelihood: "
                    + logLikelihood);
        }
        return (System.nanoTime() - start) / (double) measureIterations;
    }

    private static final class TreeLikelihoodSetup {
        final ContinuousDataLikelihoodDelegate delegate;
        final TreeDataLikelihood likelihood;
        final double[] branchScale;

        private TreeLikelihoodSetup(final ContinuousDataLikelihoodDelegate delegate,
                                    final TreeDataLikelihood likelihood,
                                    final double[] branchScale) {
            this.delegate = delegate;
            this.likelihood = likelihood;
            this.branchScale = branchScale;
        }
    }

    private static final class BenchmarkBranchRateModel extends AbstractBranchRateModel {

        private final double[] branchScale;

        private BenchmarkBranchRateModel(final String name, final double[] branchScale) {
            super(name);
            this.branchScale = branchScale;
        }

        @Override
        public double getBranchRate(final dr.evolution.tree.Tree tree, final NodeRef node) {
            return branchScale[node.getNumber()];
        }

        @Override
        protected void handleModelChangedEvent(final Model model, final Object object, final int index) { }

        @Override
        protected void handleVariableChangedEvent(final Variable variable,
                                                  final int index,
                                                  final Parameter.ChangeType type) { }

        @Override
        protected void storeState() { }

        @Override
        protected void restoreState() { }

        @Override
        protected void acceptState() { }
    }
}
