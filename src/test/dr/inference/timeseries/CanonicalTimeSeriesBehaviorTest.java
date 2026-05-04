package test.dr.inference.timeseries;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Focused guard suite for the canonical time-series Gaussian pathway.
 */
public final class CanonicalTimeSeriesBehaviorTest extends TestCase {

    public CanonicalTimeSeriesBehaviorTest(final String name) {
        super(name);
    }

    public static Test suite() {
        final TestSuite suite = new TestSuite("Canonical time-series behavior");
        addAnalyticalKalmanTests(suite);
        suite.addTest(new OUProcessModelTest("testDefaultCovarianceGradientMethodUsesStationaryLyapunovForOrthogonalBlockDrift"));
        suite.addTest(new OUProcessModelTest("testUsesOrthogonalBlockSelectionChartDetectsExpectedCharts"));
        return suite;
    }

    public static void main(final String[] args) {
        junit.textui.TestRunner.run(suite());
    }

    private static void addAnalyticalKalmanTests(final TestSuite suite) {
        final String[] testNames = {
                "testCanonicalStationaryMeanGradientMatchesFD_2D",
                "testExactOuCanonicalStateAndTransitionMatchExpectationForm",
                "testEulerCanonicalTransitionMatchesExpectationForm",
                "testCanonicalKalmanLikelihoodMatchesExpectationEngine_ExactOu",
                "testCanonicalKalmanLikelihoodMatchesExpectationEngine_Euler",
                "testForwardLikelihoodFactorySelectsEquivalentExactOuBackends",
                "testForwardLikelihoodFactorySelectsEquivalentEulerBackends",
                "testTimeSeriesLikelihoodFactorySupportsCanonicalForwardMode",
                "testCanonicalKalmanSmootherMatchesExpectationSmoother_ExactOu",
                "testCanonicalKalmanSmootherMatchesExpectationSmoother_Euler",
                "testTimeSeriesLikelihoodFactorySupportsCanonicalAnalyticalGradients",
                "testParallelTimeSeriesLikelihoodSumsIndependentCanonicalOuSeries",
                "testParallelCanonicalOrthogonalBlockFullEvaluationAfterParameterMoves",
                "testParallelCanonicalOrthogonalBlockNativeGradientsMatchFD_2D",
                "testParallelCanonicalOrthogonalBlockBatchGradientsMatchSingleGradients",
                "testCanonicalAnalyticalDriftGradientMatchesFD_2D_NonDiagonal",
                "testCanonicalAnalyticalDiffusionGradientMatchesFD_2D_NonDiagonal",
                "testCanonicalAnalyticalGradientsMatchFD_3D_Dense",
                "testCanonicalAnalyticalGradientsMatchFD_2D_WeakStability",
                "testCanonicalAnalyticalOrthogonalBlockNativeGradientsMatchFD_2D",
                "testCanonicalAnalyticalOrthogonalBlockNativeGradientsMatchFD_4D",
                "testLatentMiddleThreePointChainBranchSumMatchesCanonicalGradient",
                "testTreeBranchProviderMatchesTimeSeriesKernel",
                "testOrthogonalBlockCanonicalTransitionMatchesDense_2D",
                "testOrthogonalBlockCanonicalTransitionMatchesDense_4D",
                "testOrthogonalBlockCanonicalLocalAdjointsMatchDense_2D",
                "testOrthogonalBlockCanonicalLocalAdjointsMatchDense_4D"
        };
        for (final String testName : testNames) {
            suite.addTest(new AnalyticalKalmanGradientEngineTest(testName));
        }
    }
}
