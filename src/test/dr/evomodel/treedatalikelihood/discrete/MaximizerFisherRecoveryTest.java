package test.dr.evomodel.treedatalikelihood.discrete;

import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.ScoreContributionsWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

public final class MaximizerFisherRecoveryTest extends TestCase {

    public void testFailureTriggersAtMostConfiguredRefreshes() {
        final MaximizerWrtParameter maximizer = runQuadraticRecovery(
                MaximizerWrtParameter.FisherInformationSource.EMPIRICAL_SCORES);

        assertEquals(1, maximizer.getFisherRefreshCount());
        assertEquals(2, maximizer.getFisherRank());
        assertEquals(0, maximizer.getFisherNullity());
        assertTrue(maximizer.fisherUsedEmpiricalScores());
        assertEquals(MaximizerWrtParameter.FisherInformationSource.EMPIRICAL_SCORES,
                maximizer.getFisherMetricSource());
        assertEquals(0, maximizer.getFisherNegativeEigenvalueCount());
        assertTrue(Double.isFinite(maximizer.getLikelihood().getLogLikelihood()));
    }

    public void testFiniteDifferenceSourceIgnoresScoreRows() {
        final MaximizerWrtParameter maximizer = runQuadraticRecovery(
                MaximizerWrtParameter.FisherInformationSource.FINITE_DIFFERENCE);

        assertEquals(1, maximizer.getFisherRefreshCount());
        assertEquals(2, maximizer.getFisherRank());
        assertEquals(0, maximizer.getFisherNullity());
        assertFalse(maximizer.fisherUsedEmpiricalScores());
        assertEquals(MaximizerWrtParameter.FisherInformationSource.FINITE_DIFFERENCE,
                maximizer.getFisherMetricSource());
        assertEquals(0, maximizer.getFisherNegativeEigenvalueCount());
        // Central differences of a linear gradient recover the exact Hessian diag(1e8, 1).
        assertEquals(1.0e8, maximizer.getFisherMaximumEigenvalue(), 1.0e8 * 1.0e-6);
        assertEquals(1.0, maximizer.getFisherMinimumPositiveEigenvalue(), 1.0e-6);
        assertTrue(maximizer.getFisherLastRefreshSeconds() >= 0.0);
        assertTrue(Double.isFinite(maximizer.getLikelihood().getLogLikelihood()));
    }

    public void testParseInformationSourceNames() {
        assertEquals(MaximizerWrtParameter.FisherInformationSource.FINITE_DIFFERENCE,
                MaximizerWrtParameter.FisherInformationSource.parse("finiteDifference"));
        assertEquals(MaximizerWrtParameter.FisherInformationSource.EMPIRICAL_SCORES,
                MaximizerWrtParameter.FisherInformationSource.parse("empiricalScores"));
        try {
            MaximizerWrtParameter.FisherInformationSource.parse("expected");
            fail("unknown source accepted");
        } catch (IllegalArgumentException expected) {
            // expected
        }
    }

    private static MaximizerWrtParameter runQuadraticRecovery(
            final MaximizerWrtParameter.FisherInformationSource source) {
        final Parameter parameter = new Parameter.Default("x", new double[]{1.0, 1.0});
        final Likelihood likelihood = new Likelihood.Abstract(null) {
            @Override
            protected boolean getLikelihoodKnown() {
                return false;
            }

            @Override
            protected double calculateLogLikelihood() {
                final double x0 = parameter.getParameterValue(0);
                final double x1 = parameter.getParameterValue(1);
                return -0.5 * (1.0e8 * x0 * x0 + x1 * x1);
            }
        };
        final ScoreContributionsWrtParameterProvider gradient = new ScoreContributionsWrtParameterProvider() {
            @Override
            public Likelihood getLikelihood() {
                return likelihood;
            }

            @Override
            public Parameter getParameter() {
                return parameter;
            }

            @Override
            public int getDimension() {
                return 2;
            }

            @Override
            public double[] getGradientLogDensity() {
                return new double[]{
                        -1.0e8 * parameter.getParameterValue(0),
                        -parameter.getParameterValue(1)
                };
            }

            @Override
            public double[][] getScoreContributionsLogDensity() {
                return new double[][]{
                        {1.0e4, 0.0},
                        {0.0, 1.0}
                };
            }
        };
        final MaximizerWrtParameter.Settings settings = new MaximizerWrtParameter.Settings(
                20, 1.0e-4, true, false, false,
                1.0e-30, 0, 1.0e-30, null, null, null, 1,
                true, 1, 1.0e-8, 10.0, 1.0, Double.POSITIVE_INFINITY, source);
        final MaximizerWrtParameter maximizer = new MaximizerWrtParameter(
                likelihood, parameter, gradient, null, settings, gradient, null);

        maximizer.maximize();
        return maximizer;
    }
}
