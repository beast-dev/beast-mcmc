package test.dr.evomodel.treedatalikelihood.discrete;

import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.ScoreContributionsWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

public final class MaximizerFisherRecoveryTest extends TestCase {

    public void testFailureTriggersAtMostConfiguredRefreshes() {
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
                true, 1, 1.0e-8, 10.0, 1.0, Double.POSITIVE_INFINITY);
        final MaximizerWrtParameter maximizer = new MaximizerWrtParameter(
                likelihood, parameter, gradient, null, settings, gradient, null);

        maximizer.maximize();

        assertEquals(1, maximizer.getFisherRefreshCount());
        assertEquals(2, maximizer.getFisherRank());
        assertEquals(0, maximizer.getFisherNullity());
        assertTrue(maximizer.fisherUsedEmpiricalScores());
        assertTrue(Double.isFinite(likelihood.getLogLikelihood()));
    }
}
