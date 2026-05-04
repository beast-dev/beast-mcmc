package test.dr.inference.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

public class GradientWrtParameterProviderDirtyRestoreTest extends TestCase {

    private static final class CachedQuadraticLikelihood extends Likelihood.Abstract {

        private final Parameter parameter;

        private CachedQuadraticLikelihood(final Parameter parameter) {
            super(null);
            this.parameter = parameter;
        }

        @Override
        protected double calculateLogLikelihood() {
            final double x = parameter.getParameterValue(0);
            return x * x;
        }
    }

    private static final class QuadraticGradient implements GradientWrtParameterProvider {

        private final Parameter parameter;
        private final Likelihood likelihood;

        private QuadraticGradient(final Parameter parameter, final Likelihood likelihood) {
            this.parameter = parameter;
            this.likelihood = likelihood;
        }

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
            return parameter.getDimension();
        }

        @Override
        public double[] getGradientLogDensity() {
            return new double[]{2.0 * parameter.getParameterValue(0)};
        }
    }

    public void testNumericalGradientDirtiesLikelihoodAndRestore() {
        final Parameter.Default parameter = new Parameter.Default(2.0);
        final CachedQuadraticLikelihood likelihood = new CachedQuadraticLikelihood(parameter);
        final QuadraticGradient gradient = new QuadraticGradient(parameter, likelihood);

        GradientWrtParameterProvider.getReportAndCheckForError(
                gradient,
                Double.NEGATIVE_INFINITY,
                Double.POSITIVE_INFINITY,
                1.0e-4);

        assertEquals(2.0, parameter.getParameterValue(0), 0.0);
        assertEquals(4.0, likelihood.getLogLikelihood(), 1.0e-12);
    }
}
