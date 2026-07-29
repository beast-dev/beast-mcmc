package test.dr.inference.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.hmc.HamiltonianMonteCarloOperator;
import dr.inference.operators.hmc.MassPreconditionScheduler;
import dr.inference.operators.hmc.MassPreconditioner;
import dr.inference.operators.hmc.MassPreconditioningOptions;
import junit.framework.TestCase;

public class HamiltonianMonteCarloOperatorTest extends TestCase {

    private static final class FlatLikelihood extends Likelihood.Abstract {

        private FlatLikelihood() {
            super(null);
        }

        @Override
        protected double calculateLogLikelihood() {
            return 0.0;
        }
    }

    private static final class DomainFailingGradient implements GradientWrtParameterProvider {

        private final Parameter parameter;
        private final Likelihood likelihood;

        private DomainFailingGradient(final Parameter parameter, final Likelihood likelihood) {
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
            throw new IllegalArgumentException("domain failure");
        }
    }

    public void testIllegalArgumentDuringLeapfrogRejectsProposal() {
        final Parameter.Default parameter = new Parameter.Default(1.0);
        final Likelihood likelihood = new FlatLikelihood();
        final DomainFailingGradient gradient = new DomainFailingGradient(parameter, likelihood);
        final MassPreconditioningOptions options =
                new MassPreconditioningOptions.Default(0, 0, 0, 0, false, null, null);
        final HamiltonianMonteCarloOperator.Options runtimeOptions =
                new HamiltonianMonteCarloOperator.Options(
                        0.1, 2, 0.0, options,
                        0, 1.0e-3, 10, 0.1, 0.8,
                        HamiltonianMonteCarloOperator.InstabilityHandler.REJECT);
        final HamiltonianMonteCarloOperator operator =
                new HamiltonianMonteCarloOperator(
                        AdaptationMode.ADAPTATION_OFF, 1.0, gradient, parameter, null, null,
                        runtimeOptions,
                        MassPreconditioner.Type.NONE.factory(gradient, null, options),
                        MassPreconditionScheduler.Type.NONE);

        assertEquals(Double.NEGATIVE_INFINITY, operator.doOperation(likelihood));
    }
}
