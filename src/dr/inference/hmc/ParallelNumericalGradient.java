package dr.inference.hmc;

import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.NumericalDerivative;
import dr.math.UnivariateFunction;
import dr.math.matrixAlgebra.WrappedVector;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

/**
 * @author Marc A Suchard
 */
public class ParallelNumericalGradient implements GradientWrtParameterProvider, Reportable {

    public ParallelNumericalGradient(List<Likelihood> likelihoods,
                                     List<Parameter> parameters) {
        this.likelihoods = likelihoods;
        this.parameters = parameters;
        this.threadCount = likelihoods.size();

        pool = Executors.newFixedThreadPool(threadCount);

//        ensureState();
    }

    @SuppressWarnings("unused")
    private void ensureState() {
        double logLik = likelihoods.get(0).getLogLikelihood();

        for (int i = 1; i < likelihoods.size(); ++i) {
            double test = likelihoods.get(i).getLogLikelihood();
            if (Math.abs(logLik - test) > tolerance) {
                throw new RuntimeException("Invalid state");
            }
        }

        double[] values = parameters.get(0).getParameterValues();

        for (int i = 1; i < parameters.size(); ++i) {
            double[] test = parameters.get(i).getParameterValues();
            if (notEqual(values, test)) {
                throw new RuntimeException("Invalid state");
            }
        }
    }

    private static final double tolerance = 1E-10;

    boolean notEqual(double[] lhs, double[] rhs) {
        for (int i = 0; i < lhs.length; ++i) {
            if (Math.abs(lhs[i] - rhs[i]) > tolerance) {
                return true;
            }
        }
        return false;
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihoods.get(0);
    }

    @Override
    public Parameter getParameter() {
        return parameters.get(0);
    }

    @Override
    public int getDimension() {
        return parameters.get(0).getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        final double[] x = parameters.get(0).getParameterValues();

        final int len = x.length;

        double[] gradient = new double[len];

        final int stride = (len + 1) / threadCount;

        List<SubsetCaller> callers = new ArrayList<>();
        for (int thread = 0; thread < threadCount; ++thread) {

            int start = thread * stride;
            int end = Math.min(len, start + stride);

            callers.add(new SubsetCaller(
                    likelihoods.get(thread),
                    parameters.get(thread),
                    gradient,
                    thread, x,
                    start, end));

        }

        try {
            pool.invokeAll(callers);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

//        ensureState(); // TODO remove

        return gradient;
    }

    @Override
    public String getReport() {
        long start = System.currentTimeMillis();
        double[] result = getGradientLogDensity();
        long duration = System.currentTimeMillis() - start;

        return "Numeric gradient = " + new WrappedVector.Raw(result) +
                " (" + threadCount + " threads in " + duration + "ms)";
    }


    private static class SubsetCaller implements Callable<Integer> {

        public SubsetCaller(Likelihood likelihood,
                            Parameter parameter,
                            double[] gradient,
                            int index,
                            double[] x,
                            int start, int end) {

            this.likelihood = likelihood;
            this.parameter = parameter;
            this.gradient = gradient;
            this.index = index;
            this.start = start;
            this.end = end;

            for (int i = 0; i < parameter.getDimension(); ++i) {
                parameter.setParameterValueQuietly(i, x[i]);
            }
            parameter.fireParameterChangedEvent();
        }

        @Override
        public Integer call() throws Exception {

            if (DEBUG_PARALLEL_EVALUATION) {
                System.err.println("Invoking thread #" + index + " for " + likelihood.getId());
            }

            for (int idx = start; idx < end; ++idx) {

                final int k = idx;
                final double xk = parameter.getParameterValue(k);

                gradient[k] = NumericalDerivative.firstDerivative(new UnivariateFunction() {
                    // TODO Benchmark effect of false-sharing across gradient[]

                    @Override
                    public double evaluate(double x) {
                        double old = parameter.getParameterValue(k);
                        parameter.setParameterValue(k, x);
                        double dx = likelihood.getLogLikelihood();
                        parameter.setParameterValue(k, old);
                        return dx;
                    }

                    @Override
                    public double getLowerBound() {
                        return Double.NEGATIVE_INFINITY;
                    }

                    @Override
                    public double getUpperBound() {
                        return Double.POSITIVE_INFINITY;
                    }
                }, xk);
            }

            return 0;
        }

        private final Likelihood likelihood;
        private final Parameter parameter;
        private final double[] gradient;
        private final int index;
        private final int start;
        private final int end;
    }

    public static final boolean DEBUG_PARALLEL_EVALUATION = false;

    private final List<Likelihood> likelihoods;
    private final List<Parameter> parameters;

    private final ExecutorService pool;

    private final int threadCount;
}
