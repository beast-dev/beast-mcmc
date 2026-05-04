package dr.inference.timeseries.likelihood;

import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.GradientProvider;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Aggregate likelihood for independent time series that share model parameters.
 *
 * <p>Each child likelihood owns its own filtering/smoothing workspace. This class
 * only coordinates embarrassingly parallel evaluation across series and reduces
 * log likelihoods and per-parameter gradients by summation.</p>
 */
public final class ParallelTimeSeriesLikelihood extends AbstractModelLikelihood
        implements TimeSeriesGradientSource {

    private final List<TimeSeriesLikelihood> likelihoods;
    private final ExecutorService pool;
    private final int threadCount;

    private boolean likelihoodKnown;
    private double logLikelihood;

    private boolean storedLikelihoodKnown;
    private double storedLogLikelihood;

    public ParallelTimeSeriesLikelihood(final String name,
                                        final int threads,
                                        final Collection<TimeSeriesLikelihood> likelihoods) {
        super(name);
        if (likelihoods == null || likelihoods.isEmpty()) {
            throw new IllegalArgumentException("likelihoods must contain at least one time-series likelihood");
        }
        this.likelihoods = new ArrayList<TimeSeriesLikelihood>(likelihoods.size());
        for (final TimeSeriesLikelihood likelihood : likelihoods) {
            if (likelihood == null) {
                throw new IllegalArgumentException("likelihoods must not contain null entries");
            }
            this.likelihoods.add(likelihood);
            addModel(likelihood.getModel());
        }

        if (threads < 0 && this.likelihoods.size() > 1) {
            this.threadCount = this.likelihoods.size();
        } else if (threads > 1 && this.likelihoods.size() > 1) {
            this.threadCount = Math.min(threads, this.likelihoods.size());
        } else {
            this.threadCount = 0;
        }

        if (threadCount > 0) {
            this.pool = Executors.newFixedThreadPool(threadCount);
        } else {
            this.pool = null;
        }
        this.likelihoodKnown = false;
    }

    public ParallelTimeSeriesLikelihood(final String name,
                                        final Collection<TimeSeriesLikelihood> likelihoods) {
        this(name, 0, likelihoods);
    }

    @Override
    public Model getModel() {
        return this;
    }

    public int getTimeSeriesCount() {
        return likelihoods.size();
    }

    public int getThreadCount() {
        return threadCount;
    }

    public TimeSeriesLikelihood getTimeSeriesLikelihood(final int index) {
        return likelihoods.get(index);
    }

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = computeLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    @Override
    public GradientProvider getGradientWrt(final Parameter parameter) {
        if (parameter == null) {
            throw new IllegalArgumentException("parameter must not be null");
        }
        return new GradientProvider() {
            @Override
            public int getDimension() {
                return parameter.getDimension();
            }

            @Override
            public double[] getGradientLogDensity(final Object x) {
                return computeGradient(parameter);
            }
        };
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
        for (final TimeSeriesLikelihood likelihood : likelihoods) {
            likelihood.makeDirty();
        }
    }

    @Override
    protected void handleModelChangedEvent(final Model model, final Object object, final int index) {
        makeDirty();
    }

    @Override
    protected void handleVariableChangedEvent(final Variable variable,
                                              final int index,
                                              final Parameter.ChangeType type) {
        makeDirty();
    }

    @Override
    protected void storeState() {
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
    }

    @Override
    protected void restoreState() {
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
        for (final TimeSeriesLikelihood likelihood : likelihoods) {
            likelihood.makeDirty();
        }
    }

    @Override
    protected void acceptState() {
        // no-op
    }

    private double computeLogLikelihood() {
        if (pool == null) {
            double value = 0.0;
            for (final TimeSeriesLikelihood likelihood : likelihoods) {
                value += likelihood.getLogLikelihood();
            }
            return value;
        }

        final List<Callable<Double>> tasks = new ArrayList<Callable<Double>>(likelihoods.size());
        for (final TimeSeriesLikelihood likelihood : likelihoods) {
            tasks.add(new Callable<Double>() {
                @Override
                public Double call() {
                    return likelihood.getLogLikelihood();
                }
            });
        }
        try {
            double value = 0.0;
            for (final Future<Double> result : pool.invokeAll(tasks)) {
                value += result.get();
            }
            return value;
        } catch (final InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Interrupted while evaluating parallel time-series likelihood", e);
        } catch (final ExecutionException e) {
            throw new RuntimeException("Failed to evaluate parallel time-series likelihood", e.getCause());
        }
    }

    private double[] computeGradient(final Parameter parameter) {
        final int dimension = parameter.getDimension();
        if (pool == null) {
            final double[] gradient = new double[dimension];
            for (int i = 0; i < likelihoods.size(); ++i) {
                addGradient(gradient, gradientFor(i, parameter), dimension);
            }
            return gradient;
        }

        final List<Callable<double[]>> tasks = new ArrayList<Callable<double[]>>(likelihoods.size());
        for (int i = 0; i < likelihoods.size(); ++i) {
            final int index = i;
            tasks.add(new Callable<double[]>() {
                @Override
                public double[] call() {
                    return gradientFor(index, parameter);
                }
            });
        }
        try {
            final double[] gradient = new double[dimension];
            for (final Future<double[]> result : pool.invokeAll(tasks)) {
                addGradient(gradient, result.get(), dimension);
            }
            return gradient;
        } catch (final InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Interrupted while evaluating parallel time-series gradient", e);
        } catch (final ExecutionException e) {
            throw new RuntimeException("Failed to evaluate parallel time-series gradient", e.getCause());
        }
    }

    private double[] gradientFor(final int likelihoodIndex, final Parameter parameter) {
        final GradientProvider provider = likelihoods.get(likelihoodIndex).getGradientWrt(parameter);
        if (provider.getDimension() != parameter.getDimension()) {
            throw new IllegalArgumentException("Gradient dimension mismatch for time-series likelihood "
                    + likelihoodIndex + ": expected " + parameter.getDimension()
                    + " but got " + provider.getDimension());
        }
        return provider.getGradientLogDensity(null);
    }

    private static void addGradient(final double[] accumulator,
                                    final double[] increment,
                                    final int dimension) {
        if (increment == null || increment.length != dimension) {
            throw new IllegalArgumentException("Gradient increment must have dimension " + dimension);
        }
        for (int i = 0; i < dimension; ++i) {
            accumulator[i] += increment[i];
        }
    }
}
