package dr.inference.timeseries.likelihood;

import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.GradientProvider;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.timeseries.engine.gaussian.SharedCanonicalTimeSeriesSchedule;

import java.util.ArrayList;
import java.util.Collection;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;

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
    private final int threadCount;
    private final ParallelEvaluator parallelEvaluator;
    private final Map<Parameter, GradientProvider[]> gradientProviderCache;
    private final SharedCanonicalTimeSeriesSchedule sharedCanonicalSchedule;

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
        this.sharedCanonicalSchedule = installSharedCanonicalSchedule(this.likelihoods);

        if (threads < 0 && this.likelihoods.size() > 1) {
            this.threadCount = this.likelihoods.size();
        } else if (threads > 1 && this.likelihoods.size() > 1) {
            this.threadCount = Math.min(threads, this.likelihoods.size());
        } else {
            this.threadCount = 0;
        }

        if (threadCount > 0) {
            this.parallelEvaluator = new ParallelEvaluator(threadCount);
        } else {
            this.parallelEvaluator = null;
        }
        this.gradientProviderCache = new IdentityHashMap<Parameter, GradientProvider[]>();
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

    public SharedCanonicalTimeSeriesSchedule getSharedCanonicalSchedule() {
        return sharedCanonicalSchedule;
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
        if (sharedCanonicalSchedule != null) {
            sharedCanonicalSchedule.makeDirty();
        }
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
        if (sharedCanonicalSchedule != null) {
            sharedCanonicalSchedule.makeDirty();
        }
        for (final TimeSeriesLikelihood likelihood : likelihoods) {
            likelihood.makeDirty();
        }
    }

    @Override
    protected void acceptState() {
        // no-op
    }

    private double computeLogLikelihood() {
        if (parallelEvaluator == null) {
            double value = 0.0;
            for (final TimeSeriesLikelihood likelihood : likelihoods) {
                value += likelihood.getLogLikelihood();
            }
            return value;
        }
        return parallelEvaluator.evaluateLogLikelihood();
    }

    private double[] computeGradient(final Parameter parameter) {
        final int dimension = parameter.getDimension();
        if (parallelEvaluator == null) {
            final double[] gradient = new double[dimension];
            for (int i = 0; i < likelihoods.size(); ++i) {
                addGradient(gradient, gradientFor(i, parameter), dimension);
            }
            return gradient;
        }
        return parallelEvaluator.evaluateGradient(parameter, dimension);
    }

    private static SharedCanonicalTimeSeriesSchedule installSharedCanonicalSchedule(
            final List<TimeSeriesLikelihood> likelihoods) {
        SharedCanonicalTimeSeriesSchedule schedule = null;
        int installed = 0;
        for (final TimeSeriesLikelihood likelihood : likelihoods) {
            if (schedule == null) {
                schedule = likelihood.createSharedCanonicalSchedule();
            }
            if (schedule != null && likelihood.setSharedCanonicalSchedule(schedule)) {
                ++installed;
            }
        }
        return installed > 1 ? schedule : null;
    }

    private double[] gradientFor(final int likelihoodIndex, final Parameter parameter) {
        final GradientProvider provider = gradientProviderFor(likelihoodIndex, parameter);
        if (provider.getDimension() != parameter.getDimension()) {
            throw new IllegalArgumentException("Gradient dimension mismatch for time-series likelihood "
                    + likelihoodIndex + ": expected " + parameter.getDimension()
                    + " but got " + provider.getDimension());
        }
        return provider.getGradientLogDensity(null);
    }

    private GradientProvider gradientProviderFor(final int likelihoodIndex,
                                                 final Parameter parameter) {
        GradientProvider[] providers = gradientProviderCache.get(parameter);
        if (providers == null) {
            providers = new GradientProvider[likelihoods.size()];
            gradientProviderCache.put(parameter, providers);
        }
        GradientProvider provider = providers[likelihoodIndex];
        if (provider == null) {
            provider = likelihoods.get(likelihoodIndex).getGradientWrt(parameter);
            providers[likelihoodIndex] = provider;
        }
        return provider;
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

    private final class ParallelEvaluator {

        private static final int MODE_IDLE = 0;
        private static final int MODE_LIKELIHOOD = 1;
        private static final int MODE_GRADIENT = 2;

        private final Object monitor = new Object();
        private final Thread[] workers;

        private int generation;
        private int mode;
        private int nextIndex;
        private int completed;
        private double logLikelihoodAccumulator;
        private double[] gradientAccumulator;
        private Parameter gradientParameter;
        private RuntimeException failure;

        private ParallelEvaluator(final int workerCount) {
            this.workers = new Thread[workerCount];
            this.mode = MODE_IDLE;
            for (int i = 0; i < workerCount; ++i) {
                final Thread worker = new Thread(new Worker(), getId() + "-worker-" + (i + 1));
                worker.setDaemon(true);
                worker.start();
                workers[i] = worker;
            }
        }

        private double evaluateLogLikelihood() {
            synchronized (monitor) {
                startEvaluation(MODE_LIKELIHOOD, null, null);
                waitForCompletion("Interrupted while evaluating parallel time-series likelihood");
                if (failure != null) {
                    throw new RuntimeException("Failed to evaluate parallel time-series likelihood", failure);
                }
                return logLikelihoodAccumulator;
            }
        }

        private double[] evaluateGradient(final Parameter parameter, final int dimension) {
            final double[] gradient = new double[dimension];
            synchronized (monitor) {
                startEvaluation(MODE_GRADIENT, parameter, gradient);
                waitForCompletion("Interrupted while evaluating parallel time-series gradient");
                if (failure != null) {
                    throw new RuntimeException("Failed to evaluate parallel time-series gradient", failure);
                }
                return gradient;
            }
        }

        private void startEvaluation(final int evaluationMode,
                                     final Parameter parameter,
                                     final double[] gradient) {
            mode = evaluationMode;
            nextIndex = 0;
            completed = 0;
            logLikelihoodAccumulator = 0.0;
            gradientAccumulator = gradient;
            gradientParameter = parameter;
            failure = null;
            ++generation;
            monitor.notifyAll();
        }

        private void waitForCompletion(final String interruptedMessage) {
            while (completed < likelihoods.size() && failure == null) {
                try {
                    monitor.wait();
                } catch (final InterruptedException e) {
                    Thread.currentThread().interrupt();
                    throw new RuntimeException(interruptedMessage, e);
                }
            }
            mode = MODE_IDLE;
            gradientParameter = null;
            gradientAccumulator = null;
        }

        private final class Worker implements Runnable {

            private int lastGeneration;

            @Override
            public void run() {
                while (true) {
                    final int activeGeneration = waitForWork();
                    runEvaluation(activeGeneration);
                }
            }

            private int waitForWork() {
                synchronized (monitor) {
                    while (lastGeneration == generation) {
                        try {
                            monitor.wait();
                        } catch (final InterruptedException ignored) {
                            // Daemon worker: keep serving later likelihood evaluations.
                        }
                    }
                    lastGeneration = generation;
                    return generation;
                }
            }

            private void runEvaluation(final int activeGeneration) {
                while (true) {
                    final int index;
                    final int activeMode;
                    final Parameter activeParameter;
                    synchronized (monitor) {
                        if (activeGeneration != generation || failure != null || nextIndex >= likelihoods.size()) {
                            return;
                        }
                        index = nextIndex++;
                        activeMode = mode;
                        activeParameter = gradientParameter;
                    }
                    try {
                        if (activeMode == MODE_LIKELIHOOD) {
                            final double value = likelihoods.get(index).getLogLikelihood();
                            synchronized (monitor) {
                                if (activeGeneration == generation && failure == null) {
                                    logLikelihoodAccumulator += value;
                                    completeOne();
                                }
                            }
                        } else if (activeMode == MODE_GRADIENT) {
                            final double[] gradient = gradientFor(index, activeParameter);
                            synchronized (monitor) {
                                if (activeGeneration == generation && failure == null) {
                                    addGradient(gradientAccumulator, gradient, gradientAccumulator.length);
                                    completeOne();
                                }
                            }
                        }
                    } catch (final RuntimeException e) {
                        synchronized (monitor) {
                            if (activeGeneration == generation && failure == null) {
                                failure = e;
                                monitor.notifyAll();
                            }
                        }
                        return;
                    }
                }
            }

            private void completeOne() {
                ++completed;
                if (completed >= likelihoods.size()) {
                    monitor.notifyAll();
                }
            }
        }
    }
}
