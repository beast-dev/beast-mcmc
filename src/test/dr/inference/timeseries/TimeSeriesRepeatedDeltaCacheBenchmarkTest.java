package test.dr.inference.timeseries;

import com.sun.management.ThreadMXBean;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.inference.model.GivensRotationMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.OrthogonalBlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.timeseries.core.BasicTimeSeriesModel;
import dr.inference.timeseries.core.IrregularTimeGrid;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.core.UniformTimeGrid;
import dr.inference.timeseries.engine.gaussian.GaussianForwardComputationMode;
import dr.inference.timeseries.gaussian.GaussianObservationModel;
import dr.inference.timeseries.gaussian.OUTimeSeriesProcessAdapter;
import dr.inference.timeseries.likelihood.GaussianGradientComputationMode;
import dr.inference.timeseries.likelihood.GaussianSmootherComputationMode;
import dr.inference.timeseries.likelihood.GaussianTimeSeriesLikelihoodFactory;
import dr.inference.timeseries.likelihood.ParallelTimeSeriesLikelihood;
import dr.inference.timeseries.likelihood.TimeSeriesLikelihood;
import dr.inference.timeseries.representation.CachedGaussianTransitionRepresentation;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;
import dr.inference.timeseries.representation.RepeatedDeltaCacheStatistics;
import junit.framework.TestCase;

import java.lang.management.ManagementFactory;
import java.util.ArrayList;
import java.util.List;

/**
 * Focused performance guard for repeated-delta OU time-series workloads.
 */
public class TimeSeriesRepeatedDeltaCacheBenchmarkTest extends TestCase {

    private static final int SERIES_COUNT = 12;
    private static final int TIME_COUNT = 80;
    private static final int DIM = 2;
    private static final double TOL = 1e-8;

    public void testRepeatedDtCacheBenchmarkDenseAndOrthogonalBlock() {
        runScenario("dense", makeDenseProcess("dense.bench"));
        runScenario("orthogonalBlock", makeOrthogonalBlockProcess("orthogonal.bench"));
    }

    public void testHeterogeneousDtCacheIsSharedAcrossOuAdapters() {
        final OUProcessModel process = makeOrthogonalBlockProcess("orthogonal.heterogeneous");
        final TimeGrid grid = new IrregularTimeGrid(new double[]{
                0.0, 0.07, 0.151, 0.244, 0.349, 0.466, 0.595, 0.736, 0.889, 1.054});
        final OUTimeSeriesProcessAdapter first = new OUTimeSeriesProcessAdapter(process);
        final OUTimeSeriesProcessAdapter second = new OUTimeSeriesProcessAdapter(process);
        final CachedGaussianTransitionRepresentation firstCache =
                (CachedGaussianTransitionRepresentation)
                        first.getRepresentation(GaussianTransitionRepresentation.class);
        final CachedGaussianTransitionRepresentation secondCache =
                (CachedGaussianTransitionRepresentation)
                        second.getRepresentation(GaussianTransitionRepresentation.class);
        assertSame("Adapters wrapping the same OU process should share transition cache",
                firstCache, secondCache);

        firstCache.prepareTimeGrid(grid);
        final CanonicalGaussianTransition transition = new CanonicalGaussianTransition(DIM);
        for (int t = 0; t < grid.getTimeCount() - 1; ++t) {
            firstCache.getCanonicalTransition(t, t + 1, grid, transition);
        }
        final RepeatedDeltaCacheStatistics afterFirstSeries = firstCache.getCacheStatistics();
        assertEquals("Heterogeneous grid should keep one entry per unique dt",
                grid.getTimeCount() - 1, afterFirstSeries.entryCount);
        assertEquals("First traversal should build each unique transition once",
                grid.getTimeCount() - 1, afterFirstSeries.canonicalBuilds);

        for (int t = 0; t < grid.getTimeCount() - 1; ++t) {
            secondCache.getCanonicalTransition(t, t + 1, grid, transition);
        }
        final RepeatedDeltaCacheStatistics afterSecondSeries = secondCache.getCacheStatistics();
        assertEquals("Second adapter should reuse all heterogeneous-grid transitions",
                afterFirstSeries.canonicalBuilds, afterSecondSeries.canonicalBuilds);
        assertTrue("Shared heterogeneous-grid cache should report cross-series hits",
                afterSecondSeries.canonicalHits() >= grid.getTimeCount() - 1);
    }

    public void testAllocationGuardsAfterWarmup() {
        final OUProcessModel process = makeOrthogonalBlockProcess("orthogonal.alloc");
        final OUTimeSeriesProcessAdapter adapter = new OUTimeSeriesProcessAdapter(process);
        final TimeGrid grid = new UniformTimeGrid(TIME_COUNT, 0.0, 0.125);
        final CachedGaussianTransitionRepresentation cached =
                (CachedGaussianTransitionRepresentation)
                        adapter.getRepresentation(GaussianTransitionRepresentation.class);
        cached.prepareTimeGrid(grid);

        final CanonicalGaussianTransition transition = new CanonicalGaussianTransition(DIM);
        cached.getCanonicalTransition(0, 1, grid, transition);
        final long transitionBytes = allocatedBytesFor(new Runnable() {
            @Override
            public void run() {
                for (int i = 0; i < 5000; ++i) {
                    cached.getCanonicalTransition(0, 1, grid, transition);
                }
            }
        });
        assertAllocationBelow("warm repeated transition lookup", transitionBytes, 64 * 1024L);

        final List<TimeSeriesLikelihood> likelihoods = makeLikelihoods(
                "alloc.canonical", adapter, grid, GaussianForwardComputationMode.CANONICAL,
                GaussianSmootherComputationMode.CANONICAL,
                GaussianGradientComputationMode.CANONICAL_ANALYTICAL);
        final TimeSeriesLikelihood likelihood = likelihoods.get(0);
        likelihood.getLogLikelihood();
        likelihood.getGradientWrt(process.getDriftMatrix()).getGradientLogDensity(null);

        final long likelihoodBytes = allocatedBytesFor(new Runnable() {
            @Override
            public void run() {
                likelihood.makeDirty();
                likelihood.getLogLikelihood();
            }
        });
        assertAllocationBelow("canonical likelihood pass", likelihoodBytes, 512 * 1024L);

        final Parameter angle = ((GivensRotationMatrixParameter)
                ((OrthogonalBlockDiagonalPolarStableMatrixParameter) process.getDriftMatrix())
                        .getRotationMatrixParameter()).getOrthogonalParameter();
        final long gradientBytes = allocatedBytesFor(new Runnable() {
            @Override
            public void run() {
                likelihood.makeDirty();
                likelihood.getGradientWrt(angle).getGradientLogDensity(null);
            }
        });
        assertAllocationBelow("canonical orthogonal gradient pass", gradientBytes, 2 * 1024 * 1024L);

        final ParallelTimeSeriesLikelihood aggregate =
                new ParallelTimeSeriesLikelihood("alloc.parallel", 4, likelihoods);
        aggregate.getGradientWrt(angle).getGradientLogDensity(null);
        final long parallelGradientBytes = allocatedBytesFor(new Runnable() {
            @Override
            public void run() {
                aggregate.makeDirty();
                aggregate.getGradientWrt(angle).getGradientLogDensity(null);
            }
        });
        assertAllocationBelow("parallel aggregate gradient summation", parallelGradientBytes, 6 * 1024 * 1024L);
    }

    private static void runScenario(final String label, final OUProcessModel process) {
        final TimeGrid grid = new UniformTimeGrid(TIME_COUNT, 0.0, 0.125);

        final OUTimeSeriesProcessAdapter canonicalAdapter = new OUTimeSeriesProcessAdapter(process);
        final List<TimeSeriesLikelihood> canonicalLikelihoods = makeLikelihoods(
                label + ".canonical", canonicalAdapter, grid,
                GaussianForwardComputationMode.CANONICAL,
                GaussianSmootherComputationMode.CANONICAL,
                GaussianGradientComputationMode.CANONICAL_ANALYTICAL);
        final ParallelTimeSeriesLikelihood parallelCanonical =
                new ParallelTimeSeriesLikelihood(label + ".parallelCanonical", 4, canonicalLikelihoods);
        final TimedValue serialCanonical = timedSerial(canonicalLikelihoods);
        parallelCanonical.makeDirty();
        final TimedValue parallelCanonicalValue = timedLikelihood(parallelCanonical);
        assertEquals(serialCanonical.value, parallelCanonicalValue.value, TOL);

        final CachedGaussianTransitionRepresentation canonicalCache =
                (CachedGaussianTransitionRepresentation) canonicalAdapter
                        .getRepresentation(GaussianTransitionRepresentation.class);
        final RepeatedDeltaCacheStatistics canonicalStats = canonicalCache.getCacheStatistics();
        assertEquals("Repeated canonical grid should have one dt entry", 1, canonicalStats.entryCount);
        assertEquals("Repeated canonical grid should build one transition", 1L, canonicalStats.canonicalBuilds);
        assertTrue("Repeated canonical grid should hit cache", canonicalStats.canonicalHits() > 0L);

        final OUTimeSeriesProcessAdapter expectationAdapter = new OUTimeSeriesProcessAdapter(process);
        final List<TimeSeriesLikelihood> expectationLikelihoods = makeLikelihoods(
                label + ".expectation", expectationAdapter, grid,
                GaussianForwardComputationMode.EXPECTATION,
                GaussianSmootherComputationMode.EXPECTATION,
                GaussianGradientComputationMode.EXPECTATION_ANALYTICAL);
        final TimedValue serialExpectation = timedSerial(expectationLikelihoods);
        assertEquals(serialExpectation.value, serialCanonical.value, 1e-7);

        final RepeatedDeltaCacheStatistics expectationStats =
                ((CachedGaussianTransitionRepresentation) expectationAdapter
                        .getRepresentation(GaussianTransitionRepresentation.class)).getCacheStatistics();
        assertEquals("Repeated expectation grid should have one dt entry", 1, expectationStats.entryCount);
        assertEquals("Repeated expectation grid should build one moment triple", 1L, expectationStats.momentBuilds);
        assertTrue("Repeated expectation grid should hit cache", expectationStats.momentHits() > 0L);

        System.out.println("timeseries repeated-dt benchmark [" + label + "] serialCanonicalNs="
                + serialCanonical.elapsedNanos
                + " parallelCanonicalNs=" + parallelCanonicalValue.elapsedNanos
                + " serialExpectationNs=" + serialExpectation.elapsedNanos
                + " canonicalRequests=" + canonicalStats.canonicalRequests
                + " canonicalBuilds=" + canonicalStats.canonicalBuilds
                + " momentRequests=" + expectationStats.momentRequests
                + " momentBuilds=" + expectationStats.momentBuilds);
    }

    private static List<TimeSeriesLikelihood> makeLikelihoods(final String prefix,
                                                              final OUTimeSeriesProcessAdapter adapter,
                                                              final TimeGrid grid,
                                                              final GaussianForwardComputationMode forwardMode,
                                                              final GaussianSmootherComputationMode smootherMode,
                                                              final GaussianGradientComputationMode gradientMode) {
        final List<TimeSeriesLikelihood> likelihoods = new ArrayList<TimeSeriesLikelihood>(SERIES_COUNT);
        for (int s = 0; s < SERIES_COUNT; ++s) {
            final GaussianObservationModel obs = makeObservation(prefix + ".obs" + s, s);
            final BasicTimeSeriesModel model = new BasicTimeSeriesModel(
                    prefix + ".model" + s, adapter, obs, grid);
            likelihoods.add(GaussianTimeSeriesLikelihoodFactory.create(
                    prefix + ".likelihood" + s,
                    model,
                    forwardMode,
                    smootherMode,
                    gradientMode));
        }
        return likelihoods;
    }

    private static GaussianObservationModel makeObservation(final String name, final int seriesIndex) {
        final MatrixParameter design = makeMatrix(name + ".H", new double[][]{
                {1.0, 0.0},
                {0.0, 1.0}
        });
        final MatrixParameter noise = makeMatrix(name + ".R", new double[][]{
                {0.35, 0.02},
                {0.02, 0.42}
        });
        final double[][] y = new double[DIM][TIME_COUNT];
        for (int t = 0; t < TIME_COUNT; ++t) {
            y[0][t] = 0.2 * Math.sin(0.15 * t + 0.1 * seriesIndex)
                    + 0.03 * seriesIndex;
            y[1][t] = 0.15 * Math.cos(0.11 * t - 0.07 * seriesIndex)
                    - 0.02 * seriesIndex;
        }
        return new GaussianObservationModel(name, DIM, design, noise, makeMatrix(name + ".Y", y));
    }

    private static TimedValue timedSerial(final List<TimeSeriesLikelihood> likelihoods) {
        for (final TimeSeriesLikelihood likelihood : likelihoods) {
            likelihood.makeDirty();
        }
        final long start = System.nanoTime();
        double value = 0.0;
        for (final TimeSeriesLikelihood likelihood : likelihoods) {
            value += likelihood.getLogLikelihood();
        }
        return new TimedValue(value, System.nanoTime() - start);
    }

    private static TimedValue timedLikelihood(final ParallelTimeSeriesLikelihood likelihood) {
        final long start = System.nanoTime();
        final double value = likelihood.getLogLikelihood();
        return new TimedValue(value, System.nanoTime() - start);
    }

    private static OUProcessModel makeDenseProcess(final String name) {
        return new OUProcessModel(name, DIM,
                makeMatrix(name + ".A", new double[][]{
                        {0.62, -0.08},
                        {0.04, 0.77}
                }),
                makeMatrix(name + ".Q", new double[][]{
                        {1.10, 0.06},
                        {0.06, 0.92}
                }),
                new Parameter.Default(name + ".mu", new double[]{0.12, -0.09}),
                makeMatrix(name + ".P0", new double[][]{
                        {0.90, 0.03},
                        {0.03, 0.95}
                }));
    }

    private static OUProcessModel makeOrthogonalBlockProcess(final String name) {
        final Parameter angles = new Parameter.Default(name + ".angles", new double[]{0.23});
        final GivensRotationMatrixParameter rotation =
                new GivensRotationMatrixParameter(name + ".rotation", DIM, angles);
        final Parameter scalar = new Parameter.Default(0);
        final Parameter rho = new Parameter.Default(name + ".rho", new double[]{0.72});
        final Parameter theta = new Parameter.Default(name + ".theta", new double[]{0.31});
        final Parameter t = new Parameter.Default(name + ".t", new double[]{-0.05});
        final OrthogonalBlockDiagonalPolarStableMatrixParameter drift =
                new OrthogonalBlockDiagonalPolarStableMatrixParameter(
                        name + ".Ablock", rotation, scalar, rho, theta, t);
        return new OUProcessModel(name, DIM,
                drift,
                makeMatrix(name + ".Q", new double[][]{
                        {1.20, 0.04},
                        {0.04, 0.86}
                }),
                new Parameter.Default(name + ".mu", new double[]{0.08, -0.11}),
                makeMatrix(name + ".P0", new double[][]{
                        {0.95, 0.02},
                        {0.02, 0.90}
                }));
    }

    private static MatrixParameter makeMatrix(final String name, final double[][] values) {
        final MatrixParameter matrix = new MatrixParameter(name, values.length, values[0].length);
        for (int i = 0; i < values.length; ++i) {
            for (int j = 0; j < values[0].length; ++j) {
                matrix.setParameterValue(i, j, values[i][j]);
            }
        }
        return matrix;
    }

    private static long allocatedBytesFor(final Runnable runnable) {
        final java.lang.management.ThreadMXBean baseBean = ManagementFactory.getThreadMXBean();
        if (!(baseBean instanceof ThreadMXBean)) {
            runnable.run();
            return -1L;
        }
        final ThreadMXBean bean = (ThreadMXBean) baseBean;
        if (!bean.isThreadAllocatedMemorySupported()) {
            runnable.run();
            return -1L;
        }
        final boolean enabled = bean.isThreadAllocatedMemoryEnabled();
        if (!enabled) {
            bean.setThreadAllocatedMemoryEnabled(true);
        }
        final long threadId = Thread.currentThread().getId();
        final long before = bean.getThreadAllocatedBytes(threadId);
        runnable.run();
        final long after = bean.getThreadAllocatedBytes(threadId);
        return after - before;
    }

    private static void assertAllocationBelow(final String label,
                                              final long allocatedBytes,
                                              final long threshold) {
        if (allocatedBytes < 0L) {
            System.out.println("allocation guard [" + label + "] unavailable on this JVM");
            return;
        }
        System.out.println("allocation guard [" + label + "] bytes=" + allocatedBytes);
        assertTrue(label + " allocated " + allocatedBytes + " bytes, threshold=" + threshold,
                allocatedBytes <= threshold);
    }

    private static final class TimedValue {
        final double value;
        final long elapsedNanos;

        TimedValue(final double value, final long elapsedNanos) {
            this.value = value;
            this.elapsedNanos = elapsedNanos;
        }
    }
}
