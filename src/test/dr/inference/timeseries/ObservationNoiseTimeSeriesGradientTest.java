package test.dr.inference.timeseries;

import static test.dr.inference.timeseries.OUTimeSeriesTestSupport.latent;
import static test.dr.inference.timeseries.OUTimeSeriesTestSupport.representation;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.inference.model.DiagonalMatrix;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.model.ScaledIdentityMatrixParameter;
import dr.inference.timeseries.core.BasicTimeSeriesModel;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.core.UniformTimeGrid;
import dr.inference.timeseries.engine.kalman.GaussianForwardComputationMode;
import dr.inference.timeseries.engine.kalman.MomentKalmanLikelihoodEngine;
import dr.inference.timeseries.likelihood.GaussianGradientComputationMode;
import dr.inference.timeseries.likelihood.GaussianSmootherComputationMode;
import dr.inference.timeseries.likelihood.GaussianTimeSeriesLikelihoodFactory;
import dr.inference.timeseries.likelihood.TimeSeriesLikelihood;
import dr.inference.timeseries.model.gaussian.LinearGaussianObservationModel;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

public class ObservationNoiseTimeSeriesGradientTest extends TestCase {

    private static final double FD_STEP = 1.0e-6;
    private static final double TOLERANCE = 3.0e-5;

    public ObservationNoiseTimeSeriesGradientTest(final String name) {
        super(name);
    }

    public void testMomentAnalyticalDiagonalObservationNoiseGradientMatchesFiniteDifference() {
        final TestModel model = makeModel("obs.noise.moment");
        final TimeSeriesLikelihood likelihood = makeLikelihood(
                "ts.obs.noise.moment",
                model,
                GaussianForwardComputationMode.MOMENT,
                GaussianSmootherComputationMode.MOMENT,
                GaussianGradientComputationMode.MOMENT_ANALYTICAL);

        assertGradientMatchesFiniteDifference("moment", model, likelihood);
    }

    public void testCanonicalAnalyticalDiagonalObservationNoiseGradientMatchesFiniteDifference() {
        final TestModel model = makeModel("obs.noise.canonical");
        final TimeSeriesLikelihood likelihood = makeLikelihood(
                "ts.obs.noise.canonical",
                model,
                GaussianForwardComputationMode.CANONICAL,
                GaussianSmootherComputationMode.CANONICAL,
                GaussianGradientComputationMode.CANONICAL_ANALYTICAL);

        assertGradientMatchesFiniteDifference("canonical", model, likelihood);
    }

    public void testMomentAnalyticalScalarIdentityObservationNoiseGradientMatchesFiniteDifference() {
        final TestModel model = makeScalarIdentityModel("obs.noise.scalar.moment");
        final TimeSeriesLikelihood likelihood = makeLikelihood(
                "ts.obs.noise.scalar.moment",
                model,
                GaussianForwardComputationMode.MOMENT,
                GaussianSmootherComputationMode.MOMENT,
                GaussianGradientComputationMode.MOMENT_ANALYTICAL);

        assertGradientMatchesFiniteDifference("moment scalar identity", model, likelihood);
    }

    public void testCanonicalAnalyticalScalarIdentityObservationNoiseGradientMatchesFiniteDifference() {
        final TestModel model = makeScalarIdentityModel("obs.noise.scalar.canonical");
        final TimeSeriesLikelihood likelihood = makeLikelihood(
                "ts.obs.noise.scalar.canonical",
                model,
                GaussianForwardComputationMode.CANONICAL,
                GaussianSmootherComputationMode.CANONICAL,
                GaussianGradientComputationMode.CANONICAL_ANALYTICAL);

        assertGradientMatchesFiniteDifference("canonical scalar identity", model, likelihood);
    }

    public void testCanonicalDriftScoreContributionsSumToGradient() {
        final TestModel model = makeModel("drift.score.canonical");
        final TimeSeriesLikelihood likelihood = makeLikelihood(
                "ts.drift.score.canonical",
                model,
                GaussianForwardComputationMode.CANONICAL,
                GaussianSmootherComputationMode.CANONICAL,
                GaussianGradientComputationMode.CANONICAL_ANALYTICAL);

        final double[] gradient = likelihood.getGradientWrt(model.drift)
                .getGradientLogDensity(null);
        final double[][] contributions = likelihood.getScoreContributionsWrt(model.drift);
        assertNotNull(contributions);
        assertEquals(model.grid.getTimeCount() - 1, contributions.length);
        final double[] sum = new double[gradient.length];
        for (final double[] contribution : contributions) {
            assertEquals(gradient.length, contribution.length);
            for (int i = 0; i < gradient.length; ++i) {
                sum[i] += contribution[i];
            }
        }
        for (int i = 0; i < gradient.length; ++i) {
            assertEquals("drift score contribution " + i, gradient[i], sum[i], 1.0e-10);
        }
    }

    private static void assertGradientMatchesFiniteDifference(final String label,
                                                              final TestModel model,
                                                              final TimeSeriesLikelihood likelihood) {
        likelihood.makeDirty();
        final double[] analytic = likelihood.getGradientWrt(model.observationVariance)
                .getGradientLogDensity(null);

        assertEquals(label + " observation-noise gradient length",
                model.observationVariance.getDimension(), analytic.length);
        for (int index = 0; index < analytic.length; ++index) {
            assertEquals(label + " observation-noise gradient entry " + index,
                    numericalGradient(model, likelihood, model.observationVariance, index),
                    analytic[index],
                    TOLERANCE);
        }
    }

    private static double numericalGradient(final TestModel model,
                                            final TimeSeriesLikelihood gradientLikelihood,
                                            final Parameter parameter,
                                            final int index) {
        final double saved = parameter.getParameterValue(index);

        parameter.setParameterValue(index, saved + FD_STEP);
        makeForwardLikelihoodDirty(model);
        final double plus = model.forwardLikelihood.getLogLikelihood();

        parameter.setParameterValue(index, saved - FD_STEP);
        makeForwardLikelihoodDirty(model);
        final double minus = model.forwardLikelihood.getLogLikelihood();

        parameter.setParameterValue(index, saved);
        makeForwardLikelihoodDirty(model);
        gradientLikelihood.makeDirty();

        return (plus - minus) / (2.0 * FD_STEP);
    }

    private static void makeForwardLikelihoodDirty(final TestModel model) {
        model.forwardLikelihood.makeDirty();
    }

    private static TimeSeriesLikelihood makeLikelihood(final String name,
                                                       final TestModel model,
                                                       final GaussianForwardComputationMode forwardMode,
                                                       final GaussianSmootherComputationMode smootherMode,
                                                       final GaussianGradientComputationMode gradientMode) {
        return GaussianTimeSeriesLikelihoodFactory.create(
                name,
                new BasicTimeSeriesModel(name + ".model", latent(model.process), model.observation, model.grid),
                forwardMode,
                smootherMode,
                gradientMode);
    }

    private static TestModel makeModel(final String prefix) {
        return makeModel(prefix, false);
    }

    private static TestModel makeScalarIdentityModel(final String prefix) {
        return makeModel(prefix, true);
    }

    private static TestModel makeModel(final String prefix, final boolean scalarIdentityObservationNoise) {
        final int dimension = 3;
        final MatrixParameter drift = makeMatrix(prefix + ".A", new double[][]{
                {0.52, -0.07, 0.04},
                {0.03, 0.86, -0.06},
                {-0.02, 0.05, 1.08}
        });
        final MatrixParameter diffusion = makeMatrix(prefix + ".Q", new double[][]{
                {0.90, 0.06, 0.02},
                {0.06, 0.82, 0.04},
                {0.02, 0.04, 1.04}
        });
        final Parameter mean = new Parameter.Default(prefix + ".mu", new double[]{0.12, -0.07, 0.18});
        final MatrixParameter initialCovariance = makeMatrix(prefix + ".P0", new double[][]{
                {0.95, 0.03, 0.01},
                {0.03, 1.02, 0.04},
                {0.01, 0.04, 0.88}
        });
        final OUProcessModel process = new OUProcessModel(
                prefix + ".ou",
                dimension,
                drift,
                diffusion,
                mean,
                initialCovariance,
                OUProcessModel.CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final MatrixParameter design = makeMatrix(prefix + ".H", new double[][]{
                {1.00, 0.12, 0.00},
                {0.00, 1.00, -0.10},
                {0.08, 0.00, 1.00}
        });
        final Parameter observationVariance = scalarIdentityObservationNoise
                ? new Parameter.Default(prefix + ".obsVariance", new double[]{0.29})
                : new Parameter.Default(prefix + ".obsVariance", new double[]{0.27, 0.31, 0.29});
        final MatrixParameter observationNoise = scalarIdentityObservationNoise
                ? new ScaledIdentityMatrixParameter(prefix + ".obsNoise", dimension, observationVariance)
                : new DiagonalMatrix(observationVariance);
        final MatrixParameter observations = makeMatrix(prefix + ".Y", new double[][]{
                {0.16, -0.24, 0.31, Double.NaN, 0.47},
                {-0.10, 0.29, -0.18, Double.NaN, 0.12},
                {0.24, -0.13, 0.08, Double.NaN, -0.20}
        });
        final LinearGaussianObservationModel observation = new LinearGaussianObservationModel(
                prefix + ".obs",
                dimension,
                design,
                observationNoise,
                observations);
        final TimeGrid grid = new UniformTimeGrid(observations.getColumnDimension(), 0.0, 0.17);
        final GaussianTransitionRepresentation transitionRepresentation =
                representation(process, GaussianTransitionRepresentation.class);
        final MomentKalmanLikelihoodEngine forwardLikelihood =
                new MomentKalmanLikelihoodEngine(transitionRepresentation, observation, grid);

        return new TestModel(process, drift, observationVariance, observation, grid, forwardLikelihood);
    }

    private static MatrixParameter makeMatrix(final String name, final double[][] values) {
        final MatrixParameter matrix = new MatrixParameter(name, values.length, values[0].length);
        for (int row = 0; row < values.length; ++row) {
            for (int col = 0; col < values[0].length; ++col) {
                matrix.setParameterValue(row, col, values[row][col]);
            }
        }
        return matrix;
    }

    public static Test suite() {
        return new TestSuite(ObservationNoiseTimeSeriesGradientTest.class);
    }

    private static final class TestModel {
        final OUProcessModel process;
        final MatrixParameter drift;
        final Parameter observationVariance;
        final LinearGaussianObservationModel observation;
        final TimeGrid grid;
        final MomentKalmanLikelihoodEngine forwardLikelihood;

        private TestModel(final OUProcessModel process,
                          final MatrixParameter drift,
                          final Parameter observationVariance,
                          final LinearGaussianObservationModel observation,
                          final TimeGrid grid,
                          final MomentKalmanLikelihoodEngine forwardLikelihood) {
            this.process = process;
            this.drift = drift;
            this.observationVariance = observationVariance;
            this.observation = observation;
            this.grid = grid;
            this.forwardLikelihood = forwardLikelihood;
        }
    }
}
