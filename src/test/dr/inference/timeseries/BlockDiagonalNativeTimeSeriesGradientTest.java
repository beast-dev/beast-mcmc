package test.dr.inference.timeseries;

import static test.dr.inference.timeseries.OUTimeSeriesTestSupport.latent;
import static test.dr.inference.timeseries.OUTimeSeriesTestSupport.representation;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalFrechetHelper;
import dr.evomodel.treedatalikelihood.continuous.canonical.gradient.CanonicalSelectionGradientProjector;
import dr.inference.model.BlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.timeseries.core.BasicTimeSeriesModel;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.core.UniformTimeGrid;
import dr.inference.timeseries.engine.kalman.GaussianForwardComputationMode;
import dr.inference.timeseries.engine.kalman.KalmanLikelihoodEngine;
import dr.inference.timeseries.model.gaussian.GaussianObservationModel;
import dr.inference.timeseries.likelihood.GaussianGradientComputationMode;
import dr.inference.timeseries.likelihood.GaussianSmootherComputationMode;
import dr.inference.timeseries.likelihood.GaussianTimeSeriesLikelihoodFactory;
import dr.inference.timeseries.likelihood.TimeSeriesLikelihood;
import dr.inference.timeseries.representation.CachedGaussianTransitionRepresentation;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;
import junit.framework.TestCase;

public class BlockDiagonalNativeTimeSeriesGradientTest extends TestCase {

    private static final double FD_STEP = 1.0e-6;
    private static final double TOL = 8.0e-5;

    public BlockDiagonalNativeTimeSeriesGradientTest(final String name) {
        super(name);
    }

    public void testExpectationAnalyticalNativeGeneralBlockGradientsMatchFiniteDifference() {
        final GeneralBlockModel model = makeGeneralBlockModel("expectation.general.block", false);
        final TimeSeriesLikelihood likelihood = makeLikelihood(
                "ts.expectation.general.block",
                model,
                GaussianForwardComputationMode.EXPECTATION,
                GaussianSmootherComputationMode.EXPECTATION,
                GaussianGradientComputationMode.EXPECTATION_ANALYTICAL);

        assertNativeBlockGradientsMatchFiniteDifference("expectation", likelihood, model);
    }

    public void testCanonicalAnalyticalNativeGeneralBlockGradientsMatchFiniteDifference() {
        final GeneralBlockModel model = makeGeneralBlockModel("canonical.general.block", true);
        final TimeSeriesLikelihood likelihood = makeLikelihood(
                "ts.canonical.general.block",
                model,
                GaussianForwardComputationMode.CANONICAL,
                GaussianSmootherComputationMode.CANONICAL,
                GaussianGradientComputationMode.CANONICAL_ANALYTICAL);

        assertNativeBlockGradientsMatchFiniteDifference("canonical", likelihood, model);
    }

    public void testCanonicalAnalyticalNativeGeneralBlockDiffusionAndMeanGradientsMatchFiniteDifference() {
        final GeneralBlockModel model = makeGeneralBlockModel("canonical.general.block.qmu", true);
        final TimeSeriesLikelihood likelihood = makeLikelihood(
                "ts.canonical.general.block.qmu",
                model,
                GaussianForwardComputationMode.CANONICAL,
                GaussianSmootherComputationMode.CANONICAL,
                GaussianGradientComputationMode.CANONICAL_ANALYTICAL);

        assertGradientMatchesFiniteDifference("canonical diffusion", likelihood, model, model.diffusion);
        assertGradientMatchesFiniteDifference("canonical mean", likelihood, model, model.mean);
    }

    private static void assertNativeBlockGradientsMatchFiniteDifference(
            final String label,
            final TimeSeriesLikelihood likelihood,
            final GeneralBlockModel model) {
        CanonicalSelectionGradientProjector.resetInstrumentation();
        BlockDiagonalFrechetHelper.resetExactPlanInstrumentation();

        assertGradientMatchesFiniteDifference(label + " rho",
                likelihood, model, model.block.getRhoParameter());
        assertGradientMatchesFiniteDifference(label + " theta",
                likelihood, model, model.block.getThetaParameter());
        assertGradientMatchesFiniteDifference(label + " t",
                likelihood, model, model.block.getTParameter());
        assertGradientMatchesFiniteDifference(label + " R",
                likelihood, model, model.block.getRotationMatrixParameter());

        assertTrue(label + " native block path should use block Frechet adjoints",
                blockFrechetApplicationCount() > 0L);
        assertEquals(label + " native block path should avoid legacy dense matrix pullback",
                0L,
                CanonicalSelectionGradientProjector.getLegacyMatrixPullbackCount());
    }

    private static void assertGradientMatchesFiniteDifference(final String label,
                                                              final TimeSeriesLikelihood likelihood,
                                                              final GeneralBlockModel model,
                                                              final Parameter parameter) {
        likelihood.makeDirty();
        final double[] analytic = likelihood.getGradientWrt(parameter).getGradientLogDensity(null);
        assertEquals(label + " gradient length", parameter.getDimension(), analytic.length);
        for (int i = 0; i < analytic.length; ++i) {
            final double fd = numericalGradient(model, likelihood, parameter, i);
            assertEquals(label + " gradient entry " + i, fd, analytic[i], TOL);
        }
    }

    private static double numericalGradient(final GeneralBlockModel model,
                                            final TimeSeriesLikelihood gradientLikelihood,
                                            final Parameter parameter,
                                            final int index) {
        final double saved = parameter.getParameterValue(index);

        parameter.setParameterValue(index, saved + FD_STEP);
        makeForwardLikelihoodDirty(model);
        final double plus = model.likelihoodEngine.getLogLikelihood();

        parameter.setParameterValue(index, saved - FD_STEP);
        makeForwardLikelihoodDirty(model);
        final double minus = model.likelihoodEngine.getLogLikelihood();

        parameter.setParameterValue(index, saved);
        makeForwardLikelihoodDirty(model);
        gradientLikelihood.makeDirty();

        return (plus - minus) / (2.0 * FD_STEP);
    }

    private static void makeForwardLikelihoodDirty(final GeneralBlockModel model) {
        if (model.transitionRepresentation instanceof CachedGaussianTransitionRepresentation) {
            ((CachedGaussianTransitionRepresentation) model.transitionRepresentation).makeDirty();
        }
        model.likelihoodEngine.makeDirty();
    }

    private static TimeSeriesLikelihood makeLikelihood(final String name,
                                                       final GeneralBlockModel model,
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

    private static GeneralBlockModel makeGeneralBlockModel(final String prefix,
                                                           final boolean includeMissingObservation) {
        final int dimension = 4;
        final double dt = 0.19;
        final double[][] rotation = {
                {1.00, 0.20, 0.00, 0.10},
                {0.10, 1.10, 0.15, 0.00},
                {0.00, 0.05, 0.90, 0.25},
                {0.05, 0.00, 0.10, 1.05}
        };
        final BlockDiagonalPolarStableMatrixParameter block =
                new BlockDiagonalPolarStableMatrixParameter(
                        prefix + ".A",
                        makeMatrix(prefix + ".R", rotation),
                        new Parameter.Default(prefix + ".scalar", 0),
                        new Parameter.Default(prefix + ".rho", new double[]{0.70, 0.90}),
                        new Parameter.Default(prefix + ".theta", new double[]{0.25, -0.35}),
                        new Parameter.Default(prefix + ".t", new double[]{0.10, -0.08}));

        final MatrixParameter diffusion = makeMatrix(prefix + ".Q", new double[][]{
                {1.20, 0.10, 0.05, 0.02},
                {0.10, 0.95, 0.04, 0.03},
                {0.05, 0.04, 1.15, 0.07},
                {0.02, 0.03, 0.07, 0.90}
        });
        final Parameter mean = new Parameter.Default(prefix + ".mu", new double[]{0.30, -0.40, 0.20, 0.60});
        final MatrixParameter initialCovariance = makeMatrix(prefix + ".P0", new double[][]{
                {1.10, 0.10, 0.00, 0.00},
                {0.10, 1.20, 0.10, 0.00},
                {0.00, 0.10, 1.30, 0.10},
                {0.00, 0.00, 0.10, 1.40}
        });
        final OUProcessModel process = new OUProcessModel(
                prefix + ".ou",
                dimension,
                block,
                diffusion,
                mean,
                initialCovariance,
                OUProcessModel.CovarianceGradientMethod.STATIONARY_LYAPUNOV);

        final MatrixParameter observationMatrix = identity(prefix + ".H", dimension);
        final MatrixParameter observationNoise = makeMatrix(prefix + ".Robs", new double[][]{
                {0.28, 0.01, 0.00, 0.00},
                {0.01, 0.31, 0.01, 0.00},
                {0.00, 0.01, 0.34, 0.01},
                {0.00, 0.00, 0.01, 0.29}
        });
        final double y03 = includeMissingObservation ? Double.NaN : 0.30;
        final double y13 = includeMissingObservation ? Double.NaN : -0.22;
        final double y23 = includeMissingObservation ? Double.NaN : 0.18;
        final double y33 = includeMissingObservation ? Double.NaN : 0.07;
        final MatrixParameter observations = makeMatrix(prefix + ".Y", new double[][]{
                {0.12, -0.20, 0.35, y03, 0.41},
                {-0.08, 0.26, -0.18, y13, 0.09},
                {0.21, -0.11, 0.14, y23, -0.16},
                {-0.04, 0.19, -0.09, y33, 0.27}
        });
        final GaussianObservationModel observation = new GaussianObservationModel(
                prefix + ".obs",
                dimension,
                observationMatrix,
                observationNoise,
                observations);
        final TimeGrid grid = new UniformTimeGrid(observations.getColumnDimension(), 0.0, dt);
        final GaussianTransitionRepresentation transitionRepresentation =
                representation(process, GaussianTransitionRepresentation.class);
        final KalmanLikelihoodEngine likelihoodEngine = new KalmanLikelihoodEngine(
                transitionRepresentation,
                observation,
                grid);

        return new GeneralBlockModel(
                process,
                block,
                diffusion,
                mean,
                observation,
                grid,
                transitionRepresentation,
                likelihoodEngine);
    }

    private static MatrixParameter identity(final String name, final int dimension) {
        final double[][] values = new double[dimension][dimension];
        for (int i = 0; i < dimension; ++i) {
            values[i][i] = 1.0;
        }
        return makeMatrix(name, values);
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

    private static long blockFrechetApplicationCount() {
        return BlockDiagonalFrechetHelper.getExactPlanApplicationCount()
                + BlockDiagonalFrechetHelper.getEqualDiagonalPlanApplicationCount();
    }

    private static final class GeneralBlockModel {
        final OUProcessModel process;
        final BlockDiagonalPolarStableMatrixParameter block;
        final MatrixParameter diffusion;
        final Parameter mean;
        final GaussianObservationModel observation;
        final TimeGrid grid;
        final GaussianTransitionRepresentation transitionRepresentation;
        final KalmanLikelihoodEngine likelihoodEngine;

        private GeneralBlockModel(final OUProcessModel process,
                                  final BlockDiagonalPolarStableMatrixParameter block,
                                  final MatrixParameter diffusion,
                                  final Parameter mean,
                                  final GaussianObservationModel observation,
                                  final TimeGrid grid,
                                  final GaussianTransitionRepresentation transitionRepresentation,
                                  final KalmanLikelihoodEngine likelihoodEngine) {
            this.process = process;
            this.block = block;
            this.diffusion = diffusion;
            this.mean = mean;
            this.observation = observation;
            this.grid = grid;
            this.transitionRepresentation = transitionRepresentation;
            this.likelihoodEngine = likelihoodEngine;
        }
    }
}
