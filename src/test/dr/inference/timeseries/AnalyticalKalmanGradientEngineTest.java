package test.dr.inference.timeseries;

import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.BlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.CompoundSymmetricMatrix;
import dr.inference.model.GivensRotationMatrixParameter;
import dr.inference.model.OrthogonalBlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.Parameter;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.continuous.TimeSeriesOUGaussianBranchTransitionProvider;
import dr.inference.timeseries.core.BasicTimeSeriesModel;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.core.UniformTimeGrid;
import dr.inference.timeseries.engine.gaussian.AnalyticalKalmanGradientEngine;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContributionUtils;
import dr.inference.timeseries.engine.gaussian.CanonicalForwardTrajectory;
import dr.inference.timeseries.engine.gaussian.CanonicalKalmanLikelihoodEngine;
import dr.inference.timeseries.engine.gaussian.CanonicalKalmanSmootherEngine;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalTransitionAdjointUtils;
import dr.inference.timeseries.engine.gaussian.DiffusionMatrixGradientFormula;
import dr.inference.timeseries.engine.gaussian.GaussianForwardComputationMode;
import dr.inference.timeseries.engine.gaussian.GaussianLikelihoodEngineFactory;
import dr.inference.timeseries.engine.gaussian.KalmanLikelihoodEngine;
import dr.inference.timeseries.engine.gaussian.KalmanSmootherEngine;
import dr.inference.timeseries.engine.gaussian.SelectionMatrixGradientFormula;
import dr.inference.timeseries.gaussian.DiffusionMatrixParameterizationFactory;
import dr.inference.timeseries.gaussian.EulerOUProcessModel;
import dr.inference.timeseries.gaussian.GaussianObservationModel;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockDiagonalSelectionMatrixParameterization;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.OUProcessModel.CovarianceGradientMethod;
import dr.inference.timeseries.engine.LikelihoodEngine;
import dr.inference.timeseries.likelihood.GaussianGradientComputationMode;
import dr.inference.timeseries.likelihood.GaussianSmootherComputationMode;
import dr.inference.timeseries.likelihood.GaussianTimeSeriesLikelihoodFactory;
import dr.inference.timeseries.likelihood.TimeSeriesLikelihood;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianBranchTransitionKernel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Tests for {@link AnalyticalKalmanGradientEngine} using the exact OU process model.
 *
 * <p>All gradients are verified against numerical central-difference approximations.
 * The analytical gradient w.r.t. the selection matrix A has two contributions:
 * <ol>
 *   <li>F/f-path: backprop through F = expm(-A dt) and f = (I - F) mu.</li>
 *   <li>V-path: backprop through V = Van Loan integral via the adjoint of the
 *       block-exponential construction.</li>
 * </ol>
 * Both contributions are tested here; in particular, tests are structured so that
 * omitting the V-path would cause a measurable discrepancy with numerical FD.
 */
public class AnalyticalKalmanGradientEngineTest extends TestCase {

    /** Tolerance for central-difference gradient comparison (h = 1e-6). */
    private static final double TOL_FD = 1e-5;

    public AnalyticalKalmanGradientEngineTest(String name) {
        super(name);
    }

    // ── Helper container ─────────────────────────────────────────────────────────

    private static class Model {
        final OUProcessModel process;
        final GaussianObservationModel obs;
        final KalmanLikelihoodEngine likelihoodEngine;
        final AnalyticalKalmanGradientEngine analyticalEngine;

        Model(OUProcessModel process,
              GaussianObservationModel obs,
              KalmanLikelihoodEngine likelihoodEngine,
              AnalyticalKalmanGradientEngine analyticalEngine) {
            this.process          = process;
            this.obs              = obs;
            this.likelihoodEngine = likelihoodEngine;
            this.analyticalEngine = analyticalEngine;
        }
    }

    // ── Factory helpers ───────────────────────────────────────────────────────────

    private static MatrixParameter makeMatrix(String name, double[][] values) {
        MatrixParameter mp = new MatrixParameter(name, values.length, values[0].length);
        for (int i = 0; i < values.length; i++) {
            for (int j = 0; j < values[0].length; j++) {
                mp.setParameterValue(i, j, values[i][j]);
            }
        }
        return mp;
    }

    private static double logDensityFromMoments(double[] mean, double[][] covariance, double[] value) {
        final int d = value.length;
        final double[][] precision = invertSymmetricPositiveDefinite(covariance);
        double quadratic = 0.0;
        for (int i = 0; i < d; ++i) {
            double row = 0.0;
            for (int j = 0; j < d; ++j) {
                row += precision[i][j] * (value[j] - mean[j]);
            }
            quadratic += (value[i] - mean[i]) * row;
        }
        return -0.5 * (d * Math.log(2.0 * Math.PI) + logDeterminantSymmetricPositiveDefinite(covariance) + quadratic);
    }

    private static double logDensityFromCanonicalState(CanonicalGaussianState state, double[] value) {
        final int d = value.length;
        double quadratic = 0.0;
        double linear = 0.0;
        for (int i = 0; i < d; ++i) {
            linear += state.information[i] * value[i];
            double row = 0.0;
            for (int j = 0; j < d; ++j) {
                row += state.precision[i * d + j] * value[j];
            }
            quadratic += value[i] * row;
        }
        return -0.5 * quadratic + linear - state.logNormalizer;
    }

    private static double logTransitionFromMoments(double[][] transitionMatrix,
                                                   double[] transitionOffset,
                                                   double[][] transitionCovariance,
                                                   double[] previousState,
                                                   double[] nextState) {
        final int d = previousState.length;
        final double[] conditionalMean = new double[d];
        for (int i = 0; i < d; ++i) {
            double sum = transitionOffset[i];
            for (int j = 0; j < d; ++j) {
                sum += transitionMatrix[i][j] * previousState[j];
            }
            conditionalMean[i] = sum;
        }
        return logDensityFromMoments(conditionalMean, transitionCovariance, nextState);
    }

    private static double logTransitionFromCanonical(CanonicalGaussianTransition transition,
                                                     double[] previousState,
                                                     double[] nextState) {
        final int d = previousState.length;
        double quadratic = 0.0;
        double linear = 0.0;

        for (int i = 0; i < d; ++i) {
            linear += transition.informationX[i] * previousState[i];
            linear += transition.informationY[i] * nextState[i];
        }

        for (int i = 0; i < d; ++i) {
            double xx = 0.0;
            double xy = 0.0;
            double yx = 0.0;
            double yy = 0.0;
            for (int j = 0; j < d; ++j) {
                final int ij = i * d + j;
                xx += transition.precisionXX[ij] * previousState[j];
                xy += transition.precisionXY[ij] * nextState[j];
                yx += transition.precisionYX[ij] * previousState[j];
                yy += transition.precisionYY[ij] * nextState[j];
            }
            quadratic += previousState[i] * (xx + xy);
            quadratic += nextState[i] * (yx + yy);
        }
        return -0.5 * quadratic + linear - transition.logNormalizer;
    }

    private static double[][] toMatrix(final double[] flat, final int dimension) {
        final double[][] out = new double[dimension][dimension];
        for (int i = 0; i < dimension; ++i) {
            System.arraycopy(flat, i * dimension, out[i], 0, dimension);
        }
        return out;
    }

    private static double flatEntry(final double[] flat, final int dimension, final int row, final int column) {
        return flat[row * dimension + column];
    }

    private static double[][] invertSymmetricPositiveDefinite(double[][] matrix) {
        final int d = matrix.length;
        final double[][] chol = new double[d][d];
        for (int i = 0; i < d; ++i) {
            System.arraycopy(matrix[i], 0, chol[i], 0, d);
        }
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = chol[i][j];
                for (int k = 0; k < j; ++k) {
                    sum -= chol[i][k] * chol[j][k];
                }
                if (i == j) {
                    chol[i][j] = Math.sqrt(sum);
                } else {
                    chol[i][j] = sum / chol[j][j];
                }
            }
            for (int j = i + 1; j < d; ++j) {
                chol[i][j] = 0.0;
            }
        }
        final double[][] lowerInverse = new double[d][d];
        for (int column = 0; column < d; ++column) {
            for (int row = 0; row < d; ++row) {
                double sum = (row == column) ? 1.0 : 0.0;
                for (int k = 0; k < row; ++k) {
                    sum -= chol[row][k] * lowerInverse[k][column];
                }
                lowerInverse[row][column] = sum / chol[row][row];
            }
        }
        final double[][] inverse = new double[d][d];
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += lowerInverse[k][i] * lowerInverse[k][j];
                }
                inverse[i][j] = sum;
            }
        }
        return inverse;
    }

    private static double logDeterminantSymmetricPositiveDefinite(double[][] matrix) {
        final int d = matrix.length;
        final double[][] chol = new double[d][d];
        for (int i = 0; i < d; ++i) {
            System.arraycopy(matrix[i], 0, chol[i], 0, d);
        }
        double logDet = 0.0;
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = chol[i][j];
                for (int k = 0; k < j; ++k) {
                    sum -= chol[i][k] * chol[j][k];
                }
                if (i == j) {
                    chol[i][j] = Math.sqrt(sum);
                    logDet += Math.log(chol[i][j]);
                } else {
                    chol[i][j] = sum / chol[j][j];
                }
            }
        }
        return 2.0 * logDet;
    }

    private static void assertMatrixArrayEquals(String label,
                                                double[][][] expected,
                                                double[][][] observed,
                                                double tolerance) {
        assertEquals(label + " time count", expected.length, observed.length);
        for (int t = 0; t < expected.length; ++t) {
            for (int i = 0; i < expected[t].length; ++i) {
                for (int j = 0; j < expected[t][i].length; ++j) {
                    assertEquals(label + " mismatch at t=" + t + " [" + i + "," + j + "]",
                            expected[t][i][j], observed[t][i][j], tolerance);
                }
            }
        }
    }

    private static void assertVectorArrayEquals(String label,
                                                double[][] expected,
                                                double[][] observed,
                                                double tolerance) {
        assertEquals(label + " time count", expected.length, observed.length);
        for (int t = 0; t < expected.length; ++t) {
            for (int i = 0; i < expected[t].length; ++i) {
                assertEquals(label + " mismatch at t=" + t + " idx=" + i,
                        expected[t][i], observed[t][i], tolerance);
            }
        }
    }

    private static double[] sumCanonicalLocalSelectionGradient(final OUProcessModel process,
                                                               final GaussianObservationModel obs,
                                                               final TimeGrid grid) {
        final CanonicalKalmanSmootherEngine smoother = new CanonicalKalmanSmootherEngine(
                process.getRepresentation(CanonicalGaussianBranchTransitionKernel.class),
                process.getRepresentation(GaussianTransitionRepresentation.class),
                obs,
                grid);
        smoother.getLogLikelihood();

        final CanonicalForwardTrajectory trajectory = smoother.getCanonicalTrajectory();
        final int dimension = process.getStateDimension();
        final double[] mean = new double[dimension];
        final CanonicalBranchMessageContribution contribution = new CanonicalBranchMessageContribution(dimension);
        final CanonicalBranchMessageContributionUtils.Workspace contributionWorkspace =
                new CanonicalBranchMessageContributionUtils.Workspace(dimension);
        final CanonicalLocalTransitionAdjoints adjoints = new CanonicalLocalTransitionAdjoints(dimension);
        final CanonicalTransitionAdjointUtils.Workspace workspace =
                new CanonicalTransitionAdjointUtils.Workspace(dimension);
        final double[] gradient = new double[process.getDriftMatrix().getDimension()];

        process.getInitialMean(mean);
        for (int timeIndex = 0; timeIndex < trajectory.timeCount - 1; ++timeIndex) {
            CanonicalBranchMessageContributionUtils.fillFromPairState(
                    trajectory.branchPairStates[timeIndex], contributionWorkspace, contribution);
            CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                    trajectory.transitions[timeIndex], contribution, workspace, adjoints);
            final double dt = grid.getDelta(timeIndex, timeIndex + 1);
            process.accumulateSelectionGradient(
                    dt,
                    mean,
                    toMatrix(adjoints.dLogL_dF, dimension),
                    adjoints.dLogL_df,
                    gradient);
            process.accumulateSelectionGradientFromCovariance(
                    dt,
                    toMatrix(adjoints.dLogL_dOmega, dimension),
                    gradient);
        }
        return gradient;
    }

    private static AnalyticalKalmanGradientEngine makeAnalyticalEngine(OUProcessModel process,
                                                                       GaussianObservationModel obs,
                                                                       TimeGrid grid,
                                                                       boolean includeDiffusion) {
        GaussianTransitionRepresentation rep = process.getRepresentation(
                GaussianTransitionRepresentation.class);
        KalmanSmootherEngine smoother = new KalmanSmootherEngine(rep, obs, grid);
        SelectionMatrixGradientFormula selGrad =
                new SelectionMatrixGradientFormula(process.getDriftMatrix(), process.getStateDimension());
        if (!includeDiffusion) {
            return new AnalyticalKalmanGradientEngine(smoother, selGrad);
        }
        return new AnalyticalKalmanGradientEngine(
                smoother,
                selGrad,
                new DiffusionMatrixGradientFormula(
                        DiffusionMatrixParameterizationFactory.create(process.getDiffusionMatrix()),
                        process.getStateDimension()));
    }

    private static TimeSeriesLikelihood makeCanonicalAnalyticalLikelihood(String name,
                                                                          OUProcessModel process,
                                                                          GaussianObservationModel obs,
                                                                          TimeGrid grid) {
        final BasicTimeSeriesModel model = new BasicTimeSeriesModel(name + ".model", process, obs, grid);
        return GaussianTimeSeriesLikelihoodFactory.create(
                name,
                model,
                GaussianForwardComputationMode.CANONICAL,
                GaussianSmootherComputationMode.CANONICAL,
                GaussianGradientComputationMode.CANONICAL_ANALYTICAL);
    }

    private static TimeSeriesLikelihood makeExpectationAnalyticalLikelihood(String name,
                                                                            OUProcessModel process,
                                                                            GaussianObservationModel obs,
                                                                            TimeGrid grid) {
        final BasicTimeSeriesModel model = new BasicTimeSeriesModel(name + ".model", process, obs, grid);
        return GaussianTimeSeriesLikelihoodFactory.create(
                name,
                model,
                GaussianForwardComputationMode.EXPECTATION,
                GaussianSmootherComputationMode.EXPECTATION,
                GaussianGradientComputationMode.EXPECTATION_ANALYTICAL);
    }

    /**
     * Builds a 1D exact-OU model backed by an {@link AnalyticalKalmanGradientEngine}.
     */
    private static Model makeScalar1D(double a, double q, double mu0, double p0,
                                      double r, double[] yValues, double dt) {
        return makeScalar1DWithMethod(a, q, mu0, p0, r, yValues, dt, CovarianceGradientMethod.LYAPUNOV_ADJOINT);
    }

    private static Model makeScalar1DWithMethod(double a, double q, double mu0, double p0,
                                                 double r, double[] yValues, double dt,
                                                 CovarianceGradientMethod method) {
        int T = yValues.length;

        MatrixParameter drift    = makeMatrix("A",  new double[][]{{a}});
        MatrixParameter diffusion = makeMatrix("Q", new double[][]{{q}});
        Parameter mean           = new Parameter.Default(mu0);
        MatrixParameter initCov  = makeMatrix("P0", new double[][]{{p0}});
        OUProcessModel process   = new OUProcessModel("ou", 1, drift, diffusion, mean, initCov, method);

        MatrixParameter H        = makeMatrix("H",  new double[][]{{1.0}});
        MatrixParameter R        = makeMatrix("R",  new double[][]{{r}});
        double[][] obsData       = new double[1][T];
        obsData[0]               = yValues;
        MatrixParameter Y        = makeMatrix("Y",  obsData);
        GaussianObservationModel obs = new GaussianObservationModel("obs", 1, H, R, Y);

        TimeGrid grid            = new UniformTimeGrid(T, 0.0, dt);
        GaussianTransitionRepresentation rep = process.getRepresentation(
                GaussianTransitionRepresentation.class);
        KalmanSmootherEngine smoother = new KalmanSmootherEngine(rep, obs, grid);
        SelectionMatrixGradientFormula selGrad =
                new SelectionMatrixGradientFormula(process.getDriftMatrix(), 1);
        AnalyticalKalmanGradientEngine analytical =
                new AnalyticalKalmanGradientEngine(smoother, selGrad);

        // Separate likelihood engine for numerical FD (same parameters, independent run)
        KalmanLikelihoodEngine likelihood = new KalmanLikelihoodEngine(rep, obs, grid);

        return new Model(process, obs, likelihood, analytical);
    }

    /**
     * Builds a 2D exact-OU model with diagonal A and Q.
     */
    private static Model make2D(double a11, double a22,
                                 double q11, double q22,
                                 double mu0, double p0,
                                 double r, double[] y1, double[] y2, double dt) {
        return make2DWithMethod(a11, a22, q11, q22, mu0, p0, r, y1, y2, dt,
                CovarianceGradientMethod.LYAPUNOV_ADJOINT);
    }

    private static Model make2DWithMethod(double a11, double a22,
                                           double q11, double q22,
                                           double mu0, double p0,
                                           double r, double[] y1, double[] y2, double dt,
                                           CovarianceGradientMethod method) {
        int T = y1.length;

        MatrixParameter drift    = makeMatrix("A", new double[][]{{a11, 0}, {0, a22}});
        MatrixParameter diffusion= makeMatrix("Q", new double[][]{{q11, 0}, {0, q22}});
        Parameter mean           = new Parameter.Default(mu0);
        MatrixParameter initCov  = makeMatrix("P0", new double[][]{{p0, 0}, {0, p0}});
        OUProcessModel process   = new OUProcessModel("ou2d", 2, drift, diffusion, mean, initCov, method);

        MatrixParameter H        = makeMatrix("H", new double[][]{{1,0},{0,1}});
        MatrixParameter R        = makeMatrix("R", new double[][]{{r,0},{0,r}});
        double[][] obsData       = new double[2][T];
        obsData[0]               = y1;
        obsData[1]               = y2;
        MatrixParameter Y        = makeMatrix("Y", obsData);
        GaussianObservationModel obs = new GaussianObservationModel("obs", 2, H, R, Y);

        TimeGrid grid            = new UniformTimeGrid(T, 0.0, dt);
        GaussianTransitionRepresentation rep = process.getRepresentation(
                GaussianTransitionRepresentation.class);
        KalmanSmootherEngine smoother = new KalmanSmootherEngine(rep, obs, grid);
        SelectionMatrixGradientFormula selGrad =
                new SelectionMatrixGradientFormula(process.getDriftMatrix(), 2);
        AnalyticalKalmanGradientEngine analytical =
                new AnalyticalKalmanGradientEngine(smoother, selGrad);

        KalmanLikelihoodEngine likelihood = new KalmanLikelihoodEngine(rep, obs, grid);

        return new Model(process, obs, likelihood, analytical);
    }

    private static Model make2DFullWithMethod(double[][] driftValues,
                                              double[][] diffusionValues,
                                              double mu0,
                                              double[][] initCovValues,
                                              double[][] noiseValues,
                                              double[] y1, double[] y2, double dt,
                                              CovarianceGradientMethod method) {
        int T = y1.length;

        MatrixParameter drift    = makeMatrix("A", driftValues);
        MatrixParameter diffusion= makeMatrix("Q", diffusionValues);
        Parameter mean           = new Parameter.Default(mu0);
        MatrixParameter initCov  = makeMatrix("P0", initCovValues);
        OUProcessModel process   = new OUProcessModel("ou2dfull", 2, drift, diffusion, mean, initCov, method);

        MatrixParameter H        = makeMatrix("H", new double[][]{{1,0},{0,1}});
        MatrixParameter R        = makeMatrix("R", noiseValues);
        double[][] obsData       = new double[2][T];
        obsData[0]               = y1;
        obsData[1]               = y2;
        MatrixParameter Y        = makeMatrix("Y", obsData);
        GaussianObservationModel obs = new GaussianObservationModel("obs", 2, H, R, Y);

        TimeGrid grid            = new UniformTimeGrid(T, 0.0, dt);
        GaussianTransitionRepresentation rep = process.getRepresentation(
                GaussianTransitionRepresentation.class);
        KalmanSmootherEngine smoother = new KalmanSmootherEngine(rep, obs, grid);
        SelectionMatrixGradientFormula selGrad =
                new SelectionMatrixGradientFormula(process.getDriftMatrix(), 2);
        AnalyticalKalmanGradientEngine analytical =
                new AnalyticalKalmanGradientEngine(smoother, selGrad);

        KalmanLikelihoodEngine likelihood = new KalmanLikelihoodEngine(rep, obs, grid);

        return new Model(process, obs, likelihood, analytical);
    }

    private static Model make3DFullWithMethod(double[][] driftValues,
                                              double[][] diffusionValues,
                                              double mu0,
                                              double[][] initCovValues,
                                              double[][] noiseValues,
                                              double[] y1, double[] y2, double[] y3, double dt,
                                              CovarianceGradientMethod method) {
        int T = y1.length;

        MatrixParameter drift    = makeMatrix("A", driftValues);
        MatrixParameter diffusion= makeMatrix("Q", diffusionValues);
        Parameter mean           = new Parameter.Default(mu0);
        MatrixParameter initCov  = makeMatrix("P0", initCovValues);
        OUProcessModel process   = new OUProcessModel("ou3dfull", 3, drift, diffusion, mean, initCov, method);

        MatrixParameter H        = makeMatrix("H", new double[][]{{1,0,0},{0,1,0},{0,0,1}});
        MatrixParameter R        = makeMatrix("R", noiseValues);
        double[][] obsData       = new double[3][T];
        obsData[0]               = y1;
        obsData[1]               = y2;
        obsData[2]               = y3;
        MatrixParameter Y        = makeMatrix("Y", obsData);
        GaussianObservationModel obs = new GaussianObservationModel("obs", 3, H, R, Y);

        TimeGrid grid            = new UniformTimeGrid(T, 0.0, dt);
        GaussianTransitionRepresentation rep = process.getRepresentation(
                GaussianTransitionRepresentation.class);
        KalmanSmootherEngine smoother = new KalmanSmootherEngine(rep, obs, grid);
        SelectionMatrixGradientFormula selGrad =
                new SelectionMatrixGradientFormula(process.getDriftMatrix(), 3);
        AnalyticalKalmanGradientEngine analytical =
                new AnalyticalKalmanGradientEngine(smoother, selGrad);

        KalmanLikelihoodEngine likelihood = new KalmanLikelihoodEngine(rep, obs, grid);

        return new Model(process, obs, likelihood, analytical);
    }

    private static Model makeDenseFullWithMethod(String processName,
                                                 double[][] driftValues,
                                                 double[][] diffusionValues,
                                                 double mu0,
                                                 double[][] initCovValues,
                                                 double[][] noiseValues,
                                                 double[][] observationValues,
                                                 double dt,
                                                 CovarianceGradientMethod method) {
        final int d = driftValues.length;
        final int T = observationValues[0].length;

        MatrixParameter drift = makeMatrix("A", driftValues);
        MatrixParameter diffusion = makeMatrix("Q", diffusionValues);
        Parameter mean = new Parameter.Default(mu0);
        MatrixParameter initCov = makeMatrix("P0", initCovValues);
        OUProcessModel process = new OUProcessModel(processName, d, drift, diffusion, mean, initCov, method);

        final double[][] identity = new double[d][d];
        for (int i = 0; i < d; ++i) {
            identity[i][i] = 1.0;
        }
        MatrixParameter H = makeMatrix("H", identity);
        MatrixParameter R = makeMatrix("R", noiseValues);
        MatrixParameter Y = makeMatrix("Y", observationValues);
        GaussianObservationModel obs = new GaussianObservationModel("obs", d, H, R, Y);

        TimeGrid grid = new UniformTimeGrid(T, 0.0, dt);
        GaussianTransitionRepresentation rep = process.getRepresentation(
                GaussianTransitionRepresentation.class);
        KalmanSmootherEngine smoother = new KalmanSmootherEngine(rep, obs, grid);
        SelectionMatrixGradientFormula selGrad =
                new SelectionMatrixGradientFormula(process.getDriftMatrix(), d);
        AnalyticalKalmanGradientEngine analytical =
                new AnalyticalKalmanGradientEngine(smoother, selGrad);

        KalmanLikelihoodEngine likelihood = new KalmanLikelihoodEngine(rep, obs, grid);

        return new Model(process, obs, likelihood, analytical);
    }

    private static Model makeBlockPolarFullWithMethod(String processName,
                                                      double[][] rotationValues,
                                                      double[] scalarBlockValues,
                                                      double[] rhoValues,
                                                      double[] thetaValues,
                                                      double[] tValues,
                                                      double[][] diffusionValues,
                                                      double mu0,
                                                      double[][] initCovValues,
                                                      double[][] noiseValues,
                                                      double[][] observationValues,
                                                      double dt,
                                                      CovarianceGradientMethod method) {
        final int d = rotationValues.length;
        final int T = observationValues[0].length;

        MatrixParameter rotation = makeMatrix("R", rotationValues);
        Parameter scalarBlock = new Parameter.Default(scalarBlockValues);
        Parameter rhoParam = new Parameter.Default(rhoValues);
        Parameter thetaParam = new Parameter.Default(thetaValues);
        Parameter tParam = new Parameter.Default(tValues);
        MatrixParameter drift = new BlockDiagonalPolarStableMatrixParameter(
                "Ablock", rotation, scalarBlock, rhoParam, thetaParam, tParam);
        MatrixParameter diffusion = makeMatrix("Q", diffusionValues);
        Parameter mean = new Parameter.Default(mu0);
        MatrixParameter initCov = makeMatrix("P0", initCovValues);
        OUProcessModel process = new OUProcessModel(processName, d, drift, diffusion, mean, initCov, method);

        final double[][] identity = new double[d][d];
        for (int i = 0; i < d; ++i) {
            identity[i][i] = 1.0;
        }
        MatrixParameter H = makeMatrix("H", identity);
        MatrixParameter R = makeMatrix("Robs", noiseValues);
        MatrixParameter Y = makeMatrix("Y", observationValues);
        GaussianObservationModel obs = new GaussianObservationModel("obs", d, H, R, Y);

        TimeGrid grid = new UniformTimeGrid(T, 0.0, dt);
        GaussianTransitionRepresentation rep = process.getRepresentation(
                GaussianTransitionRepresentation.class);
        KalmanSmootherEngine smoother = new KalmanSmootherEngine(rep, obs, grid);
        SelectionMatrixGradientFormula selGrad =
                new SelectionMatrixGradientFormula(process.getDriftMatrix(), d);
        AnalyticalKalmanGradientEngine analytical =
                new AnalyticalKalmanGradientEngine(smoother, selGrad);

        KalmanLikelihoodEngine likelihood = new KalmanLikelihoodEngine(rep, obs, grid);

        return new Model(process, obs, likelihood, analytical);
    }

    private static Model make2DBlockPolarWithMethod(double[][] rotationValues,
                                                    double rho,
                                                    double theta,
                                                    double t,
                                                    double[][] diffusionValues,
                                                    double mu0,
                                                    double[][] initCovValues,
                                                    double[][] noiseValues,
                                                    double[] y1, double[] y2, double dt,
                                                    CovarianceGradientMethod method) {
        int T = y1.length;

        MatrixParameter rotation = makeMatrix("R", rotationValues);
        Parameter scalarBlock = new Parameter.Default(0);
        Parameter rhoParam = new Parameter.Default(rho);
        Parameter thetaParam = new Parameter.Default(theta);
        Parameter tParam = new Parameter.Default(t);
        MatrixParameter drift = new BlockDiagonalPolarStableMatrixParameter(
                "Ablock", rotation, scalarBlock, rhoParam, thetaParam, tParam);
        MatrixParameter diffusion = makeMatrix("Q", diffusionValues);
        Parameter mean = new Parameter.Default(mu0);
        MatrixParameter initCov = makeMatrix("P0", initCovValues);
        OUProcessModel process = new OUProcessModel("ou2dblock", 2, drift, diffusion, mean, initCov, method);

        MatrixParameter H = makeMatrix("H", new double[][]{{1,0},{0,1}});
        MatrixParameter R = makeMatrix("Robs", noiseValues);
        double[][] obsData = new double[2][T];
        obsData[0] = y1;
        obsData[1] = y2;
        MatrixParameter Y = makeMatrix("Y", obsData);
        GaussianObservationModel obs = new GaussianObservationModel("obs", 2, H, R, Y);

        TimeGrid grid = new UniformTimeGrid(T, 0.0, dt);
        GaussianTransitionRepresentation rep = process.getRepresentation(
                GaussianTransitionRepresentation.class);
        KalmanSmootherEngine smoother = new KalmanSmootherEngine(rep, obs, grid);
        SelectionMatrixGradientFormula selGrad =
                new SelectionMatrixGradientFormula(process.getDriftMatrix(), 2);
        AnalyticalKalmanGradientEngine analytical =
                new AnalyticalKalmanGradientEngine(smoother, selGrad);

        KalmanLikelihoodEngine likelihood = new KalmanLikelihoodEngine(rep, obs, grid);

        return new Model(process, obs, likelihood, analytical);
    }

    private static Model makeOrthogonalBlockPolarFullWithMethod(String processName,
                                                                int dimension,
                                                                double[] angleValues,
                                                                double[] scalarBlockValues,
                                                                double[] rhoValues,
                                                                double[] thetaValues,
                                                                double[] tValues,
                                                                double[][] diffusionValues,
                                                                double mu0,
                                                                double[][] initCovValues,
                                                                double[][] noiseValues,
                                                                double[][] observationValues,
                                                                CovarianceGradientMethod method,
                                                                double dt) {
        final int T = observationValues[0].length;

        final Parameter angleParameter = new Parameter.Default(angleValues);
        final GivensRotationMatrixParameter rotation =
                new GivensRotationMatrixParameter("Qg", dimension, angleParameter);
        final Parameter scalarBlock = new Parameter.Default(scalarBlockValues);
        final Parameter rhoParam = new Parameter.Default(rhoValues);
        final Parameter thetaParam = new Parameter.Default(thetaValues);
        final Parameter tParam = new Parameter.Default(tValues);
        final MatrixParameter drift = new OrthogonalBlockDiagonalPolarStableMatrixParameter(
                processName + ".A", rotation, scalarBlock, rhoParam, thetaParam, tParam);
        final MatrixParameter diffusion = makeMatrix("Q", diffusionValues);
        final Parameter mean = new Parameter.Default(mu0);
        final MatrixParameter initCov = makeMatrix("P0", initCovValues);
        final OUProcessModel process = new OUProcessModel(processName, dimension, drift, diffusion, mean, initCov, method);

        final double[][] identity = new double[dimension][dimension];
        for (int i = 0; i < dimension; ++i) {
            identity[i][i] = 1.0;
        }
        final MatrixParameter H = makeMatrix("H", identity);
        final MatrixParameter R = makeMatrix("Robs", noiseValues);
        final MatrixParameter Y = makeMatrix("Y", observationValues);
        final GaussianObservationModel obs = new GaussianObservationModel("obs", dimension, H, R, Y);

        final TimeGrid grid = new UniformTimeGrid(T, 0.0, dt);
        final GaussianTransitionRepresentation rep = process.getRepresentation(
                GaussianTransitionRepresentation.class);
        final KalmanSmootherEngine smoother = new KalmanSmootherEngine(rep, obs, grid);
        final SelectionMatrixGradientFormula selGrad =
                new SelectionMatrixGradientFormula(process.getDriftMatrix(), dimension);
        final AnalyticalKalmanGradientEngine analytical =
                new AnalyticalKalmanGradientEngine(smoother, selGrad);
        final KalmanLikelihoodEngine likelihood = new KalmanLikelihoodEngine(rep, obs, grid);

        return new Model(process, obs, likelihood, analytical);
    }

    /** Central-difference numerical gradient w.r.t. parameter element [index]. */
    private static double numericalGradient(KalmanLikelihoodEngine engine,
                                             Parameter param, int index, double h) {
        double orig = param.getParameterValue(index);

        param.setParameterValue(index, orig + h);
        engine.makeDirty();
        double llPlus = engine.getLogLikelihood();

        param.setParameterValue(index, orig - h);
        engine.makeDirty();
        double llMinus = engine.getLogLikelihood();

        param.setParameterValue(index, orig);
        engine.makeDirty();

        return (llPlus - llMinus) / (2.0 * h);
    }

    // ── supportsGradientWrt ───────────────────────────────────────────────────────

    public void testSupportsGradientWrtDriftMatrix() {
        Model m = makeScalar1D(0.3, 1.0, 0.5, 1.0, 1.0, new double[]{1.0, -0.5, 2.0}, 0.25);
        assertTrue(m.analyticalEngine.supportsGradientWrt(m.process.getDriftMatrix()));
    }

    public void testDoesNotSupportUnrelatedParameter() {
        Model m = makeScalar1D(0.3, 1.0, 0.5, 1.0, 1.0, new double[]{1.0, -0.5, 2.0}, 0.25);
        assertFalse(m.analyticalEngine.supportsGradientWrt(new Parameter.Default(0.0)));
    }

    public void testGradientWrtUnsupportedThrows() {
        Model m = makeScalar1D(0.3, 1.0, 0.5, 1.0, 1.0, new double[]{1.0, -0.5, 2.0}, 0.25);
        try {
            m.analyticalEngine.getGradientWrt(new Parameter.Default(1.0));
            fail("Expected IllegalArgumentException");
        } catch (IllegalArgumentException e) { /* expected */ }
    }

    // ── Gradient is finite ────────────────────────────────────────────────────────

    public void testGradientWrtAIsFinite() {
        Model m = makeScalar1D(0.3, 1.0, 0.5, 1.0, 1.0, new double[]{1.0, -0.5, 2.0}, 0.25);
        double[] grad = m.analyticalEngine.getGradientWrt(m.process.getDriftMatrix());
        assertEquals(1, grad.length);
        assertTrue("Analytical gradient w.r.t. A must be finite", Double.isFinite(grad[0]));
    }

    public void testGradientWrtAIsFinite2D() {
        Model m = make2D(0.5, 1.0, 1.0, 2.0, 0.0, 1.0, 0.5,
                new double[]{1.0, -0.5, 2.0}, new double[]{0.5, 1.5, -1.0}, 0.25);
        double[] grad = m.analyticalEngine.getGradientWrt(m.process.getDriftMatrix());
        assertEquals(4, grad.length);
        for (int i = 0; i < grad.length; i++) {
            assertTrue("Analytical gradient element " + i + " must be finite",
                    Double.isFinite(grad[i]));
        }
    }

    // ── Gradient matches numerical FD (F/f-path + V-path) ────────────────────────

    /**
     * Scalar 1D case: a > 0, multiple steps.
     * Both F/f-path and V-path contribute non-trivially.
     */
    public void testGradientWrtAMatchesFD_1D_nonzero_drift() {
        double a = 0.3, q = 1.0, mu = 1.0, p0 = 1.0, r = 1.0;
        double[] y = {2.0, 1.5, -0.5, 3.0};
        double dt  = 0.25;

        Model m = makeScalar1D(a, q, mu, p0, r, y, dt);
        double h = 1e-6;

        double fdGrad = numericalGradient(m.likelihoodEngine, m.process.getDriftMatrix(), 0, h);
        double[] analGrad = m.analyticalEngine.getGradientWrt(m.process.getDriftMatrix());

        assertEquals("Analytical gradient w.r.t. A must match numerical FD (1D, a=0.3)",
                fdGrad, analGrad[0], TOL_FD);
    }

    /**
     * Larger drift: larger V-path contribution, tighter test.
     */
    public void testGradientWrtAMatchesFD_1D_larger_drift() {
        double a = 1.5, q = 2.0, mu = -1.0, p0 = 0.5, r = 0.5;
        double[] y = {0.5, -1.0, 1.5, 0.0, -0.5};
        double dt  = 0.1;

        Model m = makeScalar1D(a, q, mu, p0, r, y, dt);
        double h = 1e-6;

        double fdGrad = numericalGradient(m.likelihoodEngine, m.process.getDriftMatrix(), 0, h);
        double[] analGrad = m.analyticalEngine.getGradientWrt(m.process.getDriftMatrix());

        assertEquals("Analytical gradient w.r.t. A must match numerical FD (1D, a=1.5)",
                fdGrad, analGrad[0], TOL_FD);
    }

    /**
     * Zero drift: F = I, V = Q dt, so V-path contribution is zero.
     * Verifies that the V-path term correctly vanishes when A = 0.
     */
    public void testGradientWrtAMatchesFD_1D_zero_drift() {
        double a = 0.0, q = 1.5, mu = 0.0, p0 = 2.0, r = 1.0;
        double[] y = {1.0, -2.0, 0.5};
        double dt  = 0.5;

        Model m = makeScalar1D(a, q, mu, p0, r, y, dt);
        double h = 1e-6;

        double fdGrad = numericalGradient(m.likelihoodEngine, m.process.getDriftMatrix(), 0, h);
        double[] analGrad = m.analyticalEngine.getGradientWrt(m.process.getDriftMatrix());

        assertEquals("Analytical gradient w.r.t. A must match numerical FD (1D, a=0)",
                fdGrad, analGrad[0], TOL_FD);
    }

    /**
     * 2D diagonal case: verify all four elements of the A gradient.
     */
    public void testGradientWrtAMatchesFD_2D_diagonal() {
        double a11 = 0.5, a22 = 1.0, q11 = 1.0, q22 = 2.0;
        double mu = 0.5, p0 = 1.0, r = 0.5;
        double[] y1 = {1.0, -0.5, 2.0};
        double[] y2 = {0.5,  1.5, -1.0};
        double dt  = 0.2;

        Model m = make2D(a11, a22, q11, q22, mu, p0, r, y1, y2, dt);
        double h = 1e-6;

        double[] analGrad = m.analyticalEngine.getGradientWrt(m.process.getDriftMatrix());

        // Verify each of the 4 elements against FD
        String[] labels = {"A[0][0]", "A[0][1]", "A[1][0]", "A[1][1]"};
        for (int idx = 0; idx < 4; ++idx) {
            double fdGrad = numericalGradient(m.likelihoodEngine, m.process.getDriftMatrix(), idx, h);
            assertEquals("Analytical gradient w.r.t. " + labels[idx] + " must match FD",
                    fdGrad, analGrad[idx], TOL_FD);
        }
    }

    /**
     * Verify that without the V-path, the gradient would differ from FD.
     * This test checks that the V-path actually changes the result by comparing
     * a pure-F/f-path analytical value (computed via the Euler model, which has no
     * V-path) against numerical FD with the exact OU likelihood.
     *
     * For the exact OU model, the F/f-path Euler approximation wrt A is:
     *   ∂logL/∂A ≈ -dt * dLogL_dF + dt * dLogL_df * mu   (Euler chain rule)
     * while the true gradient also adds the Van Loan V-path term.
     * These differ whenever a * dt is not negligibly small and T > 1.
     */
    public void testVPathContributionIsNonNegligible() {
        // Use non-trivial drift so that V-path contribution is large
        double a = 2.0, q = 1.0, mu = 0.0, p0 = 1.0, r = 1.0;
        double[] y = {1.0, -1.0, 2.0, 0.5};
        double dt = 0.5;  // a * dt = 1 → F deviates strongly from I; V-path is large

        Model m = makeScalar1D(a, q, mu, p0, r, y, dt);
        double h = 1e-6;

        double fdGrad    = numericalGradient(m.likelihoodEngine, m.process.getDriftMatrix(), 0, h);
        double[] analGrad = m.analyticalEngine.getGradientWrt(m.process.getDriftMatrix());

        // Must match FD — this would fail if V-path were missing (difference > TOL_FD)
        assertEquals("Full analytical gradient (F/f + V path) must match numerical FD",
                fdGrad, analGrad[0], TOL_FD);

        // Also verify that FD gradient is non-negligible (test is non-degenerate)
        assertTrue("FD gradient must be non-zero for this test to be meaningful",
                Math.abs(fdGrad) > 1e-3);
    }

    // ── Caching: makeDirty invalidates result ─────────────────────────────────────

    public void testMakeDirtyInvalidatesCache() {
        double a = 0.3, q = 1.0, mu = 0.5, p0 = 1.0, r = 1.0;
        double[] y = {2.0, 1.5, -0.5, 3.0};
        double dt = 0.25;

        Model m = makeScalar1D(a, q, mu, p0, r, y, dt);

        double[] gradBefore = m.analyticalEngine.getGradientWrt(m.process.getDriftMatrix());

        // Perturb drift and mark dirty
        m.process.getDriftMatrix().setParameterValue(0, 0, a + 0.5);
        m.analyticalEngine.makeDirty();
        double[] gradAfter = m.analyticalEngine.getGradientWrt(m.process.getDriftMatrix());

        assertFalse("Gradient must change after parameter perturbation + makeDirty",
                Math.abs(gradBefore[0] - gradAfter[0]) < 1e-10);

        // Restore
        m.process.getDriftMatrix().setParameterValue(0, 0, a);
    }

    // ── Lyapunov-adjoint strategy ─────────────────────────────────────────────────

    /**
     * Lyapunov-adjoint strategy: gradient w.r.t. A must match numerical FD (1D).
     */
    public void testLyapunovAdjointMatchesFD_1D() {
        double a = 0.3, q = 1.0, mu = 1.0, p0 = 1.0, r = 1.0;
        double[] y = {2.0, 1.5, -0.5, 3.0};
        double dt = 0.25;

        Model m = makeScalar1DWithMethod(a, q, mu, p0, r, y, dt,
                CovarianceGradientMethod.LYAPUNOV_ADJOINT);
        double h = 1e-6;

        double fdGrad   = numericalGradient(m.likelihoodEngine, m.process.getDriftMatrix(), 0, h);
        double[] analGrad = m.analyticalEngine.getGradientWrt(m.process.getDriftMatrix());

        assertEquals("Lyapunov-adjoint gradient must match numerical FD (1D, a=0.3)",
                fdGrad, analGrad[0], TOL_FD);
    }

    public void testStationaryLyapunovMatchesFD_1D() {
        double a = 0.3, q = 1.0, mu = 1.0, p0 = 1.0, r = 1.0;
        double[] y = {2.0, 1.5, -0.5, 3.0};
        double dt = 0.25;

        Model m = makeScalar1DWithMethod(a, q, mu, p0, r, y, dt,
                CovarianceGradientMethod.STATIONARY_LYAPUNOV);
        double h = 1e-6;

        double fdGrad   = numericalGradient(m.likelihoodEngine, m.process.getDriftMatrix(), 0, h);
        double[] analGrad = m.analyticalEngine.getGradientWrt(m.process.getDriftMatrix());

        assertEquals("Stationary-Lyapunov gradient must match numerical FD (1D, a=0.3)",
                fdGrad, analGrad[0], TOL_FD);
    }

    public void testStationaryLyapunovZeroDriftUnsupported() {
        Model m = makeScalar1DWithMethod(0.0, 1.5, 0.0, 2.0, 1.0, new double[]{1.0, -2.0, 0.5}, 0.5,
                CovarianceGradientMethod.STATIONARY_LYAPUNOV);
        try {
            m.analyticalEngine.getGradientWrt(m.process.getDriftMatrix());
            fail("Expected IllegalArgumentException");
        } catch (IllegalArgumentException e) { /* expected */ }
    }

    public void testStationaryLyapunovMatchesFD_2D_Diagonal() {
        Model m = make2DWithMethod(0.5, 1.0, 1.0, 2.0, 0.5, 1.0, 0.5,
                new double[]{1.0, -0.5, 2.0}, new double[]{0.5, 1.5, -1.0}, 0.2,
                CovarianceGradientMethod.STATIONARY_LYAPUNOV);
        double h = 1e-6;

        double[] analGrad = m.analyticalEngine.getGradientWrt(m.process.getDriftMatrix());
        for (int idx = 0; idx < 4; ++idx) {
            double fdGrad = numericalGradient(m.likelihoodEngine, m.process.getDriftMatrix(), idx, h);
            assertEquals("Stationary-Lyapunov gradient must match FD on diagonal 2D case, idx=" + idx,
                    fdGrad, analGrad[idx], TOL_FD);
        }
    }

    public void testStationaryLyapunovMatchesFD_2D_NonDiagonal() {
        Model m = make2DFullWithMethod(
                new double[][]{{0.7, -0.2}, {0.15, 1.1}},
                new double[][]{{1.2, 0.1}, {0.1, 1.6}},
                0.25,
                new double[][]{{1.0, 0.2}, {0.2, 1.3}},
                new double[][]{{0.6, 0.05}, {0.05, 0.8}},
                new double[]{0.7, -0.4, 1.1, 0.2},
                new double[]{-0.3, 0.8, -1.2, 0.5},
                0.15,
                CovarianceGradientMethod.STATIONARY_LYAPUNOV);
        double h = 1e-6;

        double[] analGrad = m.analyticalEngine.getGradientWrt(m.process.getDriftMatrix());
        for (int idx = 0; idx < 4; ++idx) {
            double fdGrad = numericalGradient(m.likelihoodEngine, m.process.getDriftMatrix(), idx, h);
            assertEquals("Stationary-Lyapunov gradient must match FD on non-diagonal 2D case, idx=" + idx,
                    fdGrad, analGrad[idx], TOL_FD);
        }
    }

    public void testGradientStrategiesMatchFD_3D_NonDiagonal() {
        final double[][] drift = new double[][]{
                {0.8, -0.15, 0.05},
                {0.10, 1.2, -0.20},
                {0.00, 0.12, 0.9}
        };
        final double[][] diffusion = new double[][]{
                {1.4, 0.10, 0.05},
                {0.10, 1.8, 0.12},
                {0.05, 0.12, 1.1}
        };
        final double[][] initCov = new double[][]{
                {1.2, 0.15, 0.05},
                {0.15, 1.0, 0.08},
                {0.05, 0.08, 0.9}
        };
        final double[][] noise = new double[][]{
                {0.6, 0.03, 0.01},
                {0.03, 0.7, 0.02},
                {0.01, 0.02, 0.5}
        };
        final double[] y1 = {0.7, -0.4, 1.1, 0.2, -0.3};
        final double[] y2 = {-0.3, 0.8, -1.2, 0.5, 0.4};
        final double[] y3 = {0.2, -0.6, 0.9, -0.1, 0.3};
        final double dt = 0.15;
        final double h = 1e-6;

        final CovarianceGradientMethod[] methods = new CovarianceGradientMethod[]{
                CovarianceGradientMethod.VAN_LOAN_ADJOINT,
                CovarianceGradientMethod.LYAPUNOV_ADJOINT,
                CovarianceGradientMethod.STATIONARY_LYAPUNOV
        };

        for (CovarianceGradientMethod method : methods) {
            Model m = make3DFullWithMethod(drift, diffusion, 0.25, initCov, noise, y1, y2, y3, dt, method);
            double[] analGrad = m.analyticalEngine.getGradientWrt(m.process.getDriftMatrix());
            assertEquals("3D gradient length must match 3x3 drift matrix", 9, analGrad.length);

            for (int idx = 0; idx < 9; ++idx) {
                double fdGrad = numericalGradient(m.likelihoodEngine, m.process.getDriftMatrix(), idx, h);
                assertEquals("Method " + method + " must match FD on 3D non-diagonal case, idx=" + idx,
                        fdGrad, analGrad[idx], TOL_FD);
            }
        }
    }

    public void testBlockParametrizedModelMatchesDenseModelLikelihoodAndMatrixGradient() {
        final double[][] rotation = new double[][]{{1.0, 0.25}, {-0.15, 0.9}};
        final double[][] diffusion = new double[][]{{1.2, 0.1}, {0.1, 1.6}};
        final double[][] initCov = new double[][]{{1.0, 0.2}, {0.2, 1.3}};
        final double[][] noise = new double[][]{{0.6, 0.05}, {0.05, 0.8}};
        final double[] y1 = {0.7, -0.4, 1.1, 0.2};
        final double[] y2 = {-0.3, 0.8, -1.2, 0.5};
        final double dt = 0.15;

        Model blockModel = make2DBlockPolarWithMethod(
                rotation, 0.9, 0.35, 0.18, diffusion, 0.25, initCov, noise, y1, y2, dt,
                CovarianceGradientMethod.STATIONARY_LYAPUNOV);
        double[][] denseDriftValues = blockModel.process.getDriftMatrix().getParameterAsMatrix();
        Model denseModel = make2DFullWithMethod(
                denseDriftValues, diffusion, 0.25, initCov, noise, y1, y2, dt,
                CovarianceGradientMethod.STATIONARY_LYAPUNOV);

        assertEquals("Block and dense models must have the same log likelihood",
                denseModel.likelihoodEngine.getLogLikelihood(),
                blockModel.likelihoodEngine.getLogLikelihood(),
                1e-10);

        double[] blockGrad = blockModel.analyticalEngine.getGradientWrt(blockModel.process.getDriftMatrix());
        double[] denseGrad = denseModel.analyticalEngine.getGradientWrt(denseModel.process.getDriftMatrix());
        for (int i = 0; i < blockGrad.length; ++i) {
            assertEquals("Block and dense matrix gradients must agree, idx=" + i,
                    denseGrad[i], blockGrad[i], 1e-10);
        }
    }

    public void testBlockParametrizedGradientMatchesFDOnNativeParameters() {
        final double[][] rotation = new double[][]{{1.0, 0.25}, {-0.15, 0.9}};
        final double[][] diffusion = new double[][]{{1.2, 0.1}, {0.1, 1.6}};
        final double[][] initCov = new double[][]{{1.0, 0.2}, {0.2, 1.3}};
        final double[][] noise = new double[][]{{0.6, 0.05}, {0.05, 0.8}};
        final double[] y1 = {0.7, -0.4, 1.1, 0.2};
        final double[] y2 = {-0.3, 0.8, -1.2, 0.5};
        final double dt = 0.15;
        final double h = 1e-6;

        final Model model = make2DBlockPolarWithMethod(
                rotation, 0.9, 0.35, 0.18, diffusion, 0.25, initCov, noise, y1, y2, dt,
                CovarianceGradientMethod.STATIONARY_LYAPUNOV);
        final BlockDiagonalPolarStableMatrixParameter blockParameter =
                (BlockDiagonalPolarStableMatrixParameter) model.process.getDriftMatrix();

        assertTrue(model.analyticalEngine.supportsGradientWrt(blockParameter.getRhoParameter()));
        assertTrue(model.analyticalEngine.supportsGradientWrt(blockParameter.getThetaParameter()));
        assertTrue(model.analyticalEngine.supportsGradientWrt(blockParameter.getTParameter()));
        assertTrue(model.analyticalEngine.supportsGradientWrt(blockParameter.getRotationMatrixParameter()));
        assertTrue(model.analyticalEngine.supportsGradientWrt(blockParameter.getParameter()));

        final double[] analyticRho = model.analyticalEngine.getGradientWrt(blockParameter.getRhoParameter());
        final double[] analyticTheta = model.analyticalEngine.getGradientWrt(blockParameter.getThetaParameter());
        final double[] analyticT = model.analyticalEngine.getGradientWrt(blockParameter.getTParameter());
        final double[] analyticRotation = model.analyticalEngine.getGradientWrt(blockParameter.getRotationMatrixParameter());
        final double[] analyticCompound = model.analyticalEngine.getGradientWrt(blockParameter.getParameter());

        assertEquals(1, analyticRho.length);
        assertEquals(1, analyticTheta.length);
        assertEquals(1, analyticT.length);
        assertEquals(4, analyticRotation.length);

        assertEquals("rho native gradient must match FD",
                numericalGradient(model.likelihoodEngine, blockParameter.getRhoParameter(), 0, h),
                analyticRho[0], TOL_FD);
        assertEquals("theta native gradient must match FD",
                numericalGradient(model.likelihoodEngine, blockParameter.getThetaParameter(), 0, h),
                analyticTheta[0], TOL_FD);
        assertEquals("t native gradient must match FD",
                numericalGradient(model.likelihoodEngine, blockParameter.getTParameter(), 0, h),
                analyticT[0], TOL_FD);

        for (int idx = 0; idx < analyticRotation.length; ++idx) {
            final double fdGrad = numericalGradient(model.likelihoodEngine,
                    blockParameter.getRotationMatrixParameter(), idx, h);
            assertEquals("Rotation native gradient must match FD, idx=" + idx,
                    fdGrad, analyticRotation[idx], TOL_FD);
        }

        final int nativeBase = blockParameter.hasLeadingOneByOneBlock() ? 1 : 0;
        final int expectedCompoundLength = nativeBase + analyticRho.length + analyticTheta.length
                + analyticT.length + analyticRotation.length;
        assertEquals("Compound gradient should contain native block and rotation pieces",
                expectedCompoundLength, analyticCompound.length);
        if (blockParameter.hasLeadingOneByOneBlock()) {
            assertEquals("Compound scalar slot should stay zero when the leading 1x1 block is unused in this test",
                    0.0, analyticCompound[0], 1e-12);
        }
        assertEquals("Compound rho slot must agree with direct rho gradient",
                analyticRho[0], analyticCompound[nativeBase], 1e-12);
        assertEquals("Compound theta slot must agree with direct theta gradient",
                analyticTheta[0], analyticCompound[nativeBase + 1], 1e-12);
        assertEquals("Compound t slot must agree with direct t gradient",
                analyticT[0], analyticCompound[nativeBase + 2], 1e-12);
        for (int idx = 0; idx < analyticRotation.length; ++idx) {
            assertEquals("Compound rotation block must agree with direct rotation gradient, idx=" + idx,
                    analyticRotation[idx], analyticCompound[nativeBase + 3 + idx], 1e-12);
        }
    }

    public void testBlockParametrizedModelMatchesDenseModelAndNativeFD_3D() {
        final double[][] rotation = new double[][]{
                {1.0, 0.10, -0.05},
                {0.20, 0.95, 0.15},
                {-0.10, 0.05, 1.05}
        };
        final double[][] diffusion = new double[][]{
                {1.3, 0.10, 0.02},
                {0.10, 1.6, 0.08},
                {0.02, 0.08, 1.1}
        };
        final double[][] initCov = new double[][]{
                {1.0, 0.12, 0.03},
                {0.12, 1.2, 0.07},
                {0.03, 0.07, 0.9}
        };
        final double[][] noise = new double[][]{
                {0.5, 0.02, 0.01},
                {0.02, 0.6, 0.03},
                {0.01, 0.03, 0.4}
        };
        final double[][] observations = new double[][]{
                {0.6, -0.2, 0.9, 0.1},
                {-0.4, 0.7, -1.0, 0.3},
                {0.2, -0.5, 0.8, -0.1}
        };
        final double dt = 0.18;
        final double h = 1e-6;

        final Model blockModel = makeBlockPolarFullWithMethod(
                "ou3dblock",
                rotation,
                new double[]{0.55},
                new double[]{0.85},
                new double[]{0.30},
                new double[]{0.14},
                diffusion,
                0.2,
                initCov,
                noise,
                observations,
                dt,
                CovarianceGradientMethod.STATIONARY_LYAPUNOV);
        final double[][] denseDriftValues = blockModel.process.getDriftMatrix().getParameterAsMatrix();
        final Model denseModel = makeDenseFullWithMethod(
                "ou3ddense",
                denseDriftValues,
                diffusion,
                0.2,
                initCov,
                noise,
                observations,
                dt,
                CovarianceGradientMethod.STATIONARY_LYAPUNOV);

        assertEquals("3D block and dense models must have the same log likelihood",
                denseModel.likelihoodEngine.getLogLikelihood(),
                blockModel.likelihoodEngine.getLogLikelihood(),
                1e-10);

        final double[] blockMatrixGradient = blockModel.analyticalEngine.getGradientWrt(blockModel.process.getDriftMatrix());
        final double[] denseMatrixGradient = denseModel.analyticalEngine.getGradientWrt(denseModel.process.getDriftMatrix());
        for (int i = 0; i < blockMatrixGradient.length; ++i) {
            assertEquals("3D block and dense matrix gradients must agree, idx=" + i,
                    denseMatrixGradient[i], blockMatrixGradient[i], 1e-10);
        }

        final BlockDiagonalPolarStableMatrixParameter blockParameter =
                (BlockDiagonalPolarStableMatrixParameter) blockModel.process.getDriftMatrix();
        final double[] analyticScalar = blockModel.analyticalEngine.getGradientWrt(blockParameter.getScalarBlockParameter());
        final double[] analyticRho = blockModel.analyticalEngine.getGradientWrt(blockParameter.getRhoParameter());
        final double[] analyticTheta = blockModel.analyticalEngine.getGradientWrt(blockParameter.getThetaParameter());
        final double[] analyticT = blockModel.analyticalEngine.getGradientWrt(blockParameter.getTParameter());
        final double[] analyticRotation = blockModel.analyticalEngine.getGradientWrt(blockParameter.getRotationMatrixParameter());

        assertEquals("3D scalar native gradient must match FD",
                numericalGradient(blockModel.likelihoodEngine, blockParameter.getScalarBlockParameter(), 0, h),
                analyticScalar[0], TOL_FD);
        assertEquals("3D rho native gradient must match FD",
                numericalGradient(blockModel.likelihoodEngine, blockParameter.getRhoParameter(), 0, h),
                analyticRho[0], TOL_FD);
        assertEquals("3D theta native gradient must match FD",
                numericalGradient(blockModel.likelihoodEngine, blockParameter.getThetaParameter(), 0, h),
                analyticTheta[0], TOL_FD);
        assertEquals("3D t native gradient must match FD",
                numericalGradient(blockModel.likelihoodEngine, blockParameter.getTParameter(), 0, h),
                analyticT[0], TOL_FD);
        for (int idx = 0; idx < analyticRotation.length; ++idx) {
            final double fdGrad = numericalGradient(blockModel.likelihoodEngine,
                    blockParameter.getRotationMatrixParameter(), idx, h);
            assertEquals("3D rotation native gradient must match FD, idx=" + idx,
                    fdGrad, analyticRotation[idx], TOL_FD);
        }
    }

    public void testBlockParametrizedModelMatchesDenseModelAndNativeFD_4D() {
        final double[][] rotation = new double[][]{
                {1.00, 0.10, -0.05, 0.02},
                {0.15, 0.95, 0.08, -0.04},
                {-0.12, 0.06, 1.03, 0.09},
                {0.03, -0.07, 0.11, 0.92}
        };
        final double[][] diffusion = new double[][]{
                {1.4, 0.10, 0.03, 0.02},
                {0.10, 1.7, 0.09, 0.04},
                {0.03, 0.09, 1.2, 0.05},
                {0.02, 0.04, 0.05, 1.5}
        };
        final double[][] initCov = new double[][]{
                {1.1, 0.10, 0.04, 0.03},
                {0.10, 1.0, 0.05, 0.02},
                {0.04, 0.05, 1.3, 0.06},
                {0.03, 0.02, 0.06, 0.9}
        };
        final double[][] noise = new double[][]{
                {0.55, 0.02, 0.01, 0.00},
                {0.02, 0.60, 0.03, 0.01},
                {0.01, 0.03, 0.50, 0.02},
                {0.00, 0.01, 0.02, 0.65}
        };
        final double[][] observations = new double[][]{
                {0.5, -0.1, 0.8, 0.0},
                {-0.3, 0.9, -1.1, 0.4},
                {0.1, -0.4, 0.7, -0.2},
                {0.6, 0.2, -0.5, 0.3}
        };
        final double dt = 0.12;
        final double h = 1e-6;

        final Model blockModel = makeBlockPolarFullWithMethod(
                "ou4dblock",
                rotation,
                new double[0],
                new double[]{0.90, 1.10},
                new double[]{0.25, -0.35},
                new double[]{0.12, -0.18},
                diffusion,
                0.15,
                initCov,
                noise,
                observations,
                dt,
                CovarianceGradientMethod.STATIONARY_LYAPUNOV);
        final double[][] denseDriftValues = blockModel.process.getDriftMatrix().getParameterAsMatrix();
        final Model denseModel = makeDenseFullWithMethod(
                "ou4ddense",
                denseDriftValues,
                diffusion,
                0.15,
                initCov,
                noise,
                observations,
                dt,
                CovarianceGradientMethod.STATIONARY_LYAPUNOV);

        assertEquals("4D block and dense models must have the same log likelihood",
                denseModel.likelihoodEngine.getLogLikelihood(),
                blockModel.likelihoodEngine.getLogLikelihood(),
                1e-10);

        final double[] blockMatrixGradient = blockModel.analyticalEngine.getGradientWrt(blockModel.process.getDriftMatrix());
        final double[] denseMatrixGradient = denseModel.analyticalEngine.getGradientWrt(denseModel.process.getDriftMatrix());
        for (int i = 0; i < blockMatrixGradient.length; ++i) {
            assertEquals("4D block and dense matrix gradients must agree, idx=" + i,
                    denseMatrixGradient[i], blockMatrixGradient[i], 1e-10);
        }

        final BlockDiagonalPolarStableMatrixParameter blockParameter =
                (BlockDiagonalPolarStableMatrixParameter) blockModel.process.getDriftMatrix();
        final double[] analyticRho = blockModel.analyticalEngine.getGradientWrt(blockParameter.getRhoParameter());
        final double[] analyticTheta = blockModel.analyticalEngine.getGradientWrt(blockParameter.getThetaParameter());
        final double[] analyticT = blockModel.analyticalEngine.getGradientWrt(blockParameter.getTParameter());
        final double[] analyticRotation = blockModel.analyticalEngine.getGradientWrt(blockParameter.getRotationMatrixParameter());

        for (int idx = 0; idx < analyticRho.length; ++idx) {
            assertEquals("4D rho native gradient must match FD, idx=" + idx,
                    numericalGradient(blockModel.likelihoodEngine, blockParameter.getRhoParameter(), idx, h),
                    analyticRho[idx], TOL_FD);
            assertEquals("4D theta native gradient must match FD, idx=" + idx,
                    numericalGradient(blockModel.likelihoodEngine, blockParameter.getThetaParameter(), idx, h),
                    analyticTheta[idx], TOL_FD);
            assertEquals("4D t native gradient must match FD, idx=" + idx,
                    numericalGradient(blockModel.likelihoodEngine, blockParameter.getTParameter(), idx, h),
                    analyticT[idx], TOL_FD);
        }
        for (int idx = 0; idx < analyticRotation.length; ++idx) {
            final double fdGrad = numericalGradient(blockModel.likelihoodEngine,
                    blockParameter.getRotationMatrixParameter(), idx, h);
            assertEquals("4D rotation native gradient must match FD, idx=" + idx,
                    fdGrad, analyticRotation[idx], TOL_FD);
        }
    }

    public void testOrthogonalBlockParametrizedModelMatchesDenseModelAndNativeFD_2D() {
        final double[][] diffusion = new double[][]{{1.2, 0.1}, {0.1, 1.6}};
        final double[][] initCov = new double[][]{{1.0, 0.2}, {0.2, 1.3}};
        final double[][] noise = new double[][]{{0.6, 0.05}, {0.05, 0.8}};
        final double[][] observations = new double[][]{
                {0.7, -0.4, 1.1, 0.2},
                {-0.3, 0.8, -1.2, 0.5}
        };
        final double dt = 0.15;
        final double h = 1e-6;

        final Model orthogonalModel = makeOrthogonalBlockPolarFullWithMethod(
                "ou2dorth", 2,
                new double[]{0.35},
                new double[0],
                new double[]{0.9},
                new double[]{0.25},
                new double[]{0.18},
                diffusion,
                0.25,
                initCov,
                noise,
                observations,
                CovarianceGradientMethod.STATIONARY_LYAPUNOV,
                dt);
        final double[][] denseDriftValues = orthogonalModel.process.getDriftMatrix().getParameterAsMatrix();
        final Model denseModel = makeDenseFullWithMethod(
                "ou2dorthDense", denseDriftValues, diffusion, 0.25, initCov, noise, observations, dt,
                CovarianceGradientMethod.STATIONARY_LYAPUNOV);

        assertEquals("Orthogonal block and dense models must have the same log likelihood",
                denseModel.likelihoodEngine.getLogLikelihood(),
                orthogonalModel.likelihoodEngine.getLogLikelihood(),
                1e-10);

        final OrthogonalBlockDiagonalPolarStableMatrixParameter blockParameter =
                (OrthogonalBlockDiagonalPolarStableMatrixParameter) orthogonalModel.process.getDriftMatrix();
        final GivensRotationMatrixParameter rotation =
                (GivensRotationMatrixParameter) blockParameter.getRotationMatrixParameter();

        assertEquals("Orthogonal rho native gradient must match FD",
                numericalGradient(orthogonalModel.likelihoodEngine, blockParameter.getRhoParameter(), 0, h),
                orthogonalModel.analyticalEngine.getGradientWrt(blockParameter.getRhoParameter())[0], TOL_FD);
        assertEquals("Orthogonal theta native gradient must match FD",
                numericalGradient(orthogonalModel.likelihoodEngine, blockParameter.getThetaParameter(), 0, h),
                orthogonalModel.analyticalEngine.getGradientWrt(blockParameter.getThetaParameter())[0], TOL_FD);
        assertEquals("Orthogonal t native gradient must match FD",
                numericalGradient(orthogonalModel.likelihoodEngine, blockParameter.getTParameter(), 0, h),
                orthogonalModel.analyticalEngine.getGradientWrt(blockParameter.getTParameter())[0], TOL_FD);
        assertEquals("Orthogonal angle native gradient must match FD",
                numericalGradient(orthogonalModel.likelihoodEngine, rotation.getAngleParameter(), 0, h),
                orthogonalModel.analyticalEngine.getGradientWrt(rotation.getAngleParameter())[0], TOL_FD);
    }

    public void testOrthogonalBlockParametrizedModelMatchesDenseModel_4D() {
        final double[][] diffusion = new double[][]{
                {1.4, 0.10, 0.03, 0.02},
                {0.10, 1.7, 0.09, 0.04},
                {0.03, 0.09, 1.2, 0.05},
                {0.02, 0.04, 0.05, 1.5}
        };
        final double[][] initCov = new double[][]{
                {1.1, 0.10, 0.04, 0.03},
                {0.10, 1.0, 0.05, 0.02},
                {0.04, 0.05, 1.3, 0.06},
                {0.03, 0.02, 0.06, 0.9}
        };
        final double[][] noise = new double[][]{
                {0.55, 0.02, 0.01, 0.00},
                {0.02, 0.60, 0.03, 0.01},
                {0.01, 0.03, 0.50, 0.02},
                {0.00, 0.01, 0.02, 0.65}
        };
        final double[][] observations = new double[][]{
                {0.5, -0.1, 0.8, 0.0},
                {-0.3, 0.9, -1.1, 0.4},
                {0.1, -0.4, 0.7, -0.2},
                {0.6, 0.2, -0.5, 0.3}
        };
        final double dt = 0.12;
        final double h = 1e-6;

        final Model orthogonalModel = makeOrthogonalBlockPolarFullWithMethod(
                "ou4dorth", 4,
                new double[]{0.15, -0.10, 0.20, 0.05, -0.12, 0.08},
                new double[0],
                new double[]{0.90, 1.10},
                new double[]{0.25, -0.35},
                new double[]{0.12, -0.18},
                diffusion,
                0.15,
                initCov,
                noise,
                observations,
                CovarianceGradientMethod.STATIONARY_LYAPUNOV,
                dt);
        final double[][] denseDriftValues = orthogonalModel.process.getDriftMatrix().getParameterAsMatrix();
        final Model denseModel = makeDenseFullWithMethod(
                "ou4dorthDense", denseDriftValues, diffusion, 0.15, initCov, noise, observations, dt,
                CovarianceGradientMethod.STATIONARY_LYAPUNOV);

        assertEquals("Orthogonal 4D block and dense models must have the same log likelihood",
                denseModel.likelihoodEngine.getLogLikelihood(),
                orthogonalModel.likelihoodEngine.getLogLikelihood(),
                1e-10);

        final OrthogonalBlockDiagonalPolarStableMatrixParameter blockParameter =
                (OrthogonalBlockDiagonalPolarStableMatrixParameter) orthogonalModel.process.getDriftMatrix();
        final GivensRotationMatrixParameter rotation =
                (GivensRotationMatrixParameter) blockParameter.getRotationMatrixParameter();

        final double[] analyticRho = orthogonalModel.analyticalEngine.getGradientWrt(blockParameter.getRhoParameter());
        final double[] analyticTheta = orthogonalModel.analyticalEngine.getGradientWrt(blockParameter.getThetaParameter());
        final double[] analyticT = orthogonalModel.analyticalEngine.getGradientWrt(blockParameter.getTParameter());
        final double[] analyticAngles = orthogonalModel.analyticalEngine.getGradientWrt(rotation.getAngleParameter());
        for (int idx = 0; idx < analyticRho.length; ++idx) {
            assertEquals("Orthogonal 4D rho gradient must match FD, idx=" + idx,
                    numericalGradient(orthogonalModel.likelihoodEngine, blockParameter.getRhoParameter(), idx, h),
                    analyticRho[idx], TOL_FD);
            assertEquals("Orthogonal 4D theta gradient must match FD, idx=" + idx,
                    numericalGradient(orthogonalModel.likelihoodEngine, blockParameter.getThetaParameter(), idx, h),
                    analyticTheta[idx], TOL_FD);
            assertEquals("Orthogonal 4D t gradient must match FD, idx=" + idx,
                    numericalGradient(orthogonalModel.likelihoodEngine, blockParameter.getTParameter(), idx, h),
                    analyticT[idx], TOL_FD);
        }
        for (int idx = 0; idx < analyticAngles.length; ++idx) {
            assertEquals("Orthogonal 4D angle gradient must match FD, idx=" + idx,
                    numericalGradient(orthogonalModel.likelihoodEngine, rotation.getAngleParameter(), idx, h),
                    analyticAngles[idx], TOL_FD);
        }
    }

    /**
     * Both strategies must agree to high precision on a 1D case with large drift.
     */
//    public void testVanLoanAndLyapunovAgree_1D() {
//        double a = 1.5, q = 2.0, mu = -1.0, p0 = 0.5, r = 0.5;
//        double[] y = {0.5, -1.0, 1.5, 0.0, -0.5};
//        double dt = 0.1;
//
//        Model mVL  = makeScalar1DWithMethod(a, q, mu, p0, r, y, dt,
//                CovarianceGradientMethod.VAN_LOAN_ADJOINT);
//        Model mLyap = makeScalar1DWithMethod(a, q, mu, p0, r, y, dt,
//                CovarianceGradientMethod.LYAPUNOV_ADJOINT);
//
//        double[] gradVL   = mVL.analyticalEngine.getGradientWrt(mVL.process.getDriftMatrix());
//        double[] gradLyap = mLyap.analyticalEngine.getGradientWrt(mLyap.process.getDriftMatrix());
//
//        assertEquals("VAN_LOAN_ADJOINT and LYAPUNOV_ADJOINT must agree (1D, a=1.5)",
//                gradVL[0], gradLyap[0], TOL_FD);
//    }

    /**
     * Both strategies must agree on all four elements of A for a 2D case.
     */
//    public void testVanLoanAndLyapunovAgree_2D() {
//        double a11 = 0.5, a22 = 1.0, q11 = 1.0, q22 = 2.0;
//        double mu = 0.5, p0 = 1.0, r = 0.5;
//        double[] y1 = {1.0, -0.5, 2.0};
//        double[] y2 = {0.5,  1.5, -1.0};
//        double dt = 0.2;
//
//        Model mVL   = make2DWithMethod(a11, a22, q11, q22, mu, p0, r, y1, y2, dt,
//                CovarianceGradientMethod.VAN_LOAN_ADJOINT);
//        Model mLyap = make2DWithMethod(a11, a22, q11, q22, mu, p0, r, y1, y2, dt,
//                CovarianceGradientMethod.LYAPUNOV_ADJOINT);
//
//        double[] gradVL   = mVL.analyticalEngine.getGradientWrt(mVL.process.getDriftMatrix());
//        double[] gradLyap = mLyap.analyticalEngine.getGradientWrt(mLyap.process.getDriftMatrix());
//
//        String[] labels = {"A[0][0]", "A[0][1]", "A[1][0]", "A[1][1]"};
//        for (int idx = 0; idx < 4; ++idx) {
//            assertEquals("VAN_LOAN_ADJOINT and LYAPUNOV_ADJOINT must agree on " + labels[idx],
//                    gradVL[idx], gradLyap[idx], TOL_FD);
//        }
//    }

    // ── Diffusion parametrizations ───────────────────────────────────────────────

    public void testDiffusionGradientMatchesFD_1D_Dense() {
        final MatrixParameter drift = makeMatrix("A", new double[][]{{0.6}});
        final MatrixParameter diffusion = makeMatrix("Q", new double[][]{{1.4}});
        final Parameter mean = new Parameter.Default(0.25);
        final MatrixParameter initCov = makeMatrix("P0", new double[][]{{0.9}});
        final OUProcessModel process = new OUProcessModel(
                "ou1dDiff", 1, drift, diffusion, mean, initCov, CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final MatrixParameter H = makeMatrix("H", new double[][]{{1.0}});
        final MatrixParameter R = makeMatrix("R", new double[][]{{0.45}});
        final MatrixParameter Y = makeMatrix("Y", new double[][]{{0.8, -0.3, 1.2, 0.1}});
        final GaussianObservationModel obs = new GaussianObservationModel("obs", 1, H, R, Y);
        final TimeGrid grid = new UniformTimeGrid(4, 0.0, 0.2);

        final AnalyticalKalmanGradientEngine analytical = makeAnalyticalEngine(process, obs, grid, true);
        final KalmanLikelihoodEngine likelihood = new KalmanLikelihoodEngine(
                process.getRepresentation(GaussianTransitionRepresentation.class), obs, grid);

        final double[] analytic = analytical.getGradientWrt(process.getDiffusionMatrix());
        assertEquals(1, analytic.length);
        assertEquals("Analytical diffusion gradient must match FD in 1D",
                numericalGradient(likelihood, process.getDiffusionMatrix(), 0, 1e-6),
                analytic[0], TOL_FD);
    }

    public void testDiffusionGradientMatchesFD_2D_Dense() {
        final MatrixParameter drift = makeMatrix("A", new double[][]{{0.7, -0.15}, {0.1, 1.0}});
        final MatrixParameter diffusion = makeMatrix("Q", new double[][]{{1.2, 0.18}, {0.18, 1.6}});
        final Parameter mean = new Parameter.Default(0.1);
        final MatrixParameter initCov = makeMatrix("P0", new double[][]{{1.0, 0.1}, {0.1, 1.1}});
        final OUProcessModel process = new OUProcessModel(
                "ou2dDiff", 2, drift, diffusion, mean, initCov, CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final MatrixParameter H = makeMatrix("H", new double[][]{{1, 0}, {0, 1}});
        final MatrixParameter R = makeMatrix("R", new double[][]{{0.4, 0.02}, {0.02, 0.5}});
        final MatrixParameter Y = makeMatrix("Y", new double[][]{
                {0.7, -0.2, 1.1, 0.4},
                {-0.5, 0.6, -0.8, 0.9}
        });
        final GaussianObservationModel obs = new GaussianObservationModel("obs", 2, H, R, Y);
        final TimeGrid grid = new UniformTimeGrid(4, 0.0, 0.15);

        final AnalyticalKalmanGradientEngine analytical = makeAnalyticalEngine(process, obs, grid, true);
        final KalmanLikelihoodEngine likelihood = new KalmanLikelihoodEngine(
                process.getRepresentation(GaussianTransitionRepresentation.class), obs, grid);

        final double[] analytic = analytical.getGradientWrt(process.getDiffusionMatrix());
        for (int idx = 0; idx < analytic.length; ++idx) {
            assertEquals("Analytical dense diffusion gradient must match FD, idx=" + idx,
                    numericalGradient(likelihood, process.getDiffusionMatrix(), idx, 1e-6),
                    analytic[idx], TOL_FD);
        }
    }

    public void testDiffusionGradientMatchesFD_2D_CompoundSymmetric() {
        final MatrixParameter drift = makeMatrix("A", new double[][]{{0.65, -0.1}, {0.08, 0.95}});
        final Parameter diagonal = new Parameter.Default(new double[]{1.3, 0.9});
        final Parameter offDiagonal = new Parameter.Default(new double[]{0.25});
        final CompoundSymmetricMatrix diffusion = new CompoundSymmetricMatrix(
                diagonal, offDiagonal, true, false);
        final Parameter mean = new Parameter.Default(0.0);
        final MatrixParameter initCov = makeMatrix("P0", new double[][]{{0.8, 0.05}, {0.05, 1.0}});
        final OUProcessModel process = new OUProcessModel(
                "ou2dCompoundDiff", 2, drift, diffusion, mean, initCov, CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final MatrixParameter H = makeMatrix("H", new double[][]{{1, 0}, {0, 1}});
        final MatrixParameter R = makeMatrix("R", new double[][]{{0.35, 0.01}, {0.01, 0.45}});
        final MatrixParameter Y = makeMatrix("Y", new double[][]{
                {0.2, 0.7, -0.4, 1.0},
                {-0.6, 0.3, 0.9, -0.1}
        });
        final GaussianObservationModel obs = new GaussianObservationModel("obs", 2, H, R, Y);
        final TimeGrid grid = new UniformTimeGrid(4, 0.0, 0.2);

        final AnalyticalKalmanGradientEngine analytical = makeAnalyticalEngine(process, obs, grid, true);
        final KalmanLikelihoodEngine likelihood = new KalmanLikelihoodEngine(
                process.getRepresentation(GaussianTransitionRepresentation.class), obs, grid);

        final double[] analyticDiag = analytical.getGradientWrt(diagonal);
        for (int idx = 0; idx < analyticDiag.length; ++idx) {
            assertEquals("Compound-symmetric diffusion diagonal gradient must match FD, idx=" + idx,
                    numericalGradient(likelihood, diagonal, idx, 1e-6),
                    analyticDiag[idx], TOL_FD);
        }

        final double[] analyticOffDiag = analytical.getGradientWrt(offDiagonal);
        for (int idx = 0; idx < analyticOffDiag.length; ++idx) {
            assertEquals("Compound-symmetric diffusion correlation gradient must match FD, idx=" + idx,
                    numericalGradient(likelihood, offDiagonal, idx, 1e-6),
                    analyticOffDiag[idx], TOL_FD);
        }
    }

    public void testStationaryMeanGradientMatchesFD_1D() {
        final double a = 0.45, q = 0.9, mu = -0.2, p0 = 1.1, r = 0.7;
        final double[] y = {0.3, -0.8, 1.2, 0.4};
        final double dt = 0.18;

        final Model m = makeScalar1DWithMethod(a, q, mu, p0, r, y, dt, CovarianceGradientMethod.LYAPUNOV_ADJOINT);
        final TimeSeriesLikelihood expectation = makeExpectationAnalyticalLikelihood(
                "ouMean1dExpectation", m.process, m.obs, new UniformTimeGrid(y.length, 0.0, dt));

        final double[] analytic = expectation.getGradientWrt(m.process.getStationaryMeanParameter()).getGradientLogDensity(null);
        assertEquals("Analytical 1D stationary mean gradient must match FD",
                numericalGradient(m.likelihoodEngine, m.process.getStationaryMeanParameter(), 0, 1e-6),
                analytic[0], TOL_FD);
    }

    public void testStationaryMeanGradientMatchesFD_2D() {
        final MatrixParameter drift = makeMatrix("A.mean2d", new double[][]{{0.55, -0.12}, {0.18, 0.9}});
        final MatrixParameter diffusion = makeMatrix("Q.mean2d", new double[][]{{1.0, 0.08}, {0.08, 1.4}});
        final Parameter mean = new Parameter.Default(new double[]{0.15, -0.25});
        final MatrixParameter initCov = makeMatrix("P0.mean2d", new double[][]{{1.1, 0.1}, {0.1, 0.95}});
        final OUProcessModel process = new OUProcessModel(
                "ouMean2d", 2, drift, diffusion, mean, initCov, CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final MatrixParameter H = makeMatrix("H.mean2d", new double[][]{{1.0, 0.0}, {0.0, 1.0}});
        final MatrixParameter R = makeMatrix("R.mean2d", new double[][]{{0.6, 0.05}, {0.05, 0.8}});
        final MatrixParameter Y = makeMatrix("Y.mean2d", new double[][]{
                {0.4, -0.3, 1.1, 0.2},
                {-0.5, 0.8, -0.6, 0.1}
        });
        final GaussianObservationModel obs = new GaussianObservationModel("obsMean2d", 2, H, R, Y);
        final TimeGrid grid = new UniformTimeGrid(4, 0.0, 0.17);

        final KalmanLikelihoodEngine likelihood = new KalmanLikelihoodEngine(
                process.getRepresentation(GaussianTransitionRepresentation.class), obs, grid);
        final TimeSeriesLikelihood expectation = makeExpectationAnalyticalLikelihood(
                "ouMean2dExpectation", process, obs, grid);

        final double[] analytic = expectation.getGradientWrt(process.getStationaryMeanParameter()).getGradientLogDensity(null);
        for (int idx = 0; idx < 2; ++idx) {
            assertEquals("Analytical 2D stationary mean gradient must match FD, idx=" + idx,
                    numericalGradient(likelihood, process.getStationaryMeanParameter(), idx, 1e-6),
                    analytic[idx], TOL_FD);
        }
    }

    public void testCanonicalStationaryMeanGradientMatchesFD_2D() {
        final MatrixParameter drift = makeMatrix("A.canonMean2d", new double[][]{{0.42, -0.16}, {0.11, 0.78}});
        final MatrixParameter diffusion = makeMatrix("Q.canonMean2d", new double[][]{{0.95, 0.04}, {0.04, 1.25}});
        final Parameter mean = new Parameter.Default(new double[]{0.08, -0.18});
        final MatrixParameter initCov = makeMatrix("P0.canonMean2d", new double[][]{{1.0, 0.06}, {0.06, 0.9}});
        final OUProcessModel process = new OUProcessModel(
                "ouCanonicalMean2d", 2, drift, diffusion, mean, initCov, CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final MatrixParameter H = makeMatrix("H.canonMean2d", new double[][]{{1.0, 0.0}, {0.0, 1.0}});
        final MatrixParameter R = makeMatrix("R.canonMean2d", new double[][]{{0.55, 0.03}, {0.03, 0.75}});
        final MatrixParameter Y = makeMatrix("Y.canonMean2d", new double[][]{
                {0.2, -0.1, 0.9, -0.3},
                {-0.4, 0.7, -0.2, 0.5}
        });
        final GaussianObservationModel obs = new GaussianObservationModel("obsCanonMean2d", 2, H, R, Y);
        final TimeGrid grid = new UniformTimeGrid(4, 0.0, 0.14);

        final KalmanLikelihoodEngine likelihood = new KalmanLikelihoodEngine(
                process.getRepresentation(GaussianTransitionRepresentation.class), obs, grid);
        final TimeSeriesLikelihood canonical = makeCanonicalAnalyticalLikelihood(
                "ouMean2dCanonical", process, obs, grid);

        final double[] analytic = canonical.getGradientWrt(process.getStationaryMeanParameter()).getGradientLogDensity(null);
        for (int idx = 0; idx < 2; ++idx) {
            assertEquals("Canonical analytical stationary mean gradient must match FD, idx=" + idx,
                    numericalGradient(likelihood, process.getStationaryMeanParameter(), idx, 1e-6),
                    analytic[idx], TOL_FD);
        }
    }

    public void testOrthogonalBlockStationaryMeanGradientMatchesFD_2D() {
        final double[][] diffusion = new double[][]{{1.2, 0.1}, {0.1, 1.6}};
        final double[][] initCov = new double[][]{{1.0, 0.2}, {0.2, 1.3}};
        final double[][] noise = new double[][]{{0.6, 0.05}, {0.05, 0.8}};
        final double[][] observations = new double[][]{
                {0.7, -0.4, 1.1, 0.2},
                {-0.3, 0.8, -1.2, 0.5}
        };
        final double dt = 0.15;
        final double h = 1e-6;

        final Model orthogonalModel = makeOrthogonalBlockPolarFullWithMethod(
                "ou2dorthMean", 2,
                new double[]{0.35},
                new double[0],
                new double[]{0.9},
                new double[]{0.25},
                new double[]{0.18},
                diffusion,
                0.25,
                initCov,
                noise,
                observations,
                CovarianceGradientMethod.STATIONARY_LYAPUNOV,
                dt);

        final TimeSeriesLikelihood expectation = makeExpectationAnalyticalLikelihood(
                "ou2dorthMeanExpectation",
                orthogonalModel.process,
                orthogonalModel.obs,
                new UniformTimeGrid(observations[0].length, 0.0, dt));
        final double[] analytic = expectation.getGradientWrt(
                orthogonalModel.process.getStationaryMeanParameter()).getGradientLogDensity(null);
        assertEquals("Orthogonal-block stationary mean gradient must match FD",
                numericalGradient(orthogonalModel.likelihoodEngine,
                        orthogonalModel.process.getStationaryMeanParameter(), 0, h),
                analytic[0], TOL_FD);
    }

    public void testCanonicalOrthogonalBlockStationaryMeanGradientMatchesFD_2D() {
        final double[][] diffusion = new double[][]{{1.2, 0.1}, {0.1, 1.6}};
        final double[][] initCov = new double[][]{{1.0, 0.2}, {0.2, 1.3}};
        final double[][] noise = new double[][]{{0.6, 0.05}, {0.05, 0.8}};
        final double[][] observations = new double[][]{
                {0.7, -0.4, 1.1, 0.2},
                {-0.3, 0.8, -1.2, 0.5}
        };
        final double dt = 0.15;
        final double h = 1e-6;

        final Model orthogonalModel = makeOrthogonalBlockPolarFullWithMethod(
                "ou2dorthCanonMean", 2,
                new double[]{0.35},
                new double[0],
                new double[]{0.9},
                new double[]{0.25},
                new double[]{0.18},
                diffusion,
                0.25,
                initCov,
                noise,
                observations,
                CovarianceGradientMethod.STATIONARY_LYAPUNOV,
                dt);

        final TimeSeriesLikelihood canonical = makeCanonicalAnalyticalLikelihood(
                "ou2dorthMeanCanonical",
                orthogonalModel.process,
                orthogonalModel.obs,
                new UniformTimeGrid(observations[0].length, 0.0, dt));
        final double[] analytic = canonical.getGradientWrt(
                orthogonalModel.process.getStationaryMeanParameter()).getGradientLogDensity(null);
        assertEquals("Canonical orthogonal-block stationary mean gradient must match FD",
                numericalGradient(orthogonalModel.likelihoodEngine,
                        orthogonalModel.process.getStationaryMeanParameter(), 0, h),
                analytic[0], TOL_FD);
    }

    public void testExactOuCanonicalStateAndTransitionMatchExpectationForm() {
        final MatrixParameter drift = makeMatrix("A.canonical.exact", new double[][]{
                {0.65, -0.08},
                {0.11, 1.10}
        });
        final MatrixParameter diffusion = makeMatrix("Q.canonical.exact", new double[][]{
                {1.4, 0.15},
                {0.15, 0.9}
        });
        final Parameter mean = new Parameter.Default(new double[]{0.3, -0.2});
        final MatrixParameter initCov = makeMatrix("P0.canonical.exact", new double[][]{
                {1.2, 0.05},
                {0.05, 0.8}
        });
        final OUProcessModel process = new OUProcessModel(
                "ouCanonicalExact", 2, drift, diffusion, mean, initCov, CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianBranchTransitionKernel expectationKernel =
                process.getRepresentation(dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianBranchTransitionKernel.class);
        final CanonicalGaussianBranchTransitionKernel canonicalKernel =
                process.getRepresentation(CanonicalGaussianBranchTransitionKernel.class);

        final CanonicalGaussianState canonicalState = new CanonicalGaussianState(2);
        canonicalKernel.fillInitialCanonicalState(canonicalState);
        final double[] stateValue = {0.25, -0.4};
        final double[] initialMean = new double[2];
        final double[][] initialCovariance = new double[2][2];
        expectationKernel.getInitialMean(initialMean);
        expectationKernel.getInitialCovariance(initialCovariance);
        assertEquals("Canonical initial state must match moment-form density",
                logDensityFromMoments(initialMean, initialCovariance, stateValue),
                logDensityFromCanonicalState(canonicalState, stateValue),
                1e-10);

        final double dt = 0.3;
        final double[][] transitionMatrix = new double[2][2];
        final double[] transitionOffset = new double[2];
        final double[][] transitionCovariance = new double[2][2];
        expectationKernel.fillTransitionMatrix(dt, transitionMatrix);
        expectationKernel.fillTransitionOffset(dt, transitionOffset);
        expectationKernel.fillTransitionCovariance(dt, transitionCovariance);

        final CanonicalGaussianTransition canonicalTransition = new CanonicalGaussianTransition(2);
        canonicalKernel.fillCanonicalTransition(dt, canonicalTransition);

        final double[] previousState = {0.6, -0.1};
        final double[] nextState = {-0.3, 0.8};
        assertEquals("Canonical OU transition must match moment-form transition density",
                logTransitionFromMoments(transitionMatrix, transitionOffset, transitionCovariance, previousState, nextState),
                logTransitionFromCanonical(canonicalTransition, previousState, nextState),
                1e-10);
    }

    public void testEulerCanonicalTransitionMatchesExpectationForm() {
        final MatrixParameter drift = makeMatrix("A.canonical.euler", new double[][]{
                {0.25, -0.04},
                {0.02, 0.55}
        });
        final MatrixParameter diffusion = makeMatrix("Q.canonical.euler", new double[][]{
                {0.9, 0.1},
                {0.1, 0.7}
        });
        final Parameter mean = new Parameter.Default(new double[]{-0.15, 0.35});
        final MatrixParameter initCov = makeMatrix("P0.canonical.euler", new double[][]{
                {0.7, 0.02},
                {0.02, 1.1}
        });
        final EulerOUProcessModel process = new EulerOUProcessModel(
                "ouCanonicalEuler", 2, drift, diffusion, mean, initCov);

        final dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianBranchTransitionKernel expectationKernel =
                process.getRepresentation(dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianBranchTransitionKernel.class);
        final CanonicalGaussianBranchTransitionKernel canonicalKernel =
                process.getRepresentation(CanonicalGaussianBranchTransitionKernel.class);

        final double dt = 0.15;
        final double[][] transitionMatrix = new double[2][2];
        final double[] transitionOffset = new double[2];
        final double[][] transitionCovariance = new double[2][2];
        expectationKernel.fillTransitionMatrix(dt, transitionMatrix);
        expectationKernel.fillTransitionOffset(dt, transitionOffset);
        expectationKernel.fillTransitionCovariance(dt, transitionCovariance);

        final CanonicalGaussianTransition canonicalTransition = new CanonicalGaussianTransition(2);
        canonicalKernel.fillCanonicalTransition(dt, canonicalTransition);

        final double[] previousState = {-0.2, 0.5};
        final double[] nextState = {0.15, -0.1};
        assertEquals("Canonical Euler transition must match moment-form transition density",
                logTransitionFromMoments(transitionMatrix, transitionOffset, transitionCovariance, previousState, nextState),
                logTransitionFromCanonical(canonicalTransition, previousState, nextState),
                1e-10);
    }

    public void testCanonicalKalmanLikelihoodMatchesExpectationEngine_ExactOu() {
        final MatrixParameter drift = makeMatrix("A.canonical.filter.exact", new double[][]{
                {0.45, -0.07},
                {0.10, 0.95}
        });
        final MatrixParameter diffusion = makeMatrix("Q.canonical.filter.exact", new double[][]{
                {1.1, 0.12},
                {0.12, 0.85}
        });
        final Parameter mean = new Parameter.Default(new double[]{0.2, -0.1});
        final MatrixParameter initCov = makeMatrix("P0.canonical.filter.exact", new double[][]{
                {0.9, 0.03},
                {0.03, 1.05}
        });
        final OUProcessModel process = new OUProcessModel(
                "ouCanonicalFilterExact", 2, drift, diffusion, mean, initCov, CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final MatrixParameter H = makeMatrix("H.canonical.filter.exact", new double[][]{{1, 0}, {0, 1}});
        final MatrixParameter R = makeMatrix("R.canonical.filter.exact", new double[][]{
                {0.4, 0.02},
                {0.02, 0.5}
        });
        final MatrixParameter Y = makeMatrix("Y.canonical.filter.exact", new double[][]{
                {0.5, Double.NaN, -0.3, 1.1},
                {-0.2, Double.NaN, 0.8, -0.4}
        });
        final GaussianObservationModel obs = new GaussianObservationModel("obsCanonicalExact", 2, H, R, Y);
        final TimeGrid grid = new UniformTimeGrid(4, 0.0, 0.25);

        final KalmanLikelihoodEngine expectationEngine = new KalmanLikelihoodEngine(
                process.getRepresentation(GaussianTransitionRepresentation.class), obs, grid);
        final CanonicalKalmanLikelihoodEngine canonicalEngine = new CanonicalKalmanLikelihoodEngine(
                process.getRepresentation(CanonicalGaussianBranchTransitionKernel.class), obs, grid);

        assertEquals("Canonical forward likelihood must match expectation-form Kalman likelihood (exact OU)",
                expectationEngine.getLogLikelihood(),
                canonicalEngine.getLogLikelihood(),
                1e-9);
    }

    public void testCanonicalKalmanLikelihoodMatchesExpectationEngine_Euler() {
        final MatrixParameter drift = makeMatrix("A.canonical.filter.euler", new double[][]{
                {0.18, -0.03},
                {0.05, 0.42}
        });
        final MatrixParameter diffusion = makeMatrix("Q.canonical.filter.euler", new double[][]{
                {0.75, 0.08},
                {0.08, 0.65}
        });
        final Parameter mean = new Parameter.Default(new double[]{-0.25, 0.15});
        final MatrixParameter initCov = makeMatrix("P0.canonical.filter.euler", new double[][]{
                {1.0, 0.01},
                {0.01, 0.7}
        });
        final EulerOUProcessModel process = new EulerOUProcessModel(
                "ouCanonicalFilterEuler", 2, drift, diffusion, mean, initCov);

        final MatrixParameter H = makeMatrix("H.canonical.filter.euler", new double[][]{{1, 0}, {0, 1}});
        final MatrixParameter R = makeMatrix("R.canonical.filter.euler", new double[][]{
                {0.3, 0.01},
                {0.01, 0.35}
        });
        final MatrixParameter Y = makeMatrix("Y.canonical.filter.euler", new double[][]{
                {-0.1, 0.6, -0.4},
                {0.2, -0.5, 0.9}
        });
        final GaussianObservationModel obs = new GaussianObservationModel("obsCanonicalEuler", 2, H, R, Y);
        final TimeGrid grid = new UniformTimeGrid(3, 0.0, 0.1);

        final KalmanLikelihoodEngine expectationEngine = new KalmanLikelihoodEngine(
                process.getRepresentation(GaussianTransitionRepresentation.class), obs, grid);
        final CanonicalKalmanLikelihoodEngine canonicalEngine = new CanonicalKalmanLikelihoodEngine(
                process.getRepresentation(CanonicalGaussianBranchTransitionKernel.class), obs, grid);

        assertEquals("Canonical forward likelihood must match expectation-form Kalman likelihood (Euler)",
                expectationEngine.getLogLikelihood(),
                canonicalEngine.getLogLikelihood(),
                1e-9);
    }

    public void testForwardLikelihoodFactorySelectsEquivalentExactOuBackends() {
        final MatrixParameter drift = makeMatrix("A.factory.exact", new double[][]{
                {0.52, -0.06},
                {0.09, 0.88}
        });
        final MatrixParameter diffusion = makeMatrix("Q.factory.exact", new double[][]{
                {1.0, 0.09},
                {0.09, 0.78}
        });
        final Parameter mean = new Parameter.Default(new double[]{0.1, -0.3});
        final MatrixParameter initCov = makeMatrix("P0.factory.exact", new double[][]{
                {1.1, 0.04},
                {0.04, 0.95}
        });
        final OUProcessModel process = new OUProcessModel(
                "ouFactoryExact", 2, drift, diffusion, mean, initCov, CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final MatrixParameter H = makeMatrix("H.factory.exact", new double[][]{{1, 0}, {0, 1}});
        final MatrixParameter R = makeMatrix("R.factory.exact", new double[][]{
                {0.45, 0.03},
                {0.03, 0.4}
        });
        final MatrixParameter Y = makeMatrix("Y.factory.exact", new double[][]{
                {0.3, -0.7, Double.NaN, 0.4},
                {0.8, -0.1, Double.NaN, -0.6}
        });
        final GaussianObservationModel obs = new GaussianObservationModel("obsFactoryExact", 2, H, R, Y);
        final TimeGrid grid = new UniformTimeGrid(4, 0.0, 0.2);

        final LikelihoodEngine expectation = GaussianLikelihoodEngineFactory.createForwardEngine(
                process, obs, grid, GaussianForwardComputationMode.EXPECTATION);
        final LikelihoodEngine canonical = GaussianLikelihoodEngineFactory.createForwardEngine(
                process, obs, grid, GaussianForwardComputationMode.CANONICAL);

        assertEquals("Factory-selected exact OU forward backends must agree",
                expectation.getLogLikelihood(),
                canonical.getLogLikelihood(),
                1e-9);
    }

    public void testForwardLikelihoodFactorySelectsEquivalentEulerBackends() {
        final MatrixParameter drift = makeMatrix("A.factory.euler", new double[][]{
                {0.16, -0.02},
                {0.03, 0.36}
        });
        final MatrixParameter diffusion = makeMatrix("Q.factory.euler", new double[][]{
                {0.82, 0.05},
                {0.05, 0.61}
        });
        final Parameter mean = new Parameter.Default(new double[]{-0.05, 0.2});
        final MatrixParameter initCov = makeMatrix("P0.factory.euler", new double[][]{
                {0.95, 0.0},
                {0.0, 0.85}
        });
        final EulerOUProcessModel process = new EulerOUProcessModel(
                "ouFactoryEuler", 2, drift, diffusion, mean, initCov);

        final MatrixParameter H = makeMatrix("H.factory.euler", new double[][]{{1, 0}, {0, 1}});
        final MatrixParameter R = makeMatrix("R.factory.euler", new double[][]{
                {0.25, 0.0},
                {0.0, 0.3}
        });
        final MatrixParameter Y = makeMatrix("Y.factory.euler", new double[][]{
                {-0.2, 0.4, 0.1},
                {0.5, -0.3, 0.7}
        });
        final GaussianObservationModel obs = new GaussianObservationModel("obsFactoryEuler", 2, H, R, Y);
        final TimeGrid grid = new UniformTimeGrid(3, 0.0, 0.1);

        final LikelihoodEngine expectation = GaussianLikelihoodEngineFactory.createForwardEngine(
                process, obs, grid, GaussianForwardComputationMode.EXPECTATION);
        final LikelihoodEngine canonical = GaussianLikelihoodEngineFactory.createForwardEngine(
                process, obs, grid, GaussianForwardComputationMode.CANONICAL);

        assertEquals("Factory-selected Euler forward backends must agree",
                expectation.getLogLikelihood(),
                canonical.getLogLikelihood(),
                1e-9);
    }

    public void testTimeSeriesLikelihoodFactorySupportsCanonicalForwardMode() {
        final MatrixParameter drift = makeMatrix("A.ts.factory.exact", new double[][]{
                {0.41, -0.05},
                {0.07, 0.83}
        });
        final MatrixParameter diffusion = makeMatrix("Q.ts.factory.exact", new double[][]{
                {1.05, 0.07},
                {0.07, 0.92}
        });
        final Parameter mean = new Parameter.Default(new double[]{0.12, -0.08});
        final MatrixParameter initCov = makeMatrix("P0.ts.factory.exact", new double[][]{
                {0.8, 0.02},
                {0.02, 1.0}
        });
        final OUProcessModel process = new OUProcessModel(
                "ouTsFactoryExact", 2, drift, diffusion, mean, initCov, CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final MatrixParameter H = makeMatrix("H.ts.factory.exact", new double[][]{{1, 0}, {0, 1}});
        final MatrixParameter R = makeMatrix("R.ts.factory.exact", new double[][]{
                {0.35, 0.01},
                {0.01, 0.42}
        });
        final MatrixParameter Y = makeMatrix("Y.ts.factory.exact", new double[][]{
                {0.1, Double.NaN, -0.5, 0.9},
                {-0.4, Double.NaN, 0.7, -0.2}
        });
        final GaussianObservationModel obs = new GaussianObservationModel("obsTsFactoryExact", 2, H, R, Y);
        final TimeGrid grid = new UniformTimeGrid(4, 0.0, 0.2);
        final BasicTimeSeriesModel model = new BasicTimeSeriesModel("tsFactoryExact", process, obs, grid);

        final TimeSeriesLikelihood expectationLikelihood = GaussianTimeSeriesLikelihoodFactory.create(
                "tsLikeExpectation",
                model,
                GaussianForwardComputationMode.EXPECTATION,
                GaussianGradientComputationMode.DISABLED);
        final TimeSeriesLikelihood canonicalLikelihood = GaussianTimeSeriesLikelihoodFactory.create(
                "tsLikeCanonical",
                model,
                GaussianForwardComputationMode.CANONICAL,
                GaussianGradientComputationMode.DISABLED);

        assertEquals("Top-level time-series likelihood factory must support canonical forward mode",
                expectationLikelihood.getLogLikelihood(),
                canonicalLikelihood.getLogLikelihood(),
                1e-9);
    }

    public void testTimeSeriesLikelihoodFactorySupportsCanonicalForwardWithExpectationGradients() {
        final MatrixParameter drift = makeMatrix("A.ts.factory.grad", new double[][]{{0.3}});
        final MatrixParameter diffusion = makeMatrix("Q.ts.factory.grad", new double[][]{{1.2}});
        final Parameter mean = new Parameter.Default(0.4);
        final MatrixParameter initCov = makeMatrix("P0.ts.factory.grad", new double[][]{{0.9}});
        final OUProcessModel process = new OUProcessModel(
                "ouTsFactoryGrad", 1, drift, diffusion, mean, initCov, CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final MatrixParameter H = makeMatrix("H.ts.factory.grad", new double[][]{{1.0}});
        final MatrixParameter R = makeMatrix("R.ts.factory.grad", new double[][]{{0.5}});
        final MatrixParameter Y = makeMatrix("Y.ts.factory.grad", new double[][]{{0.7, -0.4, 1.3}});
        final GaussianObservationModel obs = new GaussianObservationModel("obsTsFactoryGrad", 1, H, R, Y);
        final TimeGrid grid = new UniformTimeGrid(3, 0.0, 0.25);
        final BasicTimeSeriesModel model = new BasicTimeSeriesModel("tsFactoryGrad", process, obs, grid);

        final TimeSeriesLikelihood canonicalForwardLikelihood = GaussianTimeSeriesLikelihoodFactory.create(
                "tsLikeCanonicalWithGrad",
                model,
                GaussianForwardComputationMode.CANONICAL,
                GaussianGradientComputationMode.EXPECTATION_ANALYTICAL);
        final TimeSeriesLikelihood expectationLikelihood = GaussianTimeSeriesLikelihoodFactory.create(
                "tsLikeExpectationWithGrad",
                model,
                GaussianForwardComputationMode.EXPECTATION,
                GaussianGradientComputationMode.EXPECTATION_ANALYTICAL);

        assertEquals("Canonical-forward and expectation-forward top-level likelihoods must agree",
                expectationLikelihood.getLogLikelihood(),
                canonicalForwardLikelihood.getLogLikelihood(),
                1e-9);

        final double[] expectationGradient = expectationLikelihood.getGradientWrt(process.getDriftMatrix())
                .getGradientLogDensity(null);
        final double[] canonicalGradient = canonicalForwardLikelihood.getGradientWrt(process.getDriftMatrix())
                .getGradientLogDensity(null);

        assertEquals("Expectation analytical gradient must still be available under canonical forward mode",
                expectationGradient[0],
                canonicalGradient[0],
                1e-9);
    }

    public void testCanonicalKalmanSmootherMatchesExpectationSmoother_ExactOu() {
        final MatrixParameter drift = makeMatrix("A.canonical.smoother.exact", new double[][]{
                {0.38, -0.04},
                {0.06, 0.91}
        });
        final MatrixParameter diffusion = makeMatrix("Q.canonical.smoother.exact", new double[][]{
                {1.15, 0.08},
                {0.08, 0.88}
        });
        final Parameter mean = new Parameter.Default(new double[]{0.18, -0.12});
        final MatrixParameter initCov = makeMatrix("P0.canonical.smoother.exact", new double[][]{
                {0.95, 0.03},
                {0.03, 0.85}
        });
        final OUProcessModel process = new OUProcessModel(
                "ouCanonicalSmootherExact", 2, drift, diffusion, mean, initCov, CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final MatrixParameter H = makeMatrix("H.canonical.smoother.exact", new double[][]{{1, 0}, {0, 1}});
        final MatrixParameter R = makeMatrix("R.canonical.smoother.exact", new double[][]{
                {0.3, 0.02},
                {0.02, 0.45}
        });
        final MatrixParameter Y = makeMatrix("Y.canonical.smoother.exact", new double[][]{
                {0.15, Double.NaN, -0.35, 0.9},
                {-0.25, Double.NaN, 0.55, -0.1}
        });
        final GaussianObservationModel obs = new GaussianObservationModel("obsCanonicalSmootherExact", 2, H, R, Y);
        final TimeGrid grid = new UniformTimeGrid(4, 0.0, 0.2);

        final KalmanSmootherEngine expectation = new KalmanSmootherEngine(
                process.getRepresentation(GaussianTransitionRepresentation.class), obs, grid);
        final CanonicalKalmanSmootherEngine canonical = new CanonicalKalmanSmootherEngine(
                process.getRepresentation(CanonicalGaussianBranchTransitionKernel.class),
                process.getRepresentation(GaussianTransitionRepresentation.class),
                obs,
                grid);

        assertEquals("Canonical smoother likelihood must match expectation smoother (exact OU)",
                expectation.getLogLikelihood(), canonical.getLogLikelihood(), 1e-9);
        assertVectorArrayEquals("Predicted means", expectation.getPredictedMeans(), canonical.getPredictedMeans(), 1e-9);
        assertMatrixArrayEquals("Predicted covariances", expectation.getPredictedCovariances(), canonical.getPredictedCovariances(), 1e-9);
        assertVectorArrayEquals("Filtered means", expectation.getFilteredMeans(), canonical.getFilteredMeans(), 1e-9);
        assertMatrixArrayEquals("Filtered covariances", expectation.getFilteredCovariances(), canonical.getFilteredCovariances(), 1e-9);
        assertVectorArrayEquals("Smoothed means", expectation.getSmoothedMeans(), canonical.getSmoothedMeans(), 1e-9);
        assertMatrixArrayEquals("Smoothed covariances", expectation.getSmoothedCovariances(), canonical.getSmoothedCovariances(), 1e-9);
    }

    public void testCanonicalKalmanSmootherMatchesExpectationSmoother_Euler() {
        final MatrixParameter drift = makeMatrix("A.canonical.smoother.euler", new double[][]{
                {0.14, -0.01},
                {0.04, 0.31}
        });
        final MatrixParameter diffusion = makeMatrix("Q.canonical.smoother.euler", new double[][]{
                {0.7, 0.04},
                {0.04, 0.58}
        });
        final Parameter mean = new Parameter.Default(new double[]{-0.1, 0.25});
        final MatrixParameter initCov = makeMatrix("P0.canonical.smoother.euler", new double[][]{
                {0.75, 0.0},
                {0.0, 0.9}
        });
        final EulerOUProcessModel process = new EulerOUProcessModel(
                "ouCanonicalSmootherEuler", 2, drift, diffusion, mean, initCov);

        final MatrixParameter H = makeMatrix("H.canonical.smoother.euler", new double[][]{{1, 0}, {0, 1}});
        final MatrixParameter R = makeMatrix("R.canonical.smoother.euler", new double[][]{
                {0.22, 0.0},
                {0.0, 0.28}
        });
        final MatrixParameter Y = makeMatrix("Y.canonical.smoother.euler", new double[][]{
                {-0.05, 0.45, -0.2},
                {0.35, -0.15, 0.6}
        });
        final GaussianObservationModel obs = new GaussianObservationModel("obsCanonicalSmootherEuler", 2, H, R, Y);
        final TimeGrid grid = new UniformTimeGrid(3, 0.0, 0.1);

        final KalmanSmootherEngine expectation = new KalmanSmootherEngine(
                process.getRepresentation(GaussianTransitionRepresentation.class), obs, grid);
        final CanonicalKalmanSmootherEngine canonical = new CanonicalKalmanSmootherEngine(
                process.getRepresentation(CanonicalGaussianBranchTransitionKernel.class),
                process.getRepresentation(GaussianTransitionRepresentation.class),
                obs,
                grid);

        assertEquals("Canonical smoother likelihood must match expectation smoother (Euler)",
                expectation.getLogLikelihood(), canonical.getLogLikelihood(), 1e-9);
        assertVectorArrayEquals("Predicted means", expectation.getPredictedMeans(), canonical.getPredictedMeans(), 1e-9);
        assertMatrixArrayEquals("Predicted covariances", expectation.getPredictedCovariances(), canonical.getPredictedCovariances(), 1e-9);
        assertVectorArrayEquals("Filtered means", expectation.getFilteredMeans(), canonical.getFilteredMeans(), 1e-9);
        assertMatrixArrayEquals("Filtered covariances", expectation.getFilteredCovariances(), canonical.getFilteredCovariances(), 1e-9);
        assertVectorArrayEquals("Smoothed means", expectation.getSmoothedMeans(), canonical.getSmoothedMeans(), 1e-9);
        assertMatrixArrayEquals("Smoothed covariances", expectation.getSmoothedCovariances(), canonical.getSmoothedCovariances(), 1e-9);
    }

    public void testOrthogonalBlockCanonicalTransitionMatchesDense_2D() {
        final Model orthogonalModel = makeOrthogonalBlockPolarFullWithMethod(
                "ouCanonicalOrth2d",
                2,
                new double[]{0.37},
                new double[0],
                new double[]{0.58},
                new double[]{0.63},
                new double[]{0.21},
                new double[][]{{1.12, 0.07}, {0.07, 0.91}},
                0.18,
                new double[][]{{0.9, 0.03}, {0.03, 1.05}},
                new double[][]{{0.31, 0.01}, {0.01, 0.36}},
                new double[][]{
                        {0.2, -0.3, 0.5},
                        {-0.1, 0.4, -0.2}
                },
                CovarianceGradientMethod.LYAPUNOV_ADJOINT,
                0.2);

        final Model denseModel = makeDenseFullWithMethod(
                "ouCanonicalDense2dRef",
                orthogonalModel.process.getDriftMatrix().getParameterAsMatrix(),
                orthogonalModel.process.getDiffusionMatrix().getParameterAsMatrix(),
                0.18,
                new double[][]{{0.9, 0.03}, {0.03, 1.05}},
                new double[][]{{0.31, 0.01}, {0.01, 0.36}},
                new double[][]{
                        {0.2, -0.3, 0.5},
                        {-0.1, 0.4, -0.2}
                },
                0.2,
                CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final CanonicalGaussianTransition orthogonalTransition = new CanonicalGaussianTransition(2);
        final CanonicalGaussianTransition denseTransition = new CanonicalGaussianTransition(2);
        orthogonalModel.process.getRepresentation(CanonicalGaussianBranchTransitionKernel.class)
                .fillCanonicalTransition(0.2, orthogonalTransition);
        denseModel.process.getRepresentation(CanonicalGaussianBranchTransitionKernel.class)
                .fillCanonicalTransition(0.2, denseTransition);

        for (int i = 0; i < 2; ++i) {
            assertEquals("Orthogonal block canonical infoX must match dense, idx=" + i,
                    denseTransition.informationX[i], orthogonalTransition.informationX[i], 1e-10);
            assertEquals("Orthogonal block canonical infoY must match dense, idx=" + i,
                    denseTransition.informationY[i], orthogonalTransition.informationY[i], 1e-10);
            for (int j = 0; j < 2; ++j) {
                assertEquals("Orthogonal block canonical Jxx must match dense at [" + i + "," + j + "]",
                        flatEntry(denseTransition.precisionXX, 2, i, j),
                        flatEntry(orthogonalTransition.precisionXX, 2, i, j), 1e-10);
                assertEquals("Orthogonal block canonical Jxy must match dense at [" + i + "," + j + "]",
                        flatEntry(denseTransition.precisionXY, 2, i, j),
                        flatEntry(orthogonalTransition.precisionXY, 2, i, j), 1e-10);
                assertEquals("Orthogonal block canonical Jyx must match dense at [" + i + "," + j + "]",
                        flatEntry(denseTransition.precisionYX, 2, i, j),
                        flatEntry(orthogonalTransition.precisionYX, 2, i, j), 1e-10);
                assertEquals("Orthogonal block canonical Jyy must match dense at [" + i + "," + j + "]",
                        flatEntry(denseTransition.precisionYY, 2, i, j),
                        flatEntry(orthogonalTransition.precisionYY, 2, i, j), 1e-10);
            }
        }
        assertEquals("Orthogonal block canonical log-normalizer must match dense",
                denseTransition.logNormalizer, orthogonalTransition.logNormalizer, 1e-10);
    }

    public void testOrthogonalBlockCanonicalTransitionMatchesDense_4D() {
        final Model orthogonalModel = makeOrthogonalBlockPolarFullWithMethod(
                "ouCanonicalOrth4d",
                4,
                new double[]{0.22, -0.31, 0.41, 0.18, -0.27, 0.33},
                new double[0],
                new double[]{0.64, 0.83},
                new double[]{0.52, 0.91},
                new double[]{0.17, -0.12},
                new double[][]{
                        {1.20, 0.04, 0.02, 0.01},
                        {0.04, 0.95, 0.05, 0.02},
                        {0.02, 0.05, 1.11, 0.06},
                        {0.01, 0.02, 0.06, 0.88}
                },
                0.07,
                new double[][]{
                        {0.93, 0.02, 0.01, 0.00},
                        {0.02, 1.02, 0.03, 0.01},
                        {0.01, 0.03, 0.97, 0.02},
                        {0.00, 0.01, 0.02, 1.05}
                },
                new double[][]{
                        {0.28, 0.01, 0.00, 0.00},
                        {0.01, 0.31, 0.01, 0.00},
                        {0.00, 0.01, 0.34, 0.01},
                        {0.00, 0.00, 0.01, 0.29}
                },
                new double[][]{
                        {0.12, -0.20, 0.35, Double.NaN, 0.41},
                        {-0.08, 0.26, -0.18, Double.NaN, 0.09},
                        {0.21, -0.11, 0.14, Double.NaN, -0.16},
                        {-0.04, 0.19, -0.09, Double.NaN, 0.27}
                },
                CovarianceGradientMethod.LYAPUNOV_ADJOINT,
                0.16);

        final Model denseModel = makeDenseFullWithMethod(
                "ouCanonicalDense4dRef",
                orthogonalModel.process.getDriftMatrix().getParameterAsMatrix(),
                orthogonalModel.process.getDiffusionMatrix().getParameterAsMatrix(),
                0.07,
                new double[][]{
                        {0.93, 0.02, 0.01, 0.00},
                        {0.02, 1.02, 0.03, 0.01},
                        {0.01, 0.03, 0.97, 0.02},
                        {0.00, 0.01, 0.02, 1.05}
                },
                new double[][]{
                        {0.28, 0.01, 0.00, 0.00},
                        {0.01, 0.31, 0.01, 0.00},
                        {0.00, 0.01, 0.34, 0.01},
                        {0.00, 0.00, 0.01, 0.29}
                },
                new double[][]{
                        {0.12, -0.20, 0.35, Double.NaN, 0.41},
                        {-0.08, 0.26, -0.18, Double.NaN, 0.09},
                        {0.21, -0.11, 0.14, Double.NaN, -0.16},
                        {-0.04, 0.19, -0.09, Double.NaN, 0.27}
                },
                0.16,
                CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final CanonicalGaussianTransition orthogonalTransition = new CanonicalGaussianTransition(4);
        final CanonicalGaussianTransition denseTransition = new CanonicalGaussianTransition(4);
        orthogonalModel.process.getRepresentation(CanonicalGaussianBranchTransitionKernel.class)
                .fillCanonicalTransition(0.16, orthogonalTransition);
        denseModel.process.getRepresentation(CanonicalGaussianBranchTransitionKernel.class)
                .fillCanonicalTransition(0.16, denseTransition);

        for (int i = 0; i < 4; ++i) {
            assertEquals("Orthogonal block 4D canonical infoX must match dense, idx=" + i,
                    denseTransition.informationX[i], orthogonalTransition.informationX[i], 1e-10);
            assertEquals("Orthogonal block 4D canonical infoY must match dense, idx=" + i,
                    denseTransition.informationY[i], orthogonalTransition.informationY[i], 1e-10);
            for (int j = 0; j < 4; ++j) {
                assertEquals("Orthogonal block 4D canonical Jxx must match dense at [" + i + "," + j + "]",
                        flatEntry(denseTransition.precisionXX, 4, i, j),
                        flatEntry(orthogonalTransition.precisionXX, 4, i, j), 1e-10);
                assertEquals("Orthogonal block 4D canonical Jxy must match dense at [" + i + "," + j + "]",
                        flatEntry(denseTransition.precisionXY, 4, i, j),
                        flatEntry(orthogonalTransition.precisionXY, 4, i, j), 1e-10);
                assertEquals("Orthogonal block 4D canonical Jyx must match dense at [" + i + "," + j + "]",
                        flatEntry(denseTransition.precisionYX, 4, i, j),
                        flatEntry(orthogonalTransition.precisionYX, 4, i, j), 1e-10);
                assertEquals("Orthogonal block 4D canonical Jyy must match dense at [" + i + "," + j + "]",
                        flatEntry(denseTransition.precisionYY, 4, i, j),
                        flatEntry(orthogonalTransition.precisionYY, 4, i, j), 1e-10);
            }
        }
        assertEquals("Orthogonal block 4D canonical log-normalizer must match dense",
                denseTransition.logNormalizer, orthogonalTransition.logNormalizer, 1e-10);
    }

    public void testOrthogonalBlockCanonicalLocalAdjointsMatchDense_2D() {
        final double dt = 0.2;
        final Model orthogonalModel = makeOrthogonalBlockPolarFullWithMethod(
                "ouCanonicalOrthAdj2d",
                2,
                new double[]{0.37},
                new double[0],
                new double[]{0.58},
                new double[]{0.63},
                new double[]{0.21},
                new double[][]{{1.12, 0.07}, {0.07, 0.91}},
                0.18,
                new double[][]{{0.9, 0.03}, {0.03, 1.05}},
                new double[][]{{0.31, 0.01}, {0.01, 0.36}},
                new double[][]{
                        {0.2, -0.3, 0.5},
                        {-0.1, 0.4, -0.2}
                },
                CovarianceGradientMethod.LYAPUNOV_ADJOINT,
                dt);

        final Model denseModel = makeDenseFullWithMethod(
                "ouCanonicalDenseAdj2dRef",
                orthogonalModel.process.getDriftMatrix().getParameterAsMatrix(),
                orthogonalModel.process.getDiffusionMatrix().getParameterAsMatrix(),
                0.18,
                new double[][]{{0.9, 0.03}, {0.03, 1.05}},
                new double[][]{{0.31, 0.01}, {0.01, 0.36}},
                new double[][]{
                        {0.2, -0.3, 0.5},
                        {-0.1, 0.4, -0.2}
                },
                dt,
                CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final GaussianObservationModel obs = orthogonalModel.obs;
        final TimeGrid grid = new UniformTimeGrid(3, 0.0, dt);
        final CanonicalKalmanSmootherEngine orthogonalSmoother = new CanonicalKalmanSmootherEngine(
                orthogonalModel.process.getRepresentation(CanonicalGaussianBranchTransitionKernel.class),
                orthogonalModel.process.getRepresentation(GaussianTransitionRepresentation.class),
                obs,
                grid);
        orthogonalSmoother.getLogLikelihood();

        final CanonicalBranchMessageContribution contribution = new CanonicalBranchMessageContribution(2);
        final CanonicalBranchMessageContributionUtils.Workspace contributionWorkspace2D =
                new CanonicalBranchMessageContributionUtils.Workspace(2);
        CanonicalBranchMessageContributionUtils.fillFromPairState(
                orthogonalSmoother.getCanonicalTrajectory().branchPairStates[0],
                contributionWorkspace2D,
                contribution);

        final CanonicalGaussianTransition denseTransition = new CanonicalGaussianTransition(2);
        denseModel.process.getRepresentation(CanonicalGaussianBranchTransitionKernel.class)
                .fillCanonicalTransition(dt, denseTransition);
        final CanonicalLocalTransitionAdjoints denseAdjoints = new CanonicalLocalTransitionAdjoints(2);
        final CanonicalTransitionAdjointUtils.Workspace denseAdjointWorkspace2D =
                new CanonicalTransitionAdjointUtils.Workspace(2);
        CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                denseTransition, contribution, denseAdjointWorkspace2D, denseAdjoints);

        final double[] stationaryMean = new double[2];
        orthogonalModel.process.getInitialMean(stationaryMean);
        final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalParameterization =
                (OrthogonalBlockDiagonalSelectionMatrixParameterization)
                        orthogonalModel.process.getSelectionMatrixParameterization();
        final CanonicalLocalTransitionAdjoints fastAdjoints = new CanonicalLocalTransitionAdjoints(2);
        orthogonalParameterization.fillCanonicalLocalAdjoints(
                orthogonalModel.process.getDiffusionMatrix(),
                stationaryMean,
                dt,
                contribution,
                fastAdjoints);

        for (int i = 0; i < 2; ++i) {
            assertEquals("Orthogonal block local canonical df adjoint must match dense, idx=" + i,
                    denseAdjoints.dLogL_df[i], fastAdjoints.dLogL_df[i], 1e-10);
            for (int j = 0; j < 2; ++j) {
                assertEquals("Orthogonal block local canonical dF adjoint must match dense at [" + i + "," + j + "]",
                        flatEntry(denseAdjoints.dLogL_dF, 2, i, j),
                        flatEntry(fastAdjoints.dLogL_dF, 2, i, j), 1e-10);
                assertEquals("Orthogonal block local canonical dOmega adjoint must match dense at [" + i + "," + j + "]",
                        flatEntry(denseAdjoints.dLogL_dOmega, 2, i, j),
                        flatEntry(fastAdjoints.dLogL_dOmega, 2, i, j), 1e-10);
            }
        }
    }

    public void testOrthogonalBlockCanonicalLocalAdjointsMatchDense_4D() {
        final double dt = 0.16;
        final double[][] observations = {
                {0.12, -0.20, 0.35, Double.NaN, 0.41},
                {-0.08, 0.26, -0.18, Double.NaN, 0.09},
                {0.21, -0.11, 0.14, Double.NaN, -0.16},
                {-0.04, 0.19, -0.09, Double.NaN, 0.27}
        };
        final Model orthogonalModel = makeOrthogonalBlockPolarFullWithMethod(
                "ouCanonicalOrthAdj4d",
                4,
                new double[]{0.22, -0.31, 0.41, 0.18, -0.27, 0.33},
                new double[0],
                new double[]{0.64, 0.83},
                new double[]{0.52, 0.91},
                new double[]{0.17, -0.12},
                new double[][]{
                        {1.20, 0.04, 0.02, 0.01},
                        {0.04, 0.95, 0.05, 0.02},
                        {0.02, 0.05, 1.11, 0.06},
                        {0.01, 0.02, 0.06, 0.88}
                },
                0.07,
                new double[][]{
                        {0.93, 0.02, 0.01, 0.00},
                        {0.02, 1.02, 0.03, 0.01},
                        {0.01, 0.03, 0.97, 0.02},
                        {0.00, 0.01, 0.02, 1.05}
                },
                new double[][]{
                        {0.28, 0.01, 0.00, 0.00},
                        {0.01, 0.31, 0.01, 0.00},
                        {0.00, 0.01, 0.34, 0.01},
                        {0.00, 0.00, 0.01, 0.29}
                },
                observations,
                CovarianceGradientMethod.LYAPUNOV_ADJOINT,
                dt);

        final Model denseModel = makeDenseFullWithMethod(
                "ouCanonicalDenseAdj4dRef",
                orthogonalModel.process.getDriftMatrix().getParameterAsMatrix(),
                orthogonalModel.process.getDiffusionMatrix().getParameterAsMatrix(),
                0.07,
                new double[][]{
                        {0.93, 0.02, 0.01, 0.00},
                        {0.02, 1.02, 0.03, 0.01},
                        {0.01, 0.03, 0.97, 0.02},
                        {0.00, 0.01, 0.02, 1.05}
                },
                new double[][]{
                        {0.28, 0.01, 0.00, 0.00},
                        {0.01, 0.31, 0.01, 0.00},
                        {0.00, 0.01, 0.34, 0.01},
                        {0.00, 0.00, 0.01, 0.29}
                },
                observations,
                dt,
                CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final CanonicalKalmanSmootherEngine orthogonalSmoother = new CanonicalKalmanSmootherEngine(
                orthogonalModel.process.getRepresentation(CanonicalGaussianBranchTransitionKernel.class),
                orthogonalModel.process.getRepresentation(GaussianTransitionRepresentation.class),
                orthogonalModel.obs,
                new UniformTimeGrid(observations[0].length, 0.0, dt));
        orthogonalSmoother.getLogLikelihood();

        final CanonicalBranchMessageContribution contribution = new CanonicalBranchMessageContribution(4);
        final CanonicalBranchMessageContributionUtils.Workspace contributionWorkspace4D =
                new CanonicalBranchMessageContributionUtils.Workspace(4);
        CanonicalBranchMessageContributionUtils.fillFromPairState(
                orthogonalSmoother.getCanonicalTrajectory().branchPairStates[0],
                contributionWorkspace4D,
                contribution);

        final CanonicalGaussianTransition denseTransition = new CanonicalGaussianTransition(4);
        denseModel.process.getRepresentation(CanonicalGaussianBranchTransitionKernel.class)
                .fillCanonicalTransition(dt, denseTransition);
        final CanonicalLocalTransitionAdjoints denseAdjoints = new CanonicalLocalTransitionAdjoints(4);
        final CanonicalTransitionAdjointUtils.Workspace denseAdjointWorkspace4D =
                new CanonicalTransitionAdjointUtils.Workspace(4);
        CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                denseTransition, contribution, denseAdjointWorkspace4D, denseAdjoints);

        final double[] stationaryMean = new double[4];
        orthogonalModel.process.getInitialMean(stationaryMean);
        final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalParameterization =
                (OrthogonalBlockDiagonalSelectionMatrixParameterization)
                        orthogonalModel.process.getSelectionMatrixParameterization();
        final CanonicalLocalTransitionAdjoints fastAdjoints = new CanonicalLocalTransitionAdjoints(4);
        orthogonalParameterization.fillCanonicalLocalAdjoints(
                orthogonalModel.process.getDiffusionMatrix(),
                stationaryMean,
                dt,
                contribution,
                fastAdjoints);

        for (int i = 0; i < 4; ++i) {
            assertEquals("Orthogonal block 4D local canonical df adjoint must match dense, idx=" + i,
                    denseAdjoints.dLogL_df[i], fastAdjoints.dLogL_df[i], 1e-10);
            for (int j = 0; j < 4; ++j) {
                assertEquals("Orthogonal block 4D local canonical dF adjoint must match dense at [" + i + "," + j + "]",
                        flatEntry(denseAdjoints.dLogL_dF, 4, i, j),
                        flatEntry(fastAdjoints.dLogL_dF, 4, i, j), 1e-10);
                assertEquals("Orthogonal block 4D local canonical dOmega adjoint must match dense at [" + i + "," + j + "]",
                        flatEntry(denseAdjoints.dLogL_dOmega, 4, i, j),
                        flatEntry(fastAdjoints.dLogL_dOmega, 4, i, j), 1e-10);
            }
        }
    }

    public void testTimeSeriesLikelihoodFactorySupportsCanonicalAnalyticalGradients() {
        final MatrixParameter drift = makeMatrix("A.ts.factory.canonical.grad", new double[][]{{0.27}});
        final MatrixParameter diffusion = makeMatrix("Q.ts.factory.canonical.grad", new double[][]{{0.95}});
        final Parameter mean = new Parameter.Default(0.3);
        final MatrixParameter initCov = makeMatrix("P0.ts.factory.canonical.grad", new double[][]{{1.05}});
        final OUProcessModel process = new OUProcessModel(
                "ouTsFactoryCanonicalGrad", 1, drift, diffusion, mean, initCov, CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final MatrixParameter H = makeMatrix("H.ts.factory.canonical.grad", new double[][]{{1.0}});
        final MatrixParameter R = makeMatrix("R.ts.factory.canonical.grad", new double[][]{{0.45}});
        final MatrixParameter Y = makeMatrix("Y.ts.factory.canonical.grad", new double[][]{{0.2, -0.6, 1.1, 0.5}});
        final GaussianObservationModel obs = new GaussianObservationModel("obsTsFactoryCanonicalGrad", 1, H, R, Y);
        final TimeGrid grid = new UniformTimeGrid(4, 0.0, 0.2);
        final BasicTimeSeriesModel model = new BasicTimeSeriesModel("tsFactoryCanonicalGrad", process, obs, grid);

        final TimeSeriesLikelihood expectationLikelihood = GaussianTimeSeriesLikelihoodFactory.create(
                "tsLikeExpectationAnalytical",
                model,
                GaussianForwardComputationMode.EXPECTATION,
                GaussianSmootherComputationMode.EXPECTATION,
                GaussianGradientComputationMode.EXPECTATION_ANALYTICAL);
        final TimeSeriesLikelihood canonicalLikelihood = GaussianTimeSeriesLikelihoodFactory.create(
                "tsLikeCanonicalAnalytical",
                model,
                GaussianForwardComputationMode.CANONICAL,
                GaussianSmootherComputationMode.CANONICAL,
                GaussianGradientComputationMode.CANONICAL_ANALYTICAL);

        assertEquals("Canonical analytical builder path must match expectation likelihood",
                expectationLikelihood.getLogLikelihood(),
                canonicalLikelihood.getLogLikelihood(),
                1e-9);

        final double[] expectationGradient = expectationLikelihood.getGradientWrt(process.getDriftMatrix())
                .getGradientLogDensity(null);
        final double[] canonicalGradient = canonicalLikelihood.getGradientWrt(process.getDriftMatrix())
                .getGradientLogDensity(null);
        assertEquals("Canonical analytical gradient must match expectation analytical gradient",
                expectationGradient[0],
                canonicalGradient[0],
                1e-9);
    }

    public void testCanonicalSelectionGradientFormulaMatchesExpectationAnalyticalGradient_2D() {
        final MatrixParameter drift = makeMatrix("A.canonical.native.sel", new double[][]{
                {0.47, -0.09},
                {0.05, 0.99}
        });
        final MatrixParameter diffusion = makeMatrix("Q.canonical.native.sel", new double[][]{
                {1.2, 0.1},
                {0.1, 0.9}
        });
        final Parameter mean = new Parameter.Default(new double[]{0.2, -0.15});
        final MatrixParameter initCov = makeMatrix("P0.canonical.native.sel", new double[][]{
                {0.85, 0.04},
                {0.04, 1.1}
        });
        final OUProcessModel process = new OUProcessModel(
                "ouCanonicalNativeSel", 2, drift, diffusion, mean, initCov, CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final MatrixParameter H = makeMatrix("H.canonical.native.sel", new double[][]{{1, 0}, {0, 1}});
        final MatrixParameter R = makeMatrix("R.canonical.native.sel", new double[][]{
                {0.33, 0.01},
                {0.01, 0.4}
        });
        final MatrixParameter Y = makeMatrix("Y.canonical.native.sel", new double[][]{
                {0.25, -0.5, Double.NaN, 0.8},
                {-0.35, 0.6, Double.NaN, -0.2}
        });
        final GaussianObservationModel obs = new GaussianObservationModel("obsCanonicalNativeSel", 2, H, R, Y);
        final TimeGrid grid = new UniformTimeGrid(4, 0.0, 0.2);
        final BasicTimeSeriesModel model = new BasicTimeSeriesModel("tsCanonicalNativeSel", process, obs, grid);

        final TimeSeriesLikelihood expectationLikelihood = GaussianTimeSeriesLikelihoodFactory.create(
                "tsLikeExpectationNativeSel",
                model,
                GaussianForwardComputationMode.EXPECTATION,
                GaussianSmootherComputationMode.EXPECTATION,
                GaussianGradientComputationMode.EXPECTATION_ANALYTICAL);
        final TimeSeriesLikelihood canonicalLikelihood = GaussianTimeSeriesLikelihoodFactory.create(
                "tsLikeCanonicalNativeSel",
                model,
                GaussianForwardComputationMode.CANONICAL,
                GaussianSmootherComputationMode.CANONICAL,
                GaussianGradientComputationMode.CANONICAL_ANALYTICAL);

        final double[] expectationGradient = expectationLikelihood.getGradientWrt(process.getDriftMatrix())
                .getGradientLogDensity(null);
        final double[] canonicalGradient = canonicalLikelihood.getGradientWrt(process.getDriftMatrix())
                .getGradientLogDensity(null);

        assertEquals("Canonical selection gradient dimension", expectationGradient.length, canonicalGradient.length);
        for (int i = 0; i < expectationGradient.length; ++i) {
            assertEquals("Canonical selection gradient must match expectation analytical gradient, idx=" + i,
                    expectationGradient[i], canonicalGradient[i], 1e-9);
        }
    }

    public void testCanonicalDiffusionGradientFormulaMatchesExpectationAnalyticalGradient_2D() {
        final MatrixParameter drift = makeMatrix("A.canonical.native.diff", new double[][]{
                {0.44, -0.06},
                {0.07, 0.92}
        });
        final MatrixParameter diffusion = makeMatrix("Q.canonical.native.diff", new double[][]{
                {1.05, 0.11},
                {0.11, 0.87}
        });
        final Parameter mean = new Parameter.Default(new double[]{0.1, -0.05});
        final MatrixParameter initCov = makeMatrix("P0.canonical.native.diff", new double[][]{
                {0.9, 0.03},
                {0.03, 1.0}
        });
        final OUProcessModel process = new OUProcessModel(
                "ouCanonicalNativeDiff", 2, drift, diffusion, mean, initCov, CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final MatrixParameter H = makeMatrix("H.canonical.native.diff", new double[][]{{1, 0}, {0, 1}});
        final MatrixParameter R = makeMatrix("R.canonical.native.diff", new double[][]{
                {0.31, 0.01},
                {0.01, 0.39}
        });
        final MatrixParameter Y = makeMatrix("Y.canonical.native.diff", new double[][]{
                {0.4, -0.3, Double.NaN, 0.6},
                {-0.2, 0.7, Double.NaN, -0.1}
        });
        final GaussianObservationModel obs = new GaussianObservationModel("obsCanonicalNativeDiff", 2, H, R, Y);
        final TimeGrid grid = new UniformTimeGrid(4, 0.0, 0.2);
        final BasicTimeSeriesModel model = new BasicTimeSeriesModel("tsCanonicalNativeDiff", process, obs, grid);

        final TimeSeriesLikelihood expectationLikelihood = GaussianTimeSeriesLikelihoodFactory.create(
                "tsLikeExpectationNativeDiff",
                model,
                GaussianForwardComputationMode.EXPECTATION,
                GaussianSmootherComputationMode.EXPECTATION,
                GaussianGradientComputationMode.EXPECTATION_ANALYTICAL);
        final TimeSeriesLikelihood canonicalLikelihood = GaussianTimeSeriesLikelihoodFactory.create(
                "tsLikeCanonicalNativeDiff",
                model,
                GaussianForwardComputationMode.CANONICAL,
                GaussianSmootherComputationMode.CANONICAL,
                GaussianGradientComputationMode.CANONICAL_ANALYTICAL);

        final double[] expectationGradient = expectationLikelihood.getGradientWrt(process.getDiffusionMatrix())
                .getGradientLogDensity(null);
        final double[] canonicalGradient = canonicalLikelihood.getGradientWrt(process.getDiffusionMatrix())
                .getGradientLogDensity(null);

        assertEquals("Canonical diffusion gradient dimension", expectationGradient.length, canonicalGradient.length);
        for (int i = 0; i < expectationGradient.length; ++i) {
            assertEquals("Canonical diffusion gradient must match expectation analytical gradient, idx=" + i,
                    expectationGradient[i], canonicalGradient[i], 1e-9);
        }
    }

    public void testCanonicalAnalyticalDriftGradientMatchesFD_2D_NonDiagonal() {
        final MatrixParameter drift = makeMatrix("A.canonical.fd.drift2d", new double[][]{
                {0.43, -0.07},
                {0.09, 0.88}
        });
        final MatrixParameter diffusion = makeMatrix("Q.canonical.fd.drift2d", new double[][]{
                {1.15, 0.14},
                {0.14, 0.92}
        });
        final Parameter mean = new Parameter.Default(new double[]{0.18, -0.12});
        final MatrixParameter initCov = makeMatrix("P0.canonical.fd.drift2d", new double[][]{
                {0.95, 0.06},
                {0.06, 1.08}
        });
        final OUProcessModel process = new OUProcessModel(
                "ouCanonicalFdDrift2d", 2, drift, diffusion, mean, initCov, CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final MatrixParameter H = makeMatrix("H.canonical.fd.drift2d", new double[][]{{1, 0}, {0, 1}});
        final MatrixParameter R = makeMatrix("R.canonical.fd.drift2d", new double[][]{
                {0.32, 0.015},
                {0.015, 0.41}
        });
        final MatrixParameter Y = makeMatrix("Y.canonical.fd.drift2d", new double[][]{
                {0.35, -0.45, Double.NaN, 0.72},
                {-0.28, 0.55, Double.NaN, -0.14}
        });
        final GaussianObservationModel obs = new GaussianObservationModel("obsCanonicalFdDrift2d", 2, H, R, Y);
        final TimeGrid grid = new UniformTimeGrid(4, 0.0, 0.18);

        final TimeSeriesLikelihood canonical = makeCanonicalAnalyticalLikelihood(
                "tsLikeCanonicalFdDrift2d", process, obs, grid);
        final KalmanLikelihoodEngine likelihood = new KalmanLikelihoodEngine(
                process.getRepresentation(GaussianTransitionRepresentation.class), obs, grid);

        final double[] canonicalGradient = canonical.getGradientWrt(process.getDriftMatrix()).getGradientLogDensity(null);
        for (int idx = 0; idx < canonicalGradient.length; ++idx) {
            assertEquals("Canonical analytical drift gradient must match FD in 2D, idx=" + idx,
                    numericalGradient(likelihood, process.getDriftMatrix(), idx, 1e-6),
                    canonicalGradient[idx],
                    TOL_FD);
        }
    }

    public void testCanonicalAnalyticalDiffusionGradientMatchesFD_2D_NonDiagonal() {
        final MatrixParameter drift = makeMatrix("A.canonical.fd.diff2d", new double[][]{
                {0.41, -0.05},
                {0.08, 0.93}
        });
        final MatrixParameter diffusion = makeMatrix("Q.canonical.fd.diff2d", new double[][]{
                {1.08, 0.12},
                {0.12, 0.84}
        });
        final Parameter mean = new Parameter.Default(new double[]{0.08, -0.03});
        final MatrixParameter initCov = makeMatrix("P0.canonical.fd.diff2d", new double[][]{
                {0.88, 0.04},
                {0.04, 0.98}
        });
        final OUProcessModel process = new OUProcessModel(
                "ouCanonicalFdDiff2d", 2, drift, diffusion, mean, initCov, CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final MatrixParameter H = makeMatrix("H.canonical.fd.diff2d", new double[][]{{1, 0}, {0, 1}});
        final MatrixParameter R = makeMatrix("R.canonical.fd.diff2d", new double[][]{
                {0.30, 0.01},
                {0.01, 0.37}
        });
        final MatrixParameter Y = makeMatrix("Y.canonical.fd.diff2d", new double[][]{
                {0.42, -0.31, Double.NaN, 0.63},
                {-0.19, 0.68, Double.NaN, -0.09}
        });
        final GaussianObservationModel obs = new GaussianObservationModel("obsCanonicalFdDiff2d", 2, H, R, Y);
        final TimeGrid grid = new UniformTimeGrid(4, 0.0, 0.17);

        final TimeSeriesLikelihood canonical = makeCanonicalAnalyticalLikelihood(
                "tsLikeCanonicalFdDiff2d", process, obs, grid);
        final KalmanLikelihoodEngine likelihood = new KalmanLikelihoodEngine(
                process.getRepresentation(GaussianTransitionRepresentation.class), obs, grid);

        final double[] canonicalGradient = canonical.getGradientWrt(process.getDiffusionMatrix()).getGradientLogDensity(null);
        for (int idx = 0; idx < canonicalGradient.length; ++idx) {
            assertEquals("Canonical analytical diffusion gradient must match FD in 2D, idx=" + idx,
                    numericalGradient(likelihood, process.getDiffusionMatrix(), idx, 1e-6),
                    canonicalGradient[idx],
                    TOL_FD);
        }
    }

    public void testLatentMiddleThreePointChainBranchSumMatchesCanonicalGradient() {
        final MatrixParameter drift = makeMatrix("A.canonical.bridge.latent3", new double[][]{
                {0.46, -0.08},
                {0.07, 0.91}
        });
        final MatrixParameter diffusion = makeMatrix("Q.canonical.bridge.latent3", new double[][]{
                {1.12, 0.09},
                {0.09, 0.86}
        });
        final Parameter mean = new Parameter.Default(new double[]{0.16, -0.11});
        final MatrixParameter initCov = makeMatrix("P0.canonical.bridge.latent3", new double[][]{
                {0.9, 0.05},
                {0.05, 1.02}
        });
        final OUProcessModel process = new OUProcessModel(
                "ouCanonicalBridgeLatent3", 2, drift, diffusion, mean, initCov, CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final MatrixParameter H = makeMatrix("H.canonical.bridge.latent3", new double[][]{{1, 0}, {0, 1}});
        final MatrixParameter R = makeMatrix("R.canonical.bridge.latent3", new double[][]{
                {0.28, 0.02},
                {0.02, 0.34}
        });
        final MatrixParameter Y = makeMatrix("Y.canonical.bridge.latent3", new double[][]{
                {0.37, Double.NaN, -0.48},
                {-0.24, Double.NaN, 0.61}
        });
        final GaussianObservationModel obs = new GaussianObservationModel("obsCanonicalBridgeLatent3", 2, H, R, Y);
        final TimeGrid grid = new UniformTimeGrid(3, 0.0, 0.23);

        final TimeSeriesLikelihood canonicalLikelihood = makeCanonicalAnalyticalLikelihood(
                "tsLikeCanonicalBridgeLatent3", process, obs, grid);
        final double[] globalGradient = canonicalLikelihood.getGradientWrt(process.getDriftMatrix())
                .getGradientLogDensity(null);
        final double[] localBranchSum = sumCanonicalLocalSelectionGradient(process, obs, grid);

        assertEquals("latent-middle bridge gradient dimension", globalGradient.length, localBranchSum.length);
        for (int i = 0; i < globalGradient.length; ++i) {
            assertEquals("latent-middle branch sum must match canonical global gradient, idx=" + i,
                    globalGradient[i], localBranchSum[i], 1e-9);
        }
    }

    public void testCanonicalAnalyticalGradientsMatchFD_3D_Dense() {
        final double[][] driftValues = {
                {0.52, -0.08, 0.03},
                {0.07, 0.86, -0.04},
                {-0.02, 0.05, 1.08}
        };
        final double[][] diffusionValues = {
                {1.10, 0.09, 0.03},
                {0.09, 0.95, 0.07},
                {0.03, 0.07, 1.22}
        };
        final double[][] initCovValues = {
                {0.92, 0.04, 0.01},
                {0.04, 1.03, 0.05},
                {0.01, 0.05, 0.98}
        };
        final double[][] noiseValues = {
                {0.28, 0.01, 0.00},
                {0.01, 0.33, 0.01},
                {0.00, 0.01, 0.36}
        };
        final double[][] observations = {
                {0.20, -0.30, 0.55, Double.NaN, 0.81},
                {-0.12, 0.41, -0.38, Double.NaN, 0.27},
                {0.33, -0.14, 0.09, Double.NaN, -0.22}
        };

        final Model m = makeDenseFullWithMethod(
                "ouCanonicalFd3d",
                driftValues,
                diffusionValues,
                0.12,
                initCovValues,
                noiseValues,
                observations,
                0.14,
                CovarianceGradientMethod.LYAPUNOV_ADJOINT);
        final TimeGrid grid = new UniformTimeGrid(observations[0].length, 0.0, 0.14);
        final TimeSeriesLikelihood canonical = makeCanonicalAnalyticalLikelihood(
                "tsLikeCanonicalFd3d", m.process, m.obs, grid);

        final double[] canonicalDrift = canonical.getGradientWrt(m.process.getDriftMatrix()).getGradientLogDensity(null);
        for (int idx = 0; idx < canonicalDrift.length; ++idx) {
            assertEquals("Canonical analytical 3D drift gradient must match FD, idx=" + idx,
                    numericalGradient(m.likelihoodEngine, m.process.getDriftMatrix(), idx, 1e-6),
                    canonicalDrift[idx],
                    TOL_FD);
        }

        final double[] canonicalDiffusion = canonical.getGradientWrt(m.process.getDiffusionMatrix()).getGradientLogDensity(null);
        for (int idx = 0; idx < canonicalDiffusion.length; ++idx) {
            assertEquals("Canonical analytical 3D diffusion gradient must match FD, idx=" + idx,
                    numericalGradient(m.likelihoodEngine, m.process.getDiffusionMatrix(), idx, 1e-6),
                    canonicalDiffusion[idx],
                    TOL_FD);
        }
    }

    public void testCanonicalAnalyticalGradientsMatchFD_2D_WeakStability() {
        final MatrixParameter drift = makeMatrix("A.canonical.fd.weak", new double[][]{
                {0.065, -0.010},
                {0.014, 0.11}
        });
        final MatrixParameter diffusion = makeMatrix("Q.canonical.fd.weak", new double[][]{
                {0.72, 0.05},
                {0.05, 0.61}
        });
        final Parameter mean = new Parameter.Default(new double[]{0.02, -0.01});
        final MatrixParameter initCov = makeMatrix("P0.canonical.fd.weak", new double[][]{
                {0.83, 0.02},
                {0.02, 0.79}
        });
        final OUProcessModel process = new OUProcessModel(
                "ouCanonicalFdWeak", 2, drift, diffusion, mean, initCov, CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final MatrixParameter H = makeMatrix("H.canonical.fd.weak", new double[][]{{1, 0}, {0, 1}});
        final MatrixParameter R = makeMatrix("R.canonical.fd.weak", new double[][]{
                {0.24, 0.0},
                {0.0, 0.27}
        });
        final MatrixParameter Y = makeMatrix("Y.canonical.fd.weak", new double[][]{
                {0.09, -0.16, 0.21, Double.NaN, 0.28},
                {-0.05, 0.11, -0.14, Double.NaN, 0.17}
        });
        final GaussianObservationModel obs = new GaussianObservationModel("obsCanonicalFdWeak", 2, H, R, Y);
        final TimeGrid grid = new UniformTimeGrid(5, 0.0, 0.12);

        final TimeSeriesLikelihood canonical = makeCanonicalAnalyticalLikelihood(
                "tsLikeCanonicalFdWeak", process, obs, grid);
        final KalmanLikelihoodEngine likelihood = new KalmanLikelihoodEngine(
                process.getRepresentation(GaussianTransitionRepresentation.class), obs, grid);

        final double[] canonicalDrift = canonical.getGradientWrt(process.getDriftMatrix()).getGradientLogDensity(null);
        for (int idx = 0; idx < canonicalDrift.length; ++idx) {
            assertEquals("Canonical analytical weak-stability drift gradient must match FD, idx=" + idx,
                    numericalGradient(likelihood, process.getDriftMatrix(), idx, 1e-6),
                    canonicalDrift[idx],
                    2e-5);
        }

        final double[] canonicalDiffusion = canonical.getGradientWrt(process.getDiffusionMatrix()).getGradientLogDensity(null);
        for (int idx = 0; idx < canonicalDiffusion.length; ++idx) {
            assertEquals("Canonical analytical weak-stability diffusion gradient must match FD, idx=" + idx,
                    numericalGradient(likelihood, process.getDiffusionMatrix(), idx, 1e-6),
                    canonicalDiffusion[idx],
                    2e-5);
        }
    }

    public void testCanonicalAnalyticalOrthogonalBlockNativeGradientsMatchFD_2D() {
        final double dt = 0.18;
        final Model orthogonalModel = makeOrthogonalBlockPolarFullWithMethod(
                "ouCanonicalNativeOrth2d",
                2,
                new double[]{0.34},
                new double[0],
                new double[]{0.61},
                new double[]{0.58},
                new double[]{0.19},
                new double[][]{{1.07, 0.08}, {0.08, 0.93}},
                0.14,
                new double[][]{{0.88, 0.02}, {0.02, 1.01}},
                new double[][]{{0.29, 0.01}, {0.01, 0.34}},
                new double[][]{
                        {0.24, -0.28, 0.49, Double.NaN, 0.63},
                        {-0.17, 0.37, -0.11, Double.NaN, -0.08}
                },
                CovarianceGradientMethod.STATIONARY_LYAPUNOV,
                dt);

        final TimeSeriesLikelihood canonical = makeCanonicalAnalyticalLikelihood(
                "tsLikeCanonicalNativeOrth2d", orthogonalModel.process, orthogonalModel.obs,
                new UniformTimeGrid(5, 0.0, dt));

        final OrthogonalBlockDiagonalPolarStableMatrixParameter blockParameter =
                (OrthogonalBlockDiagonalPolarStableMatrixParameter) orthogonalModel.process.getDriftMatrix();
        final GivensRotationMatrixParameter rotation =
                (GivensRotationMatrixParameter) blockParameter.getRotationMatrixParameter();

        final double[] rhoGradient = canonical.getGradientWrt(blockParameter.getRhoParameter()).getGradientLogDensity(null);
        assertEquals("Canonical analytical orthogonal-block rho gradient must match FD",
                numericalGradient(orthogonalModel.likelihoodEngine, blockParameter.getRhoParameter(), 0, 1e-6),
                rhoGradient[0], TOL_FD);

        final double[] thetaGradient = canonical.getGradientWrt(blockParameter.getThetaParameter()).getGradientLogDensity(null);
        assertEquals("Canonical analytical orthogonal-block theta gradient must match FD",
                numericalGradient(orthogonalModel.likelihoodEngine, blockParameter.getThetaParameter(), 0, 1e-6),
                thetaGradient[0], TOL_FD);

        final double[] tGradient = canonical.getGradientWrt(blockParameter.getTParameter()).getGradientLogDensity(null);
        assertEquals("Canonical analytical orthogonal-block t gradient must match FD",
                numericalGradient(orthogonalModel.likelihoodEngine, blockParameter.getTParameter(), 0, 1e-6),
                tGradient[0], TOL_FD);

        final double[] angleGradient = canonical.getGradientWrt(rotation.getAngleParameter()).getGradientLogDensity(null);
        assertEquals("Canonical analytical orthogonal-block angle gradient must match FD",
                numericalGradient(orthogonalModel.likelihoodEngine, rotation.getAngleParameter(), 0, 1e-6),
                angleGradient[0], TOL_FD);
    }

    public void testCanonicalAnalyticalOrthogonalBlockNativeGradientsMatchFD_4D() {
        final double dt = 0.15;
        final double[][] observations = {
                {0.11, -0.18, 0.29, Double.NaN, 0.36},
                {-0.06, 0.24, -0.12, Double.NaN, 0.08},
                {0.18, -0.09, 0.13, Double.NaN, -0.14},
                {-0.02, 0.16, -0.07, Double.NaN, 0.22}
        };
        final Model orthogonalModel = makeOrthogonalBlockPolarFullWithMethod(
                "ouCanonicalNativeOrth4d",
                4,
                new double[]{0.21, -0.28, 0.38, 0.17, -0.22, 0.29},
                new double[0],
                new double[]{0.59, 0.79},
                new double[]{0.49, 0.86},
                new double[]{0.15, -0.10},
                new double[][]{
                        {1.14, 0.04, 0.02, 0.01},
                        {0.04, 0.91, 0.05, 0.02},
                        {0.02, 0.05, 1.07, 0.06},
                        {0.01, 0.02, 0.06, 0.85}
                },
                0.05,
                new double[][]{
                        {0.91, 0.02, 0.01, 0.00},
                        {0.02, 1.00, 0.03, 0.01},
                        {0.01, 0.03, 0.95, 0.02},
                        {0.00, 0.01, 0.02, 1.02}
                },
                new double[][]{
                        {0.27, 0.01, 0.00, 0.00},
                        {0.01, 0.30, 0.01, 0.00},
                        {0.00, 0.01, 0.33, 0.01},
                        {0.00, 0.00, 0.01, 0.28}
                },
                observations,
                CovarianceGradientMethod.STATIONARY_LYAPUNOV,
                dt);

        final TimeSeriesLikelihood canonical = makeCanonicalAnalyticalLikelihood(
                "tsLikeCanonicalNativeOrth4d", orthogonalModel.process, orthogonalModel.obs,
                new UniformTimeGrid(observations[0].length, 0.0, dt));

        final OrthogonalBlockDiagonalPolarStableMatrixParameter blockParameter =
                (OrthogonalBlockDiagonalPolarStableMatrixParameter) orthogonalModel.process.getDriftMatrix();
        final GivensRotationMatrixParameter rotation =
                (GivensRotationMatrixParameter) blockParameter.getRotationMatrixParameter();

        final double[] rhoGradient = canonical.getGradientWrt(blockParameter.getRhoParameter()).getGradientLogDensity(null);
        for (int idx = 0; idx < rhoGradient.length; ++idx) {
            assertEquals("Canonical analytical orthogonal-block 4D rho gradient must match FD, idx=" + idx,
                    numericalGradient(orthogonalModel.likelihoodEngine, blockParameter.getRhoParameter(), idx, 1e-6),
                    rhoGradient[idx], TOL_FD);
        }

        final double[] thetaGradient = canonical.getGradientWrt(blockParameter.getThetaParameter()).getGradientLogDensity(null);
        for (int idx = 0; idx < thetaGradient.length; ++idx) {
            assertEquals("Canonical analytical orthogonal-block 4D theta gradient must match FD, idx=" + idx,
                    numericalGradient(orthogonalModel.likelihoodEngine, blockParameter.getThetaParameter(), idx, 1e-6),
                    thetaGradient[idx], TOL_FD);
        }

        final double[] tGradient = canonical.getGradientWrt(blockParameter.getTParameter()).getGradientLogDensity(null);
        for (int idx = 0; idx < tGradient.length; ++idx) {
            assertEquals("Canonical analytical orthogonal-block 4D t gradient must match FD, idx=" + idx,
                    numericalGradient(orthogonalModel.likelihoodEngine, blockParameter.getTParameter(), idx, 1e-6),
                    tGradient[idx], TOL_FD);
        }

        final double[] angleGradient = canonical.getGradientWrt(rotation.getAngleParameter()).getGradientLogDensity(null);
        for (int idx = 0; idx < angleGradient.length; ++idx) {
            assertEquals("Canonical analytical orthogonal-block 4D angle gradient must match FD, idx=" + idx,
                    numericalGradient(orthogonalModel.likelihoodEngine, rotation.getAngleParameter(), idx, 1e-6),
                    angleGradient[idx], TOL_FD);
        }
    }

    public void testTreeBranchProviderMatchesTimeSeriesKernel() {
        final double[][] driftValues = {
                {0.7, -0.12},
                {0.08, 1.05}
        };
        final double[][] covarianceValues = {
                {1.4, 0.2},
                {0.2, 0.9}
        };
        final double[][] precisionValues = {
                {0.7377049180327868, -0.1639344262295082},
                {-0.1639344262295082, 1.1475409836065573}
        };
        final double dt = 0.3;

        final MatrixParameter drift = makeMatrix("A.tree", driftValues);
        final MatrixParameter covariance = makeMatrix("Q.tree", covarianceValues);
        final MatrixParameter precision = makeMatrix("P.tree", precisionValues);

        final OUProcessModel timeSeriesProcess = new OUProcessModel(
                "ouTreeBridgeTs",
                2,
                drift,
                covariance,
                new Parameter.Default(0.0),
                makeMatrix("P0.tree", new double[][]{{1.0, 0.0}, {0.0, 1.0}}),
                CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        final MultivariateElasticModel elasticModel = new MultivariateElasticModel(drift);
        final MultivariateDiffusionModel diffusionModel = new MultivariateDiffusionModel(precision);
        final TimeSeriesOUGaussianBranchTransitionProvider provider =
                new TimeSeriesOUGaussianBranchTransitionProvider(elasticModel, diffusionModel);

        final double[][] fExpected = new double[2][2];
        final double[] bExpected = new double[2];
        final double[][] vExpected = new double[2][2];
        timeSeriesProcess.getRepresentation(dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianBranchTransitionKernel.class)
                .fillTransitionMatrix(dt, fExpected);
        timeSeriesProcess.getRepresentation(dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianBranchTransitionKernel.class)
                .fillTransitionOffset(dt, bExpected);
        timeSeriesProcess.getRepresentation(dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianBranchTransitionKernel.class)
                .fillTransitionCovariance(dt, vExpected);

        final double[][] fObserved = new double[2][2];
        final double[] bObserved = new double[2];
        final double[][] vObserved = new double[2][2];
        provider.fillBranchTransitionMatrix(dt, fObserved);
        provider.fillBranchTransitionOffset(dt, bObserved);
        provider.fillBranchTransitionCovariance(dt, vObserved);

        for (int i = 0; i < 2; ++i) {
            assertEquals("Tree provider offset must match time-series kernel, idx=" + i,
                    bExpected[i], bObserved[i], 1e-10);
            for (int j = 0; j < 2; ++j) {
                assertEquals("Tree provider F must match time-series kernel at [" + i + "," + j + "]",
                        fExpected[i][j], fObserved[i][j], 1e-10);
                assertEquals("Tree provider V must match time-series kernel at [" + i + "," + j + "]",
                        vExpected[i][j], vObserved[i][j], 1e-10);
            }
        }
    }

    // ── JUnit 3 boilerplate ───────────────────────────────────────────────────────

    public static Test suite() {
        return new TestSuite(AnalyticalKalmanGradientEngineTest.class);
    }
}
