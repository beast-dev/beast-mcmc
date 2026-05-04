package test.dr.inference.timeseries;

import static test.dr.inference.timeseries.OUTimeSeriesTestSupport.representation;
import static test.dr.inference.timeseries.OUTimeSeriesTestSupport.representable;
import static test.dr.inference.timeseries.OUTimeSeriesTestSupport.latent;
import static test.dr.inference.timeseries.OUTimeSeriesTestSupport.supportsRepresentation;
import static test.dr.inference.timeseries.OUTimeSeriesTestSupport.getTransitionMatrix;
import static test.dr.inference.timeseries.OUTimeSeriesTestSupport.getTransitionOffset;
import static test.dr.inference.timeseries.OUTimeSeriesTestSupport.getTransitionCovariance;

import dr.inference.model.BlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.core.UniformTimeGrid;
import dr.inference.timeseries.engine.gaussian.AnalyticalKalmanGradientEngine;
import dr.inference.timeseries.engine.gaussian.KalmanLikelihoodEngine;
import dr.inference.timeseries.engine.gaussian.KalmanSmootherEngine;
import dr.inference.timeseries.engine.gaussian.SelectionMatrixGradientFormula;
import dr.inference.timeseries.gaussian.GaussianObservationModel;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

public class BoundarySweepBlockGradientCheck {

    private static final double FD_H = 1e-6;
    private static final double WARN_ABS = 2e-4;

    private static MatrixParameter makeMatrix(String name, double[][] values) {
        MatrixParameter mp = new MatrixParameter(name, values.length, values[0].length);
        for (int i = 0; i < values.length; i++) {
            for (int j = 0; j < values[0].length; j++) {
                mp.setParameterValue(i, j, values[i][j]);
            }
        }
        return mp;
    }

    private static final class Model {
        final BlockDiagonalPolarStableMatrixParameter drift;
        final KalmanLikelihoodEngine likelihood;
        final AnalyticalKalmanGradientEngine analytical;

        private Model(BlockDiagonalPolarStableMatrixParameter drift,
                      KalmanLikelihoodEngine likelihood,
                      AnalyticalKalmanGradientEngine analytical) {
            this.drift = drift;
            this.likelihood = likelihood;
            this.analytical = analytical;
        }
    }

    private static Model makeModel(double[][] rotationValues,
                                   double rho,
                                   double theta,
                                   double t) {
        final double[][] diffusion = new double[][]{{1.2, 0.1}, {0.1, 1.6}};
        final double[][] initCov = new double[][]{{1.0, 0.2}, {0.2, 1.3}};
        final double[][] noise = new double[][]{{0.6, 0.05}, {0.05, 0.8}};
        final double[] y1 = {0.7, -0.4, 1.1, 0.2};
        final double[] y2 = {-0.3, 0.8, -1.2, 0.5};
        final double dt = 0.15;

        MatrixParameter rotation = makeMatrix("R", rotationValues);
        Parameter scalarBlock = new Parameter.Default(0);
        Parameter rhoParam = new Parameter.Default(rho);
        Parameter thetaParam = new Parameter.Default(theta);
        Parameter tParam = new Parameter.Default(t);
        BlockDiagonalPolarStableMatrixParameter drift = new BlockDiagonalPolarStableMatrixParameter(
                "Ablock", rotation, scalarBlock, rhoParam, thetaParam, tParam);
        MatrixParameter diffusionParam = makeMatrix("Q", diffusion);
        Parameter mean = new Parameter.Default(0.25);
        MatrixParameter initCovParam = makeMatrix("P0", initCov);
        OUProcessModel process = new OUProcessModel("ou2dblock", 2, drift, diffusionParam, mean, initCovParam,
                OUProcessModel.CovarianceGradientMethod.STATIONARY_LYAPUNOV);

        MatrixParameter H = makeMatrix("H", new double[][]{{1, 0}, {0, 1}});
        MatrixParameter R = makeMatrix("Robs", noise);
        MatrixParameter Y = makeMatrix("Y", new double[][]{y1, y2});
        GaussianObservationModel obs = new GaussianObservationModel("obs", 2, H, R, Y);
        TimeGrid grid = new UniformTimeGrid(y1.length, 0.0, dt);
        GaussianTransitionRepresentation rep = representation(process, GaussianTransitionRepresentation.class);
        KalmanSmootherEngine smoother = new KalmanSmootherEngine(rep, obs, grid);
        AnalyticalKalmanGradientEngine analytical = new AnalyticalKalmanGradientEngine(
                smoother, new SelectionMatrixGradientFormula(process.getDriftMatrix(), 2));
        KalmanLikelihoodEngine likelihood = new KalmanLikelihoodEngine(rep, obs, grid);
        return new Model(drift, likelihood, analytical);
    }

    private static double numericalGradient(KalmanLikelihoodEngine engine,
                                            Parameter parameter,
                                            int index) {
        double orig = parameter.getParameterValue(index);
        try {
            parameter.setParameterValue(index, orig + FD_H);
            engine.makeDirty();
            double plus = engine.getLogLikelihood();
            parameter.setParameterValue(index, orig - FD_H);
            engine.makeDirty();
            double minus = engine.getLogLikelihood();
            parameter.setParameterValue(index, orig);
            engine.makeDirty();
            return (plus - minus) / (2.0 * FD_H);
        } finally {
            parameter.setParameterValue(index, orig);
            engine.makeDirty();
        }
    }

    public static void main(String[] args) {
        final double[] rhoValues = {1e-4, 1e-3, 0.02, 0.9, 3.0};
        final double[] thetaValues = {0.0, 1.2, 1.45, 1.55, 1.569};
        final double[] tValues = {0.0, 1e-6, 1e-3, 0.2, 1.0};
        final double[][][] rotations = {
                {{1.0, 0.25}, {-0.15, 0.9}},
                {{1.0, 0.98}, {1.0, 0.981}},
                {{1.0, 0.999}, {1.0, 0.9995}}
        };
        final String[] rotationLabels = {
                "well-conditioned",
                "mildly-ill-conditioned",
                "near-singular"
        };

        int cases = 0;
        int warnings = 0;
        int invalid = 0;
        double maxAbs = -1.0;
        String worst = "";
        int filteredCases = 0;
        int filteredWarnings = 0;
        double filteredMaxAbs = -1.0;
        String filteredWorst = "";

        for (int r = 0; r < rotations.length; ++r) {
            for (double rho : rhoValues) {
                for (double theta : thetaValues) {
                    if (rho * Math.cos(theta) <= 5e-5) {
                        continue;
                    }
                    for (double t : tValues) {
                        Model model = makeModel(rotations[r], rho, theta, t);
                        BlockDiagonalPolarStableMatrixParameter drift = model.drift;
                        try {
                            model.likelihood.getLogLikelihood();
                        } catch (RuntimeException ex) {
                            invalid++;
                            continue;
                        }

                        final double[] analyticRho;
                        final double[] analyticTheta;
                        final double[] analyticT;
                        final double[] analyticR;
                        try {
                            analyticRho = model.analytical.getGradientWrt(drift.getRhoParameter());
                            analyticTheta = model.analytical.getGradientWrt(drift.getThetaParameter());
                            analyticT = model.analytical.getGradientWrt(drift.getTParameter());
                            analyticR = model.analytical.getGradientWrt(drift.getRotationMatrixParameter());
                        } catch (RuntimeException ex) {
                            invalid++;
                            continue;
                        }

                        final double fdRho;
                        final double fdTheta;
                        final double fdT;
                        try {
                            fdRho = numericalGradient(model.likelihood, drift.getRhoParameter(), 0);
                            fdTheta = numericalGradient(model.likelihood, drift.getThetaParameter(), 0);
                            fdT = numericalGradient(model.likelihood, drift.getTParameter(), 0);
                        } catch (RuntimeException ex) {
                            invalid++;
                            continue;
                        }

                        double errRho = Math.abs(analyticRho[0] - fdRho);
                        double errTheta = Math.abs(analyticTheta[0] - fdTheta);
                        double errT = Math.abs(analyticT[0] - fdT);

                        maxAbs = updateMax(maxAbs, errRho);
                        if (maxAbs == errRho) worst = describe(rotationLabels[r], rho, theta, t, "rho", errRho);
                        maxAbs = updateMax(maxAbs, errTheta);
                        if (maxAbs == errTheta) worst = describe(rotationLabels[r], rho, theta, t, "theta", errTheta);
                        maxAbs = updateMax(maxAbs, errT);
                        if (maxAbs == errT) worst = describe(rotationLabels[r], rho, theta, t, "t", errT);

                        if (errRho > WARN_ABS || errTheta > WARN_ABS || errT > WARN_ABS) {
                            warnings++;
                        }

                        final boolean filteredCase = isFilteredCase(rotations[r], rho, theta);
                        if (filteredCase) {
                            filteredCases++;
                            filteredMaxAbs = updateMax(filteredMaxAbs, errRho);
                            if (filteredMaxAbs == errRho) filteredWorst = describe(rotationLabels[r], rho, theta, t, "rho", errRho);
                            filteredMaxAbs = updateMax(filteredMaxAbs, errTheta);
                            if (filteredMaxAbs == errTheta) filteredWorst = describe(rotationLabels[r], rho, theta, t, "theta", errTheta);
                            filteredMaxAbs = updateMax(filteredMaxAbs, errT);
                            if (filteredMaxAbs == errT) filteredWorst = describe(rotationLabels[r], rho, theta, t, "t", errT);
                            if (errRho > WARN_ABS || errTheta > WARN_ABS || errT > WARN_ABS) {
                                filteredWarnings++;
                            }
                        }

                        for (int idx = 0; idx < analyticR.length; ++idx) {
                            final double fdR;
                            try {
                                fdR = numericalGradient(model.likelihood, drift.getRotationMatrixParameter(), idx);
                            } catch (RuntimeException ex) {
                                invalid++;
                                continue;
                            }
                            double errR = Math.abs(analyticR[idx] - fdR);
                            maxAbs = updateMax(maxAbs, errR);
                            if (maxAbs == errR) worst = describe(rotationLabels[r], rho, theta, t, "R[" + idx + "]", errR);
                            if (errR > WARN_ABS) {
                                warnings++;
                            }
                            if (filteredCase) {
                                filteredMaxAbs = updateMax(filteredMaxAbs, errR);
                                if (filteredMaxAbs == errR) filteredWorst = describe(rotationLabels[r], rho, theta, t, "R[" + idx + "]", errR);
                                if (errR > WARN_ABS) {
                                    filteredWarnings++;
                                }
                            }
                        }

                        cases++;
                    }
                }
            }
        }

        System.out.println("cases = " + cases);
        System.out.println("invalid = " + invalid);
        System.out.println("warnings(abs_err > " + WARN_ABS + ") = " + warnings);
        System.out.println("max_abs_err = " + maxAbs);
        System.out.println("worst_case = " + worst);
        System.out.println("filtered_cases = " + filteredCases);
        System.out.println("filtered_warnings(abs_err > " + WARN_ABS + ") = " + filteredWarnings);
        System.out.println("filtered_max_abs_err = " + filteredMaxAbs);
        System.out.println("filtered_worst_case = " + filteredWorst);
    }

    private static double updateMax(double current, double candidate) {
        return Math.max(current, candidate);
    }

    private static String describe(String rotationLabel,
                                   double rho,
                                   double theta,
                                   double t,
                                   String parameter,
                                   double err) {
        return "rotation=" + rotationLabel
                + ", rho=" + rho
                + ", theta=" + theta
                + ", t=" + t
                + ", parameter=" + parameter
                + ", abs_err=" + err;
    }

    private static boolean isFilteredCase(double[][] rotation,
                                          double rho,
                                          double theta) {
        final double det = rotation[0][0] * rotation[1][1] - rotation[0][1] * rotation[1][0];
        return rho >= 0.02 && Math.abs(theta) <= 1.45 && Math.abs(det) >= 0.02;
    }
}
