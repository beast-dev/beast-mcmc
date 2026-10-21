package test.dr.inference.timeseries;

import static test.dr.inference.timeseries.OUTimeSeriesTestSupport.latent;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.inference.model.CompoundSymmetricMatrix;
import dr.inference.model.DiagonalFloorMatrixView;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.model.ScaledIdentityMatrixParameter;
import dr.inference.model.Variable;
import dr.inference.model.VariableListener;
import dr.inference.timeseries.core.BasicTimeSeriesModel;
import dr.inference.timeseries.core.IrregularTimeGrid;
import dr.inference.timeseries.engine.kalman.GaussianForwardComputationMode;
import dr.inference.timeseries.likelihood.GaussianGradientComputationMode;
import dr.inference.timeseries.likelihood.GaussianSmootherComputationMode;
import dr.inference.timeseries.likelihood.GaussianTimeSeriesLikelihoodFactory;
import dr.inference.timeseries.likelihood.TimeSeriesLikelihood;
import dr.inference.timeseries.model.gaussian.LinearGaussianObservationModel;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * The declared diffusion floor of {@code ouProcessModel}: the canonical likelihood and every
 * analytical gradient under {@code diffusionFloor = eps} must equal those of the same model with
 * {@code Q + eps I} supplied explicitly, a zero floor must leave the model untouched, the
 * compound-symmetric (Cholesky) diffusion chart must keep consistent gradients under the floor,
 * and the moment-form engines must refuse a floor they would not apply.
 */
public class DiffusionFloorTest extends TestCase {

    private static final double NAN = Double.NaN;
    private static final double EPS = 2.5e-3;

    public DiffusionFloorTest(final String name) {
        super(name);
    }

    public static Test suite() {
        return new TestSuite(DiffusionFloorTest.class);
    }

    public void testFlooredLikelihoodEqualsExplicitlyFlooredDiffusion() {
        final Fixture floored = denseFixture("floor.a", EPS, false);
        final Fixture explicit = denseFixture("floor.b", 0.0, true);
        final Fixture raw = denseFixture("floor.c", 0.0, false);
        assertEquals("floor eps vs explicit Q + eps I",
                explicit.likelihood.getLogLikelihood(), floored.likelihood.getLogLikelihood(), 1e-9);
        assertTrue("the floor must change the likelihood in this fixture",
                Math.abs(raw.likelihood.getLogLikelihood() - floored.likelihood.getLogLikelihood()) > 1e-4);
        final double[] viewFlat = floored.process.getDiffusionMatrix().getParameterValues();
        final double[] explicitFlat = explicit.process.getDiffusionMatrix().getParameterValues();
        assertEquals(explicitFlat.length, viewFlat.length);
        for (int i = 0; i < viewFlat.length; ++i) {
            assertEquals("view entry " + i, explicitFlat[i], viewFlat[i], 1e-15);
        }
    }

    public void testZeroFloorIsIdentity() {
        final Fixture zero = denseFixture("floor.zero", 0.0, false);
        assertSame("zero floor must hand engines the raw matrix",
                zero.diffusion, zero.process.getDiffusionMatrix());
        assertSame(zero.diffusion, zero.process.getRawDiffusionMatrix());
        assertEquals(0.0, zero.process.getDiffusionFloor(), 0.0);
        final OUProcessModel plain = new OUProcessModel("floor.plain", 3, zero.drift, zero.diffusion, zero.mean,
                matrix("floor.plain.P0", scaledIdentity(3, 0.05)));
        final TimeSeriesLikelihood plainLikelihood = canonicalLikelihood("floor.plain", plain, zero.observation, zero.grid);
        assertEquals(plainLikelihood.getLogLikelihood(), zero.likelihood.getLogLikelihood(), 0.0);
    }

    public void testFlooredGradientsMatchExplicitAndFiniteDifferences() {
        final Fixture floored = denseFixture("floor.grad", EPS, false);
        final Fixture explicit = denseFixture("floor.grad.explicit", 0.0, true);
        final double[] driftFloored = floored.likelihood.getGradientWrt(floored.drift).getGradientLogDensity(null);
        final double[] driftExplicit = explicit.likelihood.getGradientWrt(explicit.drift).getGradientLogDensity(null);
        final double[] diffusionFloored = floored.likelihood.getGradientWrt(floored.diffusion).getGradientLogDensity(null);
        final double[] diffusionExplicit = explicit.likelihood.getGradientWrt(explicit.diffusion).getGradientLogDensity(null);
        final double[] meanFloored = floored.likelihood.getGradientWrt(floored.mean).getGradientLogDensity(null);
        for (int i = 0; i < driftFloored.length; ++i) {
            assertEquals("drift gradient idx=" + i, driftExplicit[i], driftFloored[i], 1e-8);
            assertEquals("drift gradient vs finite difference idx=" + i,
                    fdFourthOrder(floored.likelihood, floored.drift, i, 1e-3), driftFloored[i], 1e-6);
        }
        for (int i = 0; i < diffusionFloored.length; ++i) {
            assertEquals("diffusion gradient idx=" + i, diffusionExplicit[i], diffusionFloored[i], 1e-8);
            assertEquals("diffusion gradient vs finite difference idx=" + i,
                    fdFourthOrder(floored.likelihood, floored.diffusion, i, 1e-4), diffusionFloored[i], 1e-6);
        }
        for (int i = 0; i < meanFloored.length; ++i) {
            assertEquals("mean gradient vs finite difference idx=" + i,
                    fdFourthOrder(floored.likelihood, floored.mean, i, 1e-3), meanFloored[i], 1e-7);
        }
    }

    public void testCompoundSymmetricDiffusionGradientsUnderFloorMatchFiniteDifferences() {
        final Fixture f = compoundSymmetricFixture("floor.cs", EPS);
        final double[] diagonalGradient = f.likelihood.getGradientWrt(f.csDiagonal).getGradientLogDensity(null);
        final double[] offDiagonalGradient = f.likelihood.getGradientWrt(f.csOffDiagonal).getGradientLogDensity(null);
        for (int i = 0; i < diagonalGradient.length; ++i) {
            assertEquals("compound-symmetric diagonal gradient idx=" + i,
                    fdFourthOrder(f.likelihood, f.csDiagonal, i, 1e-4), diagonalGradient[i], 1e-6);
        }
        for (int i = 0; i < offDiagonalGradient.length; ++i) {
            assertEquals("compound-symmetric off-diagonal gradient idx=" + i,
                    fdFourthOrder(f.likelihood, f.csOffDiagonal, i, 1e-4), offDiagonalGradient[i], 1e-6);
        }
        // The engines see Q + eps I: the view's diagonal exceeds the raw chart's by exactly eps.
        final MatrixParameterInterface view = f.process.getDiffusionMatrix();
        final MatrixParameterInterface rawMatrix = f.process.getRawDiffusionMatrix();
        for (int i = 0; i < 3; ++i) {
            assertEquals(rawMatrix.getParameterValue(i, i) + EPS, view.getParameterValue(i, i), 1e-15);
            assertEquals(rawMatrix.getParameterValue(i, (i + 1) % 3), view.getParameterValue(i, (i + 1) % 3), 0.0);
        }
    }

    public void testViewValuesAndChangeEvents() {
        final MatrixParameter raw = matrix("view.raw", new double[][]{
                {0.40, 0.03, 0.00},
                {0.03, 0.30, 0.02},
                {0.00, 0.02, 0.50}
        });
        final DiagonalFloorMatrixView view = new DiagonalFloorMatrixView("view", raw, EPS);
        assertEquals(raw.getDimension(), view.getDimension());
        final double[][] expected = raw.getParameterAsMatrix();
        for (int i = 0; i < 3; ++i) {
            expected[i][i] += EPS;
        }
        final double[][] actual = view.getParameterAsMatrix();
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                assertEquals("view[" + i + "][" + j + "]", expected[i][j], actual[i][j], 0.0);
            }
        }
        final double[] flat = view.getParameterValues();
        final double[] rawFlat = raw.getParameterValues();
        double floorMass = 0.0;
        for (int i = 0; i < flat.length; ++i) {
            floorMass += flat[i] - rawFlat[i];
        }
        assertEquals("exactly three diagonal entries carry the floor", 3.0 * EPS, floorMass, 1e-15);
        final int[] events = new int[]{0};
        view.addVariableListener(new VariableListener() {
            @Override
            public void variableChangedEvent(final Variable variable, final int index, final Parameter.ChangeType type) {
                ++events[0];
            }
        });
        raw.setParameterValue(1, 1, 0.31);
        assertTrue("view must forward the wrapped matrix's change events", events[0] > 0);
        assertEquals(0.31 + EPS, view.getParameterValue(1, 1), 0.0);
        try {
            new DiagonalFloorMatrixView("bad", raw, -1.0);
            fail("negative floor accepted");
        } catch (IllegalArgumentException rejected) {
            // expected
        }
    }

    public void testMomentModesRejectFloor() {
        final Fixture f = denseFixture("floor.moment", EPS, false);
        final BasicTimeSeriesModel model = new BasicTimeSeriesModel("floor.moment.model", latent(f.process), f.observation, f.grid);
        try {
            GaussianTimeSeriesLikelihoodFactory.create("floor.moment.like", model,
                    GaussianForwardComputationMode.MOMENT,
                    GaussianSmootherComputationMode.MOMENT,
                    GaussianGradientComputationMode.MOMENT_ANALYTICAL);
            fail("moment-form engines accepted a diffusion floor");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("diffusionFloor"));
        }
        try {
            new OUProcessModel("floor.neg", 3, f.drift, f.diffusion, f.mean,
                    matrix("floor.neg.P0", scaledIdentity(3, 0.05)), -EPS);
            fail("negative floor accepted by the model");
        } catch (IllegalArgumentException expected) {
            // expected
        }
    }

    // ------------------------------------------------------------------ fixtures

    private static final class Fixture {
        MatrixParameter drift;
        MatrixParameterInterface diffusion;
        Parameter csDiagonal;
        Parameter csOffDiagonal;
        Parameter mean;
        OUProcessModel process;
        LinearGaussianObservationModel observation;
        IrregularTimeGrid grid;
        TimeSeriesLikelihood likelihood;
    }

    private static final double[][] Q = new double[][]{
            {0.40, 0.03, 0.00},
            {0.03, 0.30, 0.02},
            {0.00, 0.02, 0.50}
    };

    private static Fixture denseFixture(final String p, final double floor, final boolean explicitFloor) {
        final Fixture f = new Fixture();
        f.drift = matrix(p + ".A", new double[][]{
                {0.90, -0.20, 0.10},
                {0.15, 1.30, -0.10},
                {0.05, 0.20, 0.70}
        });
        final double[][] q = new double[3][3];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                q[i][j] = Q[i][j] + (explicitFloor && i == j ? EPS : 0.0);
            }
        }
        f.diffusion = matrix(p + ".Q", q);
        f.mean = new Parameter.Default(p + ".mu", new double[]{0.01, -0.02, 0.00});
        f.process = new OUProcessModel(p + ".ou", 3, f.drift, f.diffusion, f.mean,
                matrix(p + ".P0", scaledIdentity(3, 0.05)),
                OUProcessModel.CovarianceGradientMethod.LYAPUNOV_ADJOINT, floor);
        finishFixture(p, f);
        return f;
    }

    private static Fixture compoundSymmetricFixture(final String p, final double floor) {
        final Fixture f = new Fixture();
        f.drift = matrix(p + ".A", new double[][]{
                {0.90, -0.20, 0.10},
                {0.15, 1.30, -0.10},
                {0.05, 0.20, 0.70}
        });
        f.csDiagonal = new Parameter.Default(p + ".Q.diag", new double[]{0.40, 0.30, 0.50});
        f.csDiagonal.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 3));
        f.csOffDiagonal = new Parameter.Default(p + ".Q.off", new double[]{0.30, -0.20, 0.10});
        f.csOffDiagonal.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 3));
        f.diffusion = new CompoundSymmetricMatrix(f.csDiagonal, f.csOffDiagonal, true, true);
        f.mean = new Parameter.Default(p + ".mu", new double[]{0.01, -0.02, 0.00});
        f.process = new OUProcessModel(p + ".ou", 3, f.drift, f.diffusion, f.mean,
                matrix(p + ".P0", scaledIdentity(3, 0.05)),
                OUProcessModel.CovarianceGradientMethod.LYAPUNOV_ADJOINT, floor);
        finishFixture(p, f);
        return f;
    }

    private static void finishFixture(final String p, final Fixture f) {
        final MatrixParameter design = matrix(p + ".H", scaledIdentity(3, 1.0));
        final Parameter noiseVariance = new Parameter.Default(p + ".obsVar", new double[]{1.0e-5});
        final ScaledIdentityMatrixParameter noise = new ScaledIdentityMatrixParameter(p + ".R", 3, noiseVariance);
        final MatrixParameter observations = matrix(p + ".Y", new double[][]{
                {NAN, 0.013, NAN, 0.011, NAN, NAN, 0.009, NAN},
                {NAN, NAN, NAN, -0.021, -0.024, NAN, NAN, -0.019},
                {NAN, NAN, 0.004, NAN, NAN, 0.006, NAN, 0.002}
        });
        f.grid = new IrregularTimeGrid(new double[]{
                0.0, 1.6389e-4, 2.0e-4, 2.0278e-4, 9.1e-4, 1.3e-3, 1.3e-3 + 2.8e-7, 4.5e-3});
        f.observation = new LinearGaussianObservationModel(p + ".obs", 3, design, noise, observations);
        f.likelihood = canonicalLikelihood(p, f.process, f.observation, f.grid);
    }

    private static TimeSeriesLikelihood canonicalLikelihood(final String name,
                                                            final OUProcessModel process,
                                                            final LinearGaussianObservationModel obs,
                                                            final IrregularTimeGrid grid) {
        final BasicTimeSeriesModel model = new BasicTimeSeriesModel(name + ".model", latent(process), obs, grid);
        return GaussianTimeSeriesLikelihoodFactory.create(
                name + ".likelihood",
                model,
                GaussianForwardComputationMode.CANONICAL,
                GaussianSmootherComputationMode.CANONICAL,
                GaussianGradientComputationMode.CANONICAL_ANALYTICAL);
    }

    // ------------------------------------------------------------------ helpers

    private static double fdFourthOrder(final Likelihood like, final Parameter param, final int index, final double h) {
        final double x = param.getParameterValue(index);
        final double p2 = evalAt(like, param, index, x + 2.0 * h);
        final double p1 = evalAt(like, param, index, x + h);
        final double m1 = evalAt(like, param, index, x - h);
        final double m2 = evalAt(like, param, index, x - 2.0 * h);
        evalAt(like, param, index, x);
        return (-p2 + 8.0 * p1 - 8.0 * m1 + m2) / (12.0 * h);
    }

    private static double evalAt(final Likelihood like, final Parameter param, final int index, final double value) {
        param.setParameterValue(index, value);
        return like.getLogLikelihood();
    }

    private static MatrixParameter matrix(final String name, final double[][] values) {
        final MatrixParameter mp = new MatrixParameter(name, values.length, values[0].length);
        for (int i = 0; i < values.length; ++i) {
            for (int j = 0; j < values[i].length; ++j) {
                mp.setParameterValue(i, j, values[i][j]);
            }
        }
        return mp;
    }

    private static double[][] scaledIdentity(final int dim, final double scale) {
        final double[][] out = new double[dim][dim];
        for (int i = 0; i < dim; ++i) {
            out[i][i] = scale;
        }
        return out;
    }
}
