package dr.evomodel.treedatalikelihood.continuous.backprop;

/**
 * Small-block numerical integrator for Fréchet derivatives in the D-basis.
 *
 * <p>The integrals are only over 1x1 and 2x2 blocks, so fixed-order Gauss-Legendre
 * quadrature is sufficient while keeping the overall cost quadratic in the number
 * of blocks.</p>
 */
final class BlockDiagonalFrechetIntegrator {

    private static final double[] NODES = {
            0.0,
            -0.5384693101056831, 0.5384693101056831,
            -0.9061798459386640, 0.9061798459386640
    };
    private static final double[] WEIGHTS = {
            0.5688888888888889,
            0.4786286704993665, 0.4786286704993665,
            0.2369268850561891, 0.2369268850561891
    };

    private BlockDiagonalFrechetIntegrator() {
    }

    static void integrate1x1_2x2(final double lambdaI,
                                 final double a, final double b, final double c, final double d,
                                 final double e1, final double e2,
                                 final double[] out2,
                                 final double[] scratch4) {
        out2[0] = 0.0;
        out2[1] = 0.0;
        for (int k = 0; k < NODES.length; ++k) {
            final double s = 0.5 * (NODES[k] + 1.0);
            final double w = 0.5 * WEIGHTS[k];
            final double left = Math.exp((1.0 - s) * lambdaI);
            fillExpOfScaled2x2(a, b, c, d, s, scratch4);
            out2[0] += w * left * (e1 * scratch4[0] + e2 * scratch4[2]);
            out2[1] += w * left * (e1 * scratch4[1] + e2 * scratch4[3]);
        }
    }

    static void integrate2x2_1x1(final double a, final double b, final double c, final double d,
                                 final double lambdaJ,
                                 final double e1, final double e2,
                                 final double[] out2,
                                 final double[] scratch4) {
        out2[0] = 0.0;
        out2[1] = 0.0;
        for (int k = 0; k < NODES.length; ++k) {
            final double s = 0.5 * (NODES[k] + 1.0);
            final double w = 0.5 * WEIGHTS[k];
            fillExpOfScaled2x2(a, b, c, d, 1.0 - s, scratch4);
            final double right = Math.exp(s * lambdaJ);
            out2[0] += w * right * (scratch4[0] * e1 + scratch4[1] * e2);
            out2[1] += w * right * (scratch4[2] * e1 + scratch4[3] * e2);
        }
    }

    static void integrate2x2_2x2(final double ai, final double bi, final double ci, final double di,
                                 final double aj, final double bj, final double cj, final double dj,
                                 final double e11, final double e12, final double e21, final double e22,
                                 final double[] out4,
                                 final double[] left4,
                                 final double[] right4) {
        out4[0] = 0.0;
        out4[1] = 0.0;
        out4[2] = 0.0;
        out4[3] = 0.0;

        for (int k = 0; k < NODES.length; ++k) {
            final double s = 0.5 * (NODES[k] + 1.0);
            final double w = 0.5 * WEIGHTS[k];
            fillExpOfScaled2x2(ai, bi, ci, di, 1.0 - s, left4);
            fillExpOfScaled2x2(aj, bj, cj, dj, s, right4);

            final double m11 = left4[0] * e11 + left4[1] * e21;
            final double m12 = left4[0] * e12 + left4[1] * e22;
            final double m21 = left4[2] * e11 + left4[3] * e21;
            final double m22 = left4[2] * e12 + left4[3] * e22;

            out4[0] += w * (m11 * right4[0] + m12 * right4[2]);
            out4[1] += w * (m11 * right4[1] + m12 * right4[3]);
            out4[2] += w * (m21 * right4[0] + m22 * right4[2]);
            out4[3] += w * (m21 * right4[1] + m22 * right4[3]);
        }
    }

    private static void fillExpOfScaled2x2(final double a,
                                           final double b,
                                           final double c,
                                           final double d,
                                           final double scale,
                                           final double[] out4) {
        BlockDiagonalExpSolver.fillExp2x2GeneralBlock(-a, -b, -c, -d, scale, out4, 2, 0);
    }
}
