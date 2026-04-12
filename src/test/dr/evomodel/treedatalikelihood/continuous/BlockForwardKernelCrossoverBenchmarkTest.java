package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalExpSolver;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Random;

public class BlockForwardKernelCrossoverBenchmarkTest extends TestCase {

    private static final int WARMUP = 2_000;
    private static final int ITERS = 20_000;

    public BlockForwardKernelCrossoverBenchmarkTest(final String name) {
        super(name);
    }

    public void testCrossoverReportIsFinite() {
        final String report = buildReport();
        System.out.println(report);
        assertTrue(report.contains("d=2"));
        assertTrue(report.contains("d=12"));
        assertFalse(report.contains("NaN"));
        assertFalse(report.contains("Infinity"));
    }

    public static Test suite() {
        return new TestSuite(BlockForwardKernelCrossoverBenchmarkTest.class);
    }

    public static void main(final String[] args) {
        junit.textui.TestRunner.run(suite());
    }

    private static String buildReport() {
        final StringBuilder sb = new StringBuilder();
        for (int d = 2; d <= 12; ++d) {
            final Timing timing = runCase(d);
            sb.append(String.format("d=%d  generic=%.1f ns  small=%.1f ns  speedup=%.2fx%n",
                    d, timing.genericNs, timing.smallNs, timing.genericNs / timing.smallNs));
        }
        return sb.toString();
    }

    private static Timing runCase(final int d) {
        final Random random = new Random(17L + d);
        final double[] blockDParams = randomBlockParams(d, random);
        final double[] r = randomOrthogonal(d, random);
        final double[] sigmaEig = randomSPD(d, random);
        final double[] optimum = randomVector(d, random);
        final double dt = 0.17;

        final KernelState generic = new KernelState(d);
        final KernelState small = new KernelState(d);

        for (int i = 0; i < WARMUP; ++i) {
            runGeneric(blockDParams, r, sigmaEig, optimum, dt, generic);
            runSmall(blockDParams, r, sigmaEig, optimum, dt, small);
        }

        long start = System.nanoTime();
        for (int i = 0; i < ITERS; ++i) {
            runGeneric(blockDParams, r, sigmaEig, optimum, dt, generic);
        }
        final double genericNs = (System.nanoTime() - start) / (double) ITERS;

        start = System.nanoTime();
        for (int i = 0; i < ITERS; ++i) {
            runSmall(blockDParams, r, sigmaEig, optimum, dt, small);
        }
        final double smallNs = (System.nanoTime() - start) / (double) ITERS;

        return new Timing(genericNs, smallNs);
    }

    private static void runGeneric(final double[] blockDParams,
                                   final double[] r,
                                   final double[] sigmaEig,
                                   final double[] optimum,
                                   final double dt,
                                   final KernelState s) {
        final int d = s.d;
        BlockDiagonalExpSolver.expBlockMatrix(blockDParams, d, dt, s.qEig);
        System.arraycopy(r, 0, s.r.data, 0, d * d);

        CommonOps.mult(s.r, s.qEig, s.temp);
        CommonOps.multTransB(s.temp, s.r, s.actualization);

        s.vEig.set(sigmaEigMatrix(sigmaEig, d, s.sigmaEig));
        CommonOps.mult(s.qEig, s.sigmaEig, s.temp);
        CommonOps.multTransB(s.temp, s.qEig, s.branchVariance);
        CommonOps.subtract(s.vEig, s.branchVariance, s.vEig);
        CommonOps.mult(s.r, s.vEig, s.temp);
        CommonOps.multTransB(s.temp, s.r, s.branchVariance);
        CommonOps.multTransA(s.r, wrapVector(optimum, d, s.optimum), s.optimumEig);
        CommonOps.mult(s.qEig, s.optimumEig, s.displacement);
        CommonOps.mult(s.r, s.displacement, s.optimumEig);
        CommonOps.subtract(wrapVector(optimum, d, s.optimum), s.optimumEig, s.displacement);
    }

    private static void runSmall(final double[] blockDParams,
                                 final double[] r,
                                 final double[] sigmaEig,
                                 final double[] optimum,
                                 final double dt,
                                 final KernelState s) {
        final int d = s.d;
        final double[] qEig = s.qEig.data;
        final double[] temp = s.temp.data;
        final double[] actualization = s.actualization.data;
        final double[] variance = s.branchVariance.data;
        final double[] vEig = s.vEig.data;
        final double[] optimumEig = s.optimumEig.data;
        final double[] displacement = s.displacement.data;

        BlockDiagonalExpSolver.expBlockMatrix(blockDParams, d, dt, s.qEig);
        conjugateOrthogonalSmall(r, qEig, d, temp, actualization);

        multiplySmall(qEig, sigmaEig, d, temp);
        multiplyTransposedRightSmall(temp, qEig, d, variance);
        for (int i = 0; i < d * d; ++i) {
            vEig[i] = sigmaEig[i] - variance[i];
        }
        conjugateOrthogonalSmall(r, vEig, d, temp, variance);

        multiplyTransposeMatrixVectorSmall(r, optimum, 0, d, optimumEig);
        multiplyMatrixVectorSmall(qEig, optimumEig, d, displacement);
        multiplyMatrixVectorSmall(r, displacement, d, optimumEig);
        for (int i = 0; i < d; ++i) {
            displacement[i] = optimum[i] - optimumEig[i];
        }
    }

    private static DenseMatrix64F sigmaEigMatrix(final double[] sigmaEig, final int d, final DenseMatrix64F out) {
        System.arraycopy(sigmaEig, 0, out.data, 0, d * d);
        return out;
    }

    private static DenseMatrix64F wrapVector(final double[] vector, final int d, final DenseMatrix64F out) {
        System.arraycopy(vector, 0, out.data, 0, d);
        return out;
    }

    private static double[] randomBlockParams(final int d, final Random random) {
        final double[] out = new double[3 * d - 2];
        int i = 0;
        while (i < d) {
            if (i + 1 < d) {
                out[i] = 0.6 + 0.5 * random.nextDouble();
                out[i + 1] = 0.6 + 0.5 * random.nextDouble();
                out[d + i / 2] = 0.2 * (random.nextDouble() - 0.5);
                i += 2;
            } else {
                out[i] = 0.6 + 0.5 * random.nextDouble();
                i += 1;
            }
        }
        final int rotOffset = 2 * d - 1;
        for (int j = rotOffset; j < out.length; ++j) {
            out[j] = 0.5 * (random.nextDouble() - 0.5);
        }
        return out;
    }

    private static double[] randomOrthogonal(final int d, final Random random) {
        final double[] a = new double[d * d];
        for (int i = 0; i < a.length; ++i) {
            a[i] = random.nextGaussian();
        }
        final double[] q = new double[d * d];
        final double[] v = new double[d];
        for (int j = 0; j < d; ++j) {
            for (int i = 0; i < d; ++i) {
                v[i] = a[i * d + j];
            }
            for (int k = 0; k < j; ++k) {
                double dot = 0.0;
                for (int i = 0; i < d; ++i) {
                    dot += q[i * d + k] * v[i];
                }
                for (int i = 0; i < d; ++i) {
                    v[i] -= dot * q[i * d + k];
                }
            }
            double norm = 0.0;
            for (int i = 0; i < d; ++i) {
                norm += v[i] * v[i];
            }
            norm = Math.sqrt(norm);
            for (int i = 0; i < d; ++i) {
                q[i * d + j] = v[i] / norm;
            }
        }
        return q;
    }

    private static double[] randomSPD(final int d, final Random random) {
        final double[] a = new double[d * d];
        final double[] out = new double[d * d];
        for (int i = 0; i < a.length; ++i) {
            a[i] = 0.3 * random.nextGaussian();
        }
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += a[i * d + k] * a[j * d + k];
                }
                out[i * d + j] = sum;
            }
            out[i * d + i] += 1.0;
        }
        return out;
    }

    private static double[] randomVector(final int d, final Random random) {
        final double[] out = new double[d];
        for (int i = 0; i < d; ++i) {
            out[i] = random.nextGaussian();
        }
        return out;
    }

    private static void conjugateOrthogonalSmall(final double[] r,
                                                 final double[] matrix,
                                                 final int d,
                                                 final double[] temp,
                                                 final double[] out) {
        multiplySmall(r, matrix, d, temp);
        multiplyTransposedRightSmall(temp, r, d, out);
    }

    private static void multiplySmall(final double[] left,
                                      final double[] right,
                                      final int d,
                                      final double[] out) {
        for (int i = 0; i < d; ++i) {
            final int rowOffset = i * d;
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += left[rowOffset + k] * right[k * d + j];
                }
                out[rowOffset + j] = sum;
            }
        }
    }

    private static void multiplyTransposedRightSmall(final double[] left,
                                                     final double[] right,
                                                     final int d,
                                                     final double[] out) {
        for (int i = 0; i < d; ++i) {
            final int rowOffset = i * d;
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                final int rightRowOffset = j * d;
                for (int k = 0; k < d; ++k) {
                    sum += left[rowOffset + k] * right[rightRowOffset + k];
                }
                out[rowOffset + j] = sum;
            }
        }
    }

    private static void multiplyTransposeMatrixVectorSmall(final double[] matrix,
                                                           final double[] vector,
                                                           final int vectorOffset,
                                                           final int d,
                                                           final double[] out) {
        for (int j = 0; j < d; ++j) {
            double sum = 0.0;
            for (int i = 0; i < d; ++i) {
                sum += matrix[i * d + j] * vector[vectorOffset + i];
            }
            out[j] = sum;
        }
    }

    private static void multiplyMatrixVectorSmall(final double[] matrix,
                                                  final double[] vector,
                                                  final int d,
                                                  final double[] out) {
        for (int i = 0; i < d; ++i) {
            final int rowOffset = i * d;
            double sum = 0.0;
            for (int j = 0; j < d; ++j) {
                sum += matrix[rowOffset + j] * vector[j];
            }
            out[i] = sum;
        }
    }

    private static final class Timing {
        final double genericNs;
        final double smallNs;

        private Timing(final double genericNs, final double smallNs) {
            this.genericNs = genericNs;
            this.smallNs = smallNs;
        }
    }

    private static final class KernelState {
        final int d;
        final DenseMatrix64F r;
        final DenseMatrix64F qEig;
        final DenseMatrix64F sigmaEig;
        final DenseMatrix64F vEig;
        final DenseMatrix64F temp;
        final DenseMatrix64F actualization;
        final DenseMatrix64F branchVariance;
        final DenseMatrix64F optimum;
        final DenseMatrix64F optimumEig;
        final DenseMatrix64F displacement;

        private KernelState(final int d) {
            this.d = d;
            this.r = new DenseMatrix64F(d, d);
            this.qEig = new DenseMatrix64F(d, d);
            this.sigmaEig = new DenseMatrix64F(d, d);
            this.vEig = new DenseMatrix64F(d, d);
            this.temp = new DenseMatrix64F(d, d);
            this.actualization = new DenseMatrix64F(d, d);
            this.branchVariance = new DenseMatrix64F(d, d);
            this.optimum = new DenseMatrix64F(d, 1);
            this.optimumEig = new DenseMatrix64F(d, 1);
            this.displacement = new DenseMatrix64F(d, 1);
        }
    }
}
