package dr.evomodel.treedatalikelihood.continuous.backprop;

import org.ejml.data.DenseMatrix64F;

import java.util.Arrays;

/**
 * Computes exp(-D t) for a block-diagonal matrix D with 1x1 and 2x2 real blocks.
 *
 * <p>This class is intentionally independent of any particular parameterization
 * that produced D.  It only consumes:</p>
 *
 * <ul>
 *   <li>a fixed block topology ({@link BlockStructure})</li>
 *   <li>the raw compressed entries of D in the layout<br>
 *       [diag_0..diag_{d-1}, upper_0..upper_{d-2}, lower_0..lower_{d-2}]</li>
 * </ul>
 *
 * <p>The efficient usage pattern is to construct one solver with a fixed block
 * structure and reuse it across calls as the raw entries change.</p>
 */
public final class BlockDiagonalExpSolver {

    public static final double DEFAULT_EPS_BLOCK = 1.0e-12;

    private final BlockStructure structure;
    private final int dim;

//    public BlockDiagonalExpSolver(final int dim,
//                                  final int[] blockStarts,
//                                  final int[] blockSizes) {
//        this(new BlockStructure(dim, blockStarts, blockSizes));
//    }

    public BlockDiagonalExpSolver(final BlockStructure structure) {
        this.structure = structure;
        this.dim = structure.getDim();
    }

    public BlockStructure getStructure() {
        return structure;
    }

    /**
     * Convenience wrapper that infers the 1x1 / 2x2 block topology from the raw
     * compressed entries and performs a one-off exponential computation.
     */
    public static void expBlockMatrix(final double[] blockDParams,
                                      final int dim,
                                      final double t,
                                      final DenseMatrix64F out) {
        final BlockStructure structure = BlockStructure.inferFromCompressed(blockDParams, dim, DEFAULT_EPS_BLOCK);
        new BlockDiagonalExpSolver(structure).compute(blockDParams, t, out);
    }

    /**
     * Compute out = exp(-D t).
     */
    public void compute(final double[] blockDParams,
                        final double t,
                        final DenseMatrix64F out) {
        validateCompressedParams(blockDParams, dim);
        validateSquareMatrix(out, dim, "out");

        final double[] outData = out.data;
        Arrays.fill(outData, 0.0);

        for (int b = 0; b < structure.getNumBlocks(); b++) {
            final int start = structure.getBlockStart(b);
            final int size = structure.getBlockSize(b);
            if (size == 1) {
                fillExp1x1Block(blockDParams[start], t, outData, dim, start);
            } else if (size == 2) {
                fillExp2x2Block(blockDParams, t, outData, start);
            } else {
                throw new IllegalStateException("Only block sizes 1 and 2 are supported, got " + size);
            }
        }
    }

    /**
     * Compute out = exp(-D t) in compressed block-diagonal layout:
     * [diag(0..dim-1), activeUpper(0..num2x2-1), activeLower(0..num2x2-1)].
     */
    public void computeCompressed(final double[] blockDParams,
                                  final double t,
                                  final double[] out) {
        validateCompressedParams(blockDParams, dim);
        final int num2x2Blocks = countTwoByTwoBlocks();
        final int expected = dim + 2 * num2x2Blocks;
        if (out.length != expected) {
            throw new IllegalArgumentException(
                    "Expected compressed exponential length " + expected + " but got " + out.length);
        }
        Arrays.fill(out, 0.0);

        int block2x2Index = 0;
        for (int b = 0; b < structure.getNumBlocks(); b++) {
            final int start = structure.getBlockStart(b);
            final int size = structure.getBlockSize(b);
            if (size == 1) {
                out[start] = Math.exp(-blockDParams[start] * t);
            } else if (size == 2) {
                fillExp2x2Compressed(blockDParams, t, out, start, block2x2Index, num2x2Blocks);
                block2x2Index++;
            } else {
                throw new IllegalStateException("Only block sizes 1 and 2 are supported, got " + size);
            }
        }
    }

    private static void fillExp1x1Block(final double a,
                                        final double t,
                                        final double[] outData,
                                        final int dim,
                                        final int i) {
        outData[i * dim + i] = Math.exp(-a * t);
    }

    private void fillExp2x2Block(final double[] p,
                                 final double t,
                                 final double[] outData,
                                 final int start) {
        final int upperOffset = dim;
        final int lowerOffset = dim + (dim - 1);

        final double a = p[start];
        final double d = p[start + 1];
        final double b = p[upperOffset + start];
        final double c = p[lowerOffset + start];

        fillExp2x2GeneralBlock(a, b, c, d, t, outData, dim, start);
    }

    private void fillExp2x2Compressed(final double[] p,
                                      final double t,
                                      final double[] outData,
                                      final int start,
                                      final int block2x2Index,
                                      final int num2x2Blocks) {
        final int upperOffset = dim;
        final int lowerOffset = dim + (dim - 1);

        final double a = p[start];
        final double d = p[start + 1];
        final double b = p[upperOffset + start];
        final double c = p[lowerOffset + start];

        fillExp2x2CompressedBlock(a, b, c, d, t, outData, start, block2x2Index, num2x2Blocks);
    }

    /**
     * Fill exp(-t B) for B = [[a,b],[c,d]].
     */
    static void fillExp2x2GeneralBlock(final double a,
                                       final double b,
                                       final double c,
                                       final double d,
                                       final double t,
                                       final double[] outData,
                                       final int dim,
                                       final int start) {
        final double tau = 0.5 * (a + d);
        final double diff = 0.5 * (a - d);
        final double delta2 = diff * diff + b * c;

        final double alpha;
        final double beta;

        if (Math.abs(delta2) < DEFAULT_EPS_BLOCK) {
            alpha = 1.0;
            beta = -t;
        } else if (delta2 > 0.0) {
            final double delta = Math.sqrt(delta2);
            final double z = delta * t;
            final double expZ    = Math.exp(z);
            final double invExpZ = 1.0 / expZ;
            alpha = 0.5 * (expZ + invExpZ);
            beta  = -0.5 * (expZ - invExpZ) / delta;
        } else {
            final double delta = Math.sqrt(-delta2);
            final double z = delta * t;
            alpha = Math.cos(z);
            beta = -Math.sin(z) / delta;
        }

        final double factor = Math.exp(-tau * t);
        final double k11 = a - tau;
        final double k22 = d - tau;

        final double m11 = factor * (alpha + beta * k11);
        final double m12 = factor * (beta * b);
        final double m21 = factor * (beta * c);
        final double m22 = factor * (alpha + beta * k22);

        final int i = start;
        final int j = start + 1;
        outData[i * dim + i] = m11;
        outData[i * dim + j] = m12;
        outData[j * dim + i] = m21;
        outData[j * dim + j] = m22;
    }

    private static void fillExp2x2CompressedBlock(final double a,
                                                  final double b,
                                                  final double c,
                                                  final double d,
                                                  final double t,
                                                  final double[] outData,
                                                  final int start,
                                                  final int block2x2Index,
                                                  final int num2x2Blocks) {
        final double tau = 0.5 * (a + d);
        final double diff = 0.5 * (a - d);
        final double delta2 = diff * diff + b * c;

        final double alpha;
        final double beta;

        if (Math.abs(delta2) < DEFAULT_EPS_BLOCK) {
            alpha = 1.0;
            beta = -t;
        } else if (delta2 > 0.0) {
            final double delta = Math.sqrt(delta2);
            final double z = delta * t;
            final double expZ    = Math.exp(z);
            final double invExpZ = 1.0 / expZ;
            alpha = 0.5 * (expZ + invExpZ);
            beta  = -0.5 * (expZ - invExpZ) / delta;
        } else {
            final double delta = Math.sqrt(-delta2);
            final double z = delta * t;
            alpha = Math.cos(z);
            beta = -Math.sin(z) / delta;
        }

        final double factor = Math.exp(-tau * t);
        final double k11 = a - tau;
        final double k22 = d - tau;

        outData[start] = factor * (alpha + beta * k11);
        outData[start + 1] = factor * (alpha + beta * k22);
        outData[outData.length - 2 * num2x2Blocks + block2x2Index] = factor * (beta * b);
        outData[outData.length - num2x2Blocks + block2x2Index] = factor * (beta * c);
    }

    private int countTwoByTwoBlocks() {
        int count = 0;
        for (int b = 0; b < structure.getNumBlocks(); ++b) {
            if (structure.getBlockSize(b) == 2) {
                count++;
            }
        }
        return count;
    }

    static void validateCompressedParams(final double[] blockDParams, final int dim) {
        final int expected = 3 * dim - 2;
        if (blockDParams.length != expected) {
            throw new IllegalArgumentException(
                    "Expected compressed block matrix length " + expected + " but got " + blockDParams.length);
        }
    }

    private static void validateSquareMatrix(final DenseMatrix64F matrix,
                                             final int dim,
                                             final String name) {
        if (matrix.numRows != dim || matrix.numCols != dim) {
            throw new IllegalArgumentException(
                    name + " must be " + dim + "x" + dim + " but is " + matrix.numRows + "x" + matrix.numCols);
        }
    }

    /**
     * Fixed block topology independent of the raw numerical entries.
     */
    public static final class BlockStructure {
        private final int dim;
        private final int numBlocks;
        private final int[] blockStarts;
        private final int[] blockSizes;

        public BlockStructure(final int dim,
                              final int[] blockStarts,
                              final int[] blockSizes) {
            if (dim <= 0) {
                throw new IllegalArgumentException("dim must be positive");
            }
            if (blockStarts == null || blockSizes == null) {
                throw new IllegalArgumentException("blockStarts and blockSizes must be non-null");
            }
            if (blockStarts.length != blockSizes.length) {
                throw new IllegalArgumentException("blockStarts and blockSizes must have the same length");
            }
            this.dim = dim;
            this.numBlocks = blockStarts.length;
            this.blockStarts = blockStarts.clone();
            this.blockSizes = blockSizes.clone();
            validate();
        }

        public static BlockStructure inferFromCompressed(final double[] blockDParams,
                                                         final int dim,
                                                         final double eps) {
            validateCompressedParams(blockDParams, dim);

            final int upperOffset = dim;
            final int lowerOffset = dim + (dim - 1);
            final int[] starts = new int[dim];
            final int[] sizes = new int[dim];

            int count = 0;
            int i = 0;
            while (i < dim) {
                starts[count] = i;
                if (i == dim - 1) {
                    sizes[count] = 1;
                    i += 1;
                } else {
                    final double u = blockDParams[upperOffset + i];
                    final double l = blockDParams[lowerOffset + i];
                    if (Math.abs(u) < eps && Math.abs(l) < eps) {
                        sizes[count] = 1;
                        i += 1;
                    } else {
                        sizes[count] = 2;
                        i += 2;
                    }
                }
                count += 1;
            }

            return new BlockStructure(dim,
                    Arrays.copyOf(starts, count),
                    Arrays.copyOf(sizes, count));
        }

        public int getDim() {
            return dim;
        }

        public int getNumBlocks() {
            return numBlocks;
        }

        public int getBlockStart(final int index) {
            return blockStarts[index];
        }

        public int getBlockSize(final int index) {
            return blockSizes[index];
        }

        public int[] getBlockStarts() {
            return blockStarts.clone();
        }

        public int[] getBlockSizes() {
            return blockSizes.clone();
        }

        private void validate() {
            int covered = 0;
            for (int b = 0; b < numBlocks; b++) {
                final int start = blockStarts[b];
                final int size = blockSizes[b];
                if (size != 1 && size != 2) {
                    throw new IllegalArgumentException("Only block sizes 1 and 2 are supported, got " + size);
                }
                if (start != covered) {
                    throw new IllegalArgumentException("Blocks must be contiguous and ordered; expected start "
                            + covered + " but got " + start);
                }
                covered += size;
            }
            if (covered != dim) {
                throw new IllegalArgumentException("Block structure covers " + covered + " dimensions, expected " + dim);
            }
        }
    }
}

//package dr.evomodel.treedatalikelihood.continuous.backprop;
//
//import org.ejml.data.DenseMatrix64F;
//
//public class BlockDiagonalExpSolver {
//    static double EPS_BLOCK = 1e-12;
//    /**
//     * Compute Q_eig(t) = exp(-D t) for a block-diagonal D whose blocks
//     * are 1×1 and general 2×2.
////     *
//     * Works directly with blockDParams array, avoiding D matrix construction.
//     *
//     * Contiguous nonzero off-diagonal pairs (i, i+1) define 2×2 blocks:
//     *   B_i = [ d_i      u_i     ]
//     *         [ l_i   d_{i+1}    ].
//     * Otherwise we treat d_i as a 1×1 block.
//     */
//    public static void expBlockMatrix(double[] blockDParams, int dim, double t,
//                                      DenseMatrix64F Q) {
//
//        final int upperOffset = dim;
//        final int lowerOffset = dim + (dim - 1);
//
//        double[] qData = Q.data;
//        java.util.Arrays.fill(qData, 0.0);
//
//        int i = 0;
//        while (i < dim) {
//            double u = (i < dim - 1) ? blockDParams[upperOffset + i] : 0.0;
//            double l = (i < dim - 1) ? blockDParams[lowerOffset + i] : 0.0;
//
//            if (i == dim - 1 || (Math.abs(u) < EPS_BLOCK && Math.abs(l) < EPS_BLOCK)) {
//                // 1×1 block
//                double a = blockDParams[i];
//                int idx = i * dim + i;
//                qData[idx] = Math.exp(-a * t);
//                i += 1;
//            } else {
//                // General 2×2 block:
//                //   [ a  b ]
//                //   [ c  d ]
//                double a = blockDParams[i];
//                double d = blockDParams[i + 1];
//                double b = u;
//                double c = l;
//
//                fillExp2x2GeneralBlock(a, b, c, d, t, qData, dim, i);
//                i += 2;
//            }
//        }
//    }
//
//    /**
//     * Fill the 2×2 block of Q corresponding to exp(-t B), where
//     *
//     *   B = [ a  b ]
//     *       [ c  d ].
//     *
//     * We use the standard 2×2 formula based on
//     *   τ  = (a + d)/2,
//     *   Δ² = ((a - d)/2)² + b c,
//     *   K  = B - τ I (tr(K) = 0, K² = Δ² I),
//     *
//     * and
//     *   exp(-t B) = e^{-τ t} [ α I + β K ],
//     *
//     * with
//     *   (Δ² > 0): α = cosh(δ t), β = -sinh(δ t)/δ
//     *   (Δ² < 0): α = cos(δ t),  β = -sin(δ t)/δ
//     *   (Δ² = 0): α = 1,         β = -t.
//     *
//     * Writes directly into the flat qData array in row-major order.
//     */
//    private static void fillExp2x2GeneralBlock(double a,
//                                               double b,
//                                               double c,
//                                               double d,
//                                               double t,
//                                               double[] qData,
//                                               int dim,
//                                               int offsetRowCol) {
//
//        double tau    = 0.5 * (a + d);
//        double diff   = 0.5 * (a - d);
//        double delta2 = diff * diff + b * c;
//
//        double alpha;
//        double beta;
//
//        if (Math.abs(delta2) < EPS_BLOCK) {
//            // Δ² ≈ 0  ⇒  K² ≈ 0 and exp(-t B) ≈ e^{-τ t}(I - t K)
//            alpha = 1.0;
//            beta  = -t;
//        } else if (delta2 > 0.0) {
//            double delta = Math.sqrt(delta2);
//            double dtv   = delta * t;
//            alpha = Math.cosh(dtv);
//            beta  = -Math.sinh(dtv) / delta;
//        } else {
//            double delta = Math.sqrt(-delta2);
//            double dtv   = delta * t;
//            alpha = Math.cos(dtv);
//            beta  = -Math.sin(dtv) / delta;
//        }
//
//        double factor = Math.exp(-tau * t);
//
//        // K = B - τ I
//        double k11 = a - tau;
//        double k22 = d - tau;
//
//        double m11 = alpha + beta * k11;
//        double m12 = beta * b;
//        double m21 = beta * c;
//        double m22 = alpha + beta * k22;
//
//        int i = offsetRowCol;
//        int j = offsetRowCol + 1;
//
//        int idx00 = i       * dim + i;
//        int idx01 = i       * dim + j;
//        int idx10 = (i + 1) * dim + i;
//        int idx11 = (i + 1) * dim + j;
//
//        qData[idx00] = factor * m11;
//        qData[idx01] = factor * m12;
//        qData[idx10] = factor * m21;
//        qData[idx11] = factor * m22;
//    }
//}
