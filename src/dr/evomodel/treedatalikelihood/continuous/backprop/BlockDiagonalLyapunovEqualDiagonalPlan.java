package dr.evomodel.treedatalikelihood.continuous.backprop;

import org.ejml.data.DenseMatrix64F;

/**
 * Cached equal-diagonal Lyapunov plan for
 *
 *     D * X + X * D^T = H
 *
 * with an optional leading scalar block followed by equal-diagonal 2x2 blocks.
 * The hot apply path assumes symmetric H, solves only upper block pairs, and
 * mirrors the result.  No allocation occurs after construction.
 */
final class BlockDiagonalLyapunovEqualDiagonalPlan {

    private static final double EPS = 1.0e-12;

    private BlockDiagonalLyapunovEqualDiagonalPlan() {
    }

    static Plan createPlan(final BlockDiagonalExpSolver.BlockStructure structure) {
        return new Plan(structure);
    }

    static final class Plan {
        private final int dim;
        private final int upperOffset;
        private final int lowerOffset;
        private final boolean hasLeadingScalar;
        private final int block2Count;
        private final int[] block2Starts;

        private final double[] diag;
        private final double[] upper;
        private final double[] lower;

        // Row-major coefficient maps:
        // scalarBlock: [y0,y1] <- 2x2 * [h0,h1]
        // diagonalBlock: [y00,y01,y11] <- 3x3 * [h00,h01,h11]
        // offDiagonalBlock: [y00,y01,y10,y11] <- 4x4 * [h00,h01,h10,h11]
        private final double[] scalarBlockCoefficients;
        private final double[] diagonalBlockCoefficients;
        private final double[] offDiagonalBlockCoefficients;

        private double scalar;
        private double scalarScalarCoefficient;
        private boolean valid;

        private Plan(final BlockDiagonalExpSolver.BlockStructure structure) {
            this.dim = structure.getDim();
            this.upperOffset = dim;
            this.lowerOffset = dim + dim - 1;

            final int blockCount = structure.getNumBlocks();
            this.hasLeadingScalar = blockCount > 0 && structure.getBlockSize(0) == 1;

            int count = 0;
            int firstTwoBlock = 0;
            if (hasLeadingScalar) {
                firstTwoBlock = 1;
            }
            for (int b = firstTwoBlock; b < blockCount; ++b) {
                if (structure.getBlockSize(b) != 2) {
                    throw new IllegalArgumentException(
                            "Equal-diagonal Lyapunov plan supports only an optional leading scalar followed by 2x2 blocks");
                }
                count++;
            }

            this.block2Count = count;
            this.block2Starts = new int[count];
            for (int b = firstTwoBlock, k = 0; b < blockCount; ++b, ++k) {
                block2Starts[k] = structure.getBlockStart(b);
            }

            this.diag = new double[count];
            this.upper = new double[count];
            this.lower = new double[count];
            this.scalarBlockCoefficients = hasLeadingScalar ? new double[4 * count] : null;
            this.diagonalBlockCoefficients = new double[9 * count];
            this.offDiagonalBlockCoefficients = new double[16 * count * count];
            this.valid = false;
        }

        boolean updateParameters(final double[] blockDParams) {
            return updateParameters(blockDParams, false);
        }

        boolean updateTransposedParameters(final double[] blockDParams) {
            return updateParameters(blockDParams, true);
        }

        private boolean updateParameters(final double[] blockDParams,
                                         final boolean transpose) {
            if (blockDParams.length != 3 * dim - 2) {
                throw new IllegalArgumentException("Expected compressed block matrix length " + (3 * dim - 2));
            }

            if (hasLeadingScalar) {
                scalar = blockDParams[0];
                final double denom = 2.0 * scalar;
                if (!isNonSingular(denom)) {
                    valid = false;
                    return false;
                }
                scalarScalarCoefficient = 1.0 / denom;
            }

            for (int k = 0; k < block2Count; ++k) {
                final int start = block2Starts[k];
                final double a = blockDParams[start];
                final double u = blockDParams[(transpose ? lowerOffset : upperOffset) + start];
                final double v = blockDParams[(transpose ? upperOffset : lowerOffset) + start];

                if (Math.abs(a - blockDParams[start + 1]) > EPS) {
                    valid = false;
                    return false;
                }

                diag[k] = a;
                upper[k] = u;
                lower[k] = v;

                if (!fillDiagonalBlockCoefficients(a, u, v, diagonalBlockCoefficients, 9 * k)) {
                    valid = false;
                    return false;
                }

                if (hasLeadingScalar &&
                        !fillOneTwoCoefficients(scalar, a, u, v, scalarBlockCoefficients, 4 * k)) {
                    valid = false;
                    return false;
                }
            }

            for (int left = 0; left < block2Count; ++left) {
                final double ai = diag[left];
                final double ui = upper[left];
                final double vi = lower[left];

                for (int right = left + 1; right < block2Count; ++right) {
                    if (!fillOffDiagonalBlockCoefficients(
                            ai, ui, vi,
                            diag[right], upper[right], lower[right],
                            offDiagonalBlockCoefficients,
                            offDiagonalOffset(left, right))) {
                        valid = false;
                        return false;
                    }
                }
            }

            valid = true;
            return true;
        }

        boolean tryApplySymmetric(final DenseMatrix64F rhs,
                                  final DenseMatrix64F out) {
            checkSquare(rhs, "rhs");
            checkSquare(out, "out");
            return tryApplySymmetric(rhs.data, out.data);
        }

        boolean tryApplySymmetric(final double[] rhs,
                                  final double[] out) {
            if (!valid) {
                return false;
            }

            final int n = dim;

            if (hasLeadingScalar) {
                out[0] = scalarScalarCoefficient * rhs[0];

                final double[] sc = scalarBlockCoefficients;
                for (int k = 0; k < block2Count; ++k) {
                    final int j = block2Starts[k];
                    final int upperBase = j;
                    final double h0 = rhs[upperBase];
                    final double h1 = rhs[upperBase + 1];
                    final int c = 4 * k;
                    final double y0 = sc[c] * h0 + sc[c + 1] * h1;
                    final double y1 = sc[c + 2] * h0 + sc[c + 3] * h1;

                    out[upperBase] = y0;
                    out[upperBase + 1] = y1;
                    out[j * n] = y0;
                    out[(j + 1) * n] = y1;
                }
            }

            final double[] dc = diagonalBlockCoefficients;
            final double[] oc = offDiagonalBlockCoefficients;

            for (int left = 0; left < block2Count; ++left) {
                final int i = block2Starts[left];
                final int ii = i * n + i;

                final double h00 = rhs[ii];
                final double h01 = rhs[ii + 1];
                final double h11 = rhs[ii + n + 1];
                int c = 9 * left;

                final double y00 = dc[c] * h00 + dc[c + 1] * h01 + dc[c + 2] * h11;
                final double y01 = dc[c + 3] * h00 + dc[c + 4] * h01 + dc[c + 5] * h11;
                final double y11 = dc[c + 6] * h00 + dc[c + 7] * h01 + dc[c + 8] * h11;

                out[ii] = y00;
                out[ii + 1] = y01;
                out[ii + n] = y01;
                out[ii + n + 1] = y11;

                for (int right = left + 1; right < block2Count; ++right) {
                    final int j = block2Starts[right];
                    final int ij = i * n + j;

                    final double e00 = rhs[ij];
                    final double e01 = rhs[ij + 1];
                    final double e10 = rhs[ij + n];
                    final double e11 = rhs[ij + n + 1];

                    c = offDiagonalOffset(left, right);
                    final double z00 = oc[c] * e00 + oc[c + 1] * e01 + oc[c + 2] * e10 + oc[c + 3] * e11;
                    final double z01 = oc[c + 4] * e00 + oc[c + 5] * e01 + oc[c + 6] * e10 + oc[c + 7] * e11;
                    final double z10 = oc[c + 8] * e00 + oc[c + 9] * e01 + oc[c + 10] * e10 + oc[c + 11] * e11;
                    final double z11 = oc[c + 12] * e00 + oc[c + 13] * e01 + oc[c + 14] * e10 + oc[c + 15] * e11;

                    out[ij] = z00;
                    out[ij + 1] = z01;
                    out[ij + n] = z10;
                    out[ij + n + 1] = z11;

                    final int ji = j * n + i;
                    out[ji] = z00;
                    out[ji + 1] = z10;
                    out[ji + n] = z01;
                    out[ji + n + 1] = z11;
                }
            }

            return true;
        }

        private int offDiagonalOffset(final int left,
                                      final int right) {
            return 16 * (left * block2Count + right);
        }

        private void checkSquare(final DenseMatrix64F matrix,
                                 final String name) {
            if (matrix.numRows != dim || matrix.numCols != dim) {
                throw new IllegalArgumentException(name + " must be " + dim + " x " + dim);
            }
        }
    }

    private static boolean fillOneTwoCoefficients(final double scalar,
                                                  final double a,
                                                  final double u,
                                                  final double v,
                                                  final double[] out,
                                                  final int offset) {
        final double s = scalar + a;
        final double det = s * s - u * v;
        if (!isNonSingular(det)) {
            return false;
        }
        final double inv = 1.0 / det;
        out[offset] = s * inv;
        out[offset + 1] = -u * inv;
        out[offset + 2] = -v * inv;
        out[offset + 3] = s * inv;
        return true;
    }

    private static boolean fillDiagonalBlockCoefficients(final double a,
                                                         final double u,
                                                         final double v,
                                                         final double[] out,
                                                         final int offset) {
        final double a11 = 2.0 * a;
        final double a12 = 2.0 * u;
        final double a21 = v;
        final double a22 = 2.0 * a;
        final double a23 = u;
        final double a32 = 2.0 * v;
        final double a33 = 2.0 * a;

        final double c00 = a22 * a33 - a23 * a32;
        final double c01 = -a21 * a33;
        final double c02 = a21 * a32;
        final double c10 = -a12 * a33;
        final double c11 = a11 * a33;
        final double c12 = -a11 * a32;
        final double c20 = a12 * a23;
        final double c21 = -a11 * a23;
        final double c22 = a11 * a22 - a12 * a21;

        final double det = a11 * c00 + a12 * c01;
        if (!isNonSingular(det)) {
            return false;
        }
        final double inv = 1.0 / det;

        out[offset] = c00 * inv;
        out[offset + 1] = c10 * inv;
        out[offset + 2] = c20 * inv;
        out[offset + 3] = c01 * inv;
        out[offset + 4] = c11 * inv;
        out[offset + 5] = c21 * inv;
        out[offset + 6] = c02 * inv;
        out[offset + 7] = c12 * inv;
        out[offset + 8] = c22 * inv;
        return true;
    }

    private static boolean fillOffDiagonalBlockCoefficients(final double ai,
                                                            final double ui,
                                                            final double vi,
                                                            final double aj,
                                                            final double uj,
                                                            final double vj,
                                                            final double[] out,
                                                            final int offset) {
        final double s = ai + aj;
        final double c1 = s * s + uj * vj - ui * vi;
        final double c2 = 2.0 * s * uj;
        final double c3 = 2.0 * s * vj;
        final double det = c1 * c1 - c2 * c3;

        if (!isNonSingular(det)) {
            return false;
        }
        final double inv = 1.0 / det;

        out[offset] = (c1 * s - c2 * vj) * inv;
        out[offset + 1] = (c1 * uj - c2 * s) * inv;
        out[offset + 2] = (-c1 * ui) * inv;
        out[offset + 3] = (c2 * ui) * inv;

        out[offset + 4] = (-c3 * s + c1 * vj) * inv;
        out[offset + 5] = (-c3 * uj + c1 * s) * inv;
        out[offset + 6] = (c3 * ui) * inv;
        out[offset + 7] = (-c1 * ui) * inv;

        out[offset + 8] = (-c1 * vi) * inv;
        out[offset + 9] = (c2 * vi) * inv;
        out[offset + 10] = (c1 * s - c2 * vj) * inv;
        out[offset + 11] = (c1 * uj - c2 * s) * inv;

        out[offset + 12] = (c3 * vi) * inv;
        out[offset + 13] = (-c1 * vi) * inv;
        out[offset + 14] = (-c3 * s + c1 * vj) * inv;
        out[offset + 15] = (-c3 * uj + c1 * s) * inv;
        return true;
    }

    private static boolean isNonSingular(final double value) {
        return Math.abs(value) >= EPS;
    }
}
