package dr.evomodel.treedatalikelihood.continuous.backprop;

import org.ejml.data.DenseMatrix64F;

// Specialized Fréchet plan for block-diagonal matrices whose 2x2 blocks all have equal diagonals.
// Avoids the generic Sylvester 4x4 solve on every evaluation.
final class BlockDiagonalFrechetEqualDiagonalPlan {

    private static final ThreadLocal<EqualDiagonalFrechet2x2Kernel.Workspace> KERNEL_WORKSPACE =
            ThreadLocal.withInitial(EqualDiagonalFrechet2x2Kernel.Workspace::new);

    private BlockDiagonalFrechetEqualDiagonalPlan() {
    }

    static Plan createPlan(final BlockDiagonalExpSolver.BlockStructure structure) {
        return new Plan(structure);
    }

    static final class Plan {
        private final int dim;
        private final boolean hasLeadingScalar;
        private final int block2Count;
        private final int[] block2Starts;

        // Per-block parameter caches (one entry per 2x2 block)
        private double scalar;  // only used if hasLeadingScalar
        private final double[] diag;
        private final double[] upper;
        private final double[] lower;
        private final double[] product;
        private final double[] rootMagnitude;
        private final boolean[] productNonNegative;

        // Time-scaled versions (recomputed in evaluate)
        private final double[] scaledDiag;
        private final double[] scaledUpper;
        private final double[] scaledLower;

        // Coefficient storage
        private final double[] blockBlockCoefficients;  // block2Count * block2Count * 16
        private double scalarScalarCoefficient;         // only used if hasLeadingScalar
        private final double[] scalarBlockCoefficients; // block2Count * 4, or null
        private final double[] blockScalarCoefficients; // block2Count * 4, or null

        private Plan(final BlockDiagonalExpSolver.BlockStructure structure) {
            this.dim = structure.getDim();

            final int blockCount = structure.getNumBlocks();
            this.hasLeadingScalar = blockCount > 0 && structure.getBlockSize(0) == 1;

            int count = 0;
            for (int b = 0; b < blockCount; ++b) {
                if (structure.getBlockSize(b) == 2) {
                    count++;
                }
            }
            this.block2Count = count;
            this.block2Starts = new int[count];
            int idx = 0;
            for (int b = 0; b < blockCount; ++b) {
                if (structure.getBlockSize(b) == 2) {
                    block2Starts[idx++] = structure.getBlockStart(b);
                }
            }

            this.diag = new double[block2Count];
            this.upper = new double[block2Count];
            this.lower = new double[block2Count];
            this.product = new double[block2Count];
            this.rootMagnitude = new double[block2Count];
            this.productNonNegative = new boolean[block2Count];

            this.scaledDiag = new double[block2Count];
            this.scaledUpper = new double[block2Count];
            this.scaledLower = new double[block2Count];

            this.blockBlockCoefficients = new double[block2Count * block2Count * 16];

            if (hasLeadingScalar) {
                this.scalarBlockCoefficients = block2Count > 0 ? new double[block2Count * 4] : null;
                this.blockScalarCoefficients = block2Count > 0 ? new double[block2Count * 4] : null;
            } else {
                this.scalarBlockCoefficients = null;
                this.blockScalarCoefficients = null;
            }
        }

        void updateParameters(final double[] blockDParams) {
            final int upperOffset = dim;
            final int lowerOffset = dim + (dim - 1);

            if (hasLeadingScalar) {
                scalar = blockDParams[0];
            }

            for (int k = 0; k < block2Count; ++k) {
                final int s = block2Starts[k];
                diag[k] = blockDParams[s];
                upper[k] = blockDParams[upperOffset + s];
                lower[k] = blockDParams[lowerOffset + s];
                product[k] = upper[k] * lower[k];
                productNonNegative[k] = product[k] >= 0.0;
                rootMagnitude[k] = Math.sqrt(Math.abs(product[k]));
            }
        }

        void evaluate(final double t) {
            final double scale = -t;
            final double tAbs = Math.abs(t);
            final double t2 = t * t;

            for (int k = 0; k < block2Count; ++k) {
                scaledDiag[k] = scale * diag[k];
                scaledUpper[k] = scale * upper[k];
                scaledLower[k] = scale * lower[k];
            }

            final EqualDiagonalFrechet2x2Kernel.Workspace workspace = KERNEL_WORKSPACE.get();

            for (int left = 0; left < block2Count; ++left) {
                final double a = scaledDiag[left];
                final double scaledLeftProduct = t2 * product[left];
                final double rAbs = tAbs * rootMagnitude[left];

                for (int right = 0; right < block2Count; ++right) {
                    final double b = scaledDiag[right];
                    final double scaledRightProduct = t2 * product[right];
                    final double qAbs = tAbs * rootMagnitude[right];
                    final int coeffOffset = 16 * (left * block2Count + right);

                    EqualDiagonalFrechet2x2Kernel.fillPrepared(
                            a, b,
                            scaledLeftProduct, scaledRightProduct,
                            productNonNegative[left], productNonNegative[right],
                            rAbs, qAbs,
                            workspace.prepared,
                            workspace);
                    EqualDiagonalFrechet2x2Kernel.fillCoefficientMatrix(
                            workspace.prepared,
                            scaledUpper[left], scaledLower[left],
                            scaledUpper[right], scaledLower[right],
                            blockBlockCoefficients, coeffOffset);
                }
            }

            if (hasLeadingScalar) {
                final double scaledScalar = scale * scalar;
                // same eigenvalue on both sides — divided difference collapses to exp
                scalarScalarCoefficient = Math.exp(scaledScalar);

                for (int k = 0; k < block2Count; ++k) {
                    EqualDiagonalFrechet2x2Kernel.fillOneByTwoCoefficients(
                            scaledScalar, scaledDiag[k], scaledUpper[k], scaledLower[k],
                            scalarBlockCoefficients, 4 * k, workspace);
                    EqualDiagonalFrechet2x2Kernel.fillTwoByOneCoefficients(
                            scaledDiag[k], scaledUpper[k], scaledLower[k], scaledScalar,
                            blockScalarCoefficients, 4 * k, workspace);
                }
            }
        }

        void apply(final DenseMatrix64F input,
                   final DenseMatrix64F out) {
            final double[] inData = input.data;
            final double[] outData = out.data;

            // Block-block pairs
            for (int left = 0; left < block2Count; ++left) {
                final int leftStart = block2Starts[left];
                for (int right = 0; right < block2Count; ++right) {
                    final int rightStart = block2Starts[right];
                    final int base = leftStart * dim + rightStart;
                    final double e00 = inData[base];
                    final double e01 = inData[base + 1];
                    final double e10 = inData[base + dim];
                    final double e11 = inData[base + dim + 1];
                    final int c = 16 * (left * block2Count + right);
                    final double[] cc = blockBlockCoefficients;
                    outData[base]               = cc[c]      * e00 + cc[c + 1]  * e01 + cc[c + 2]  * e10 + cc[c + 3]  * e11;
                    outData[base + 1]           = cc[c + 4]  * e00 + cc[c + 5]  * e01 + cc[c + 6]  * e10 + cc[c + 7]  * e11;
                    outData[base + dim]         = cc[c + 8]  * e00 + cc[c + 9]  * e01 + cc[c + 10] * e10 + cc[c + 11] * e11;
                    outData[base + dim + 1]     = cc[c + 12] * e00 + cc[c + 13] * e01 + cc[c + 14] * e10 + cc[c + 15] * e11;
                }
            }

            if (hasLeadingScalar) {
                // Scalar-scalar: row 0, col 0
                outData[0] = scalarScalarCoefficient * inData[0];

                // Scalar-block: row 0, each 2x2 block column
                for (int k = 0; k < block2Count; ++k) {
                    final int rightStart = block2Starts[k];
                    final int base = rightStart;  // row 0 * dim + rightStart
                    final double e0 = inData[base];
                    final double e1 = inData[base + 1];
                    final double[] sc = scalarBlockCoefficients;
                    final int c = 4 * k;
                    outData[base]     = sc[c]     * e0 + sc[c + 1] * e1;
                    outData[base + 1] = sc[c + 2] * e0 + sc[c + 3] * e1;
                }

                // Block-scalar: each 2x2 block row, col 0
                for (int k = 0; k < block2Count; ++k) {
                    final int leftStart = block2Starts[k];
                    final int base = leftStart * dim;  // leftStart*dim + col 0
                    final double e0 = inData[base];
                    final double e1 = inData[base + dim];
                    final double[] bc = blockScalarCoefficients;
                    final int c = 4 * k;
                    outData[base]       = bc[c]     * e0 + bc[c + 1] * e1;
                    outData[base + dim] = bc[c + 2] * e0 + bc[c + 3] * e1;
                }
            }
        }
    }
}
