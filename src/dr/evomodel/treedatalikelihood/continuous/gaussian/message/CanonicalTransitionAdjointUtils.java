package dr.evomodel.treedatalikelihood.continuous.gaussian.message;

import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianTransition;

/**
 * Direct pullback from canonical transition blocks to moment-form transition
 * adjoints.
 *
 * <p>All matrix fields use flat row-major storage: element {@code (i, j)} of
 * a {@code dim × dim} matrix is at index {@code i * dim + j}. Callers are
 * expected to retain and reuse a {@link Workspace}; the hot-path methods in
 * this helper do not allocate.
 */
public final class CanonicalTransitionAdjointUtils {

    public static final class Workspace {
        final double[] omega;
        public final double[] transitionMatrix;
        final double[] transitionOffset;
        final double[] gXxPlusTranspose;
        final double[] temp1;
        final double[] temp2;
        final double[] gP;
        final double[] tempv;
        private final int dimension;

        public Workspace(final int dimension) {
            this.dimension = dimension;
            final int d2 = dimension * dimension;
            this.omega             = new double[d2];
            this.transitionMatrix  = new double[d2];
            this.transitionOffset  = new double[dimension];
            this.gXxPlusTranspose  = new double[d2];
            this.temp1             = new double[d2];
            this.temp2             = new double[d2];
            this.gP                = new double[d2];
            this.tempv             = new double[dimension];
        }

        public int getDimension() { return dimension; }
    }

    private CanonicalTransitionAdjointUtils() {
        // no instances
    }

    public static void fillFromCanonicalTransition(final CanonicalGaussianTransition transition,
                                                   final CanonicalBranchMessageContribution contribution,
                                                   final Workspace workspace,
                                                   final CanonicalLocalTransitionAdjoints out) {
        final int d = transition.getDimension();
        ensureWorkspaceDimension(workspace, d);
        invertPositiveDefiniteFlat(transition.precisionYY, workspace, workspace.omega, d);

        multiplyFlat(workspace.omega, transition.precisionYX, workspace.transitionMatrix, d);
        scaleInPlaceFlat(workspace.transitionMatrix, -1.0, d * d);

        multiplyMatVecFlat(workspace.omega, transition.informationY, workspace.transitionOffset, d);

        fillFromMomentsFlat(
                transition.precisionYY,
                workspace.omega,
                workspace.transitionMatrix,
                workspace.transitionOffset,
                contribution,
                workspace,
                out);
    }

    public static void fillFromMoments(final double[] precision,
                                       final double[] omega,
                                       final double[] transitionMatrix,
                                       final double[] transitionOffset,
                                       final CanonicalBranchMessageContribution contribution,
                                       final Workspace workspace,
                                       final CanonicalLocalTransitionAdjoints out) {
        fillFromMomentsFlat(precision, omega, transitionMatrix, transitionOffset,
                contribution, workspace, out);
    }

    private static void fillFromMomentsFlat(final double[] precision,
                                            final double[] omega,
                                            final double[] transitionMatrix,
                                            final double[] transitionOffset,
                                            final CanonicalBranchMessageContribution contribution,
                                            final Workspace workspace,
                                            final CanonicalLocalTransitionAdjoints out) {
        final int d = transitionOffset.length;
        ensureWorkspaceDimension(workspace, d);
        final int d2 = d * d;

        final double[] gXx = contribution.dLogL_dPrecisionXX;
        final double[] gXy = contribution.dLogL_dPrecisionXY;
        final double[] gYx = contribution.dLogL_dPrecisionYX;
        final double[] gYy = contribution.dLogL_dPrecisionYY;
        final double[] gX  = contribution.dLogL_dInformationX;
        final double[] gY  = contribution.dLogL_dInformationY;
        final double   g0  = contribution.dLogL_dLogNormalizer;

        final double[] gXxPlusTranspose = workspace.gXxPlusTranspose;
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                gXxPlusTranspose[i * d + j] = gXx[i * d + j] + gXx[j * d + i];
            }
        }

        final double[] temp1 = workspace.temp1;
        final double[] temp2 = workspace.temp2;

        // dLogL_dF: precision × (transitionMatrix × (gXx + gXx^T))^T - precision × gXy^T - precision × gYx - precision × (f ⊗ gX)^T
        multiplyFlat(transitionMatrix, gXxPlusTranspose, temp1, d);
        multiplyFlat(precision, temp1, out.dLogL_dF, d);

        transposeIntoFlat(gXy, temp1, d);
        multiplyFlat(precision, temp1, temp2, d);
        subtractInPlaceFlat(out.dLogL_dF, temp2, d2);

        multiplyFlat(precision, gYx, temp2, d);
        subtractInPlaceFlat(out.dLogL_dF, temp2, d2);

        outerProductFlat(transitionOffset, gX, temp1, d);
        multiplyFlat(precision, temp1, temp2, d);
        subtractInPlaceFlat(out.dLogL_dF, temp2, d2);

        final double[] tempv = workspace.tempv;
        multiplyMatVecFlat(transitionMatrix, gX, tempv, d);
        for (int i = 0; i < d; ++i) {
            tempv[i] = -tempv[i] + gY[i] + g0 * transitionOffset[i];
        }
        multiplyMatVecFlat(precision, tempv, out.dLogL_df, d);

        final double[] gP = workspace.gP;
        multiplyFlat(transitionMatrix, gXx, temp1, d);
        multiplyTransposedRightFlat(temp1, transitionMatrix, gP, d);

        multiplyFlat(transitionMatrix, gXy, temp1, d);
        subtractInPlaceFlat(gP, temp1, d2);

        multiplyTransposedRightFlat(gYx, transitionMatrix, temp1, d);
        subtractInPlaceFlat(gP, temp1, d2);

        addInPlaceFlat(gP, gYy, d2);

        outerProductFlat(gX, transitionOffset, temp1, d);
        multiplyFlat(transitionMatrix, temp1, temp2, d);
        subtractInPlaceFlat(gP, temp2, d2);

        outerProductFlat(gY, transitionOffset, temp1, d);
        addInPlaceFlat(gP, temp1, d2);

        outerProductFlat(transitionOffset, transitionOffset, temp1, d);
        for (int k = 0; k < d2; ++k) {
            gP[k] += 0.5 * g0 * (temp1[k] - omega[k]);
        }

        multiplyFlat(precision, gP, temp1, d);
        multiplyFlat(temp1, precision, out.dLogL_dOmega, d);
        scaleInPlaceFlat(out.dLogL_dOmega, -1.0, d2);
        GaussianMatrixOps.symmetrizeFlat(out.dLogL_dOmega, d);
    }

    private static void ensureWorkspaceDimension(final Workspace workspace, final int dimension) {
        if (workspace.getDimension() != dimension) {
            throw new IllegalArgumentException("Workspace dimension mismatch");
        }
    }

    private static void invertPositiveDefiniteFlat(final double[] matrix,
                                                   final Workspace workspace,
                                                   final double[] out,
                                                   final int d) {
        final double[] chol = workspace.temp1;
        final double[] lowerInverse = workspace.temp2;
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = matrix[i * d + j];
                for (int k = 0; k < j; ++k) {
                    sum -= chol[i * d + k] * chol[j * d + k];
                }
                if (i == j) {
                    chol[i * d + j] = Math.sqrt(sum);
                } else {
                    chol[i * d + j] = sum / chol[j * d + j];
                }
            }
        }
        for (int col = 0; col < d; ++col) {
            for (int row = 0; row < d; ++row) {
                double sum = (row == col) ? 1.0 : 0.0;
                for (int k = 0; k < row; ++k) {
                    sum -= chol[row * d + k] * lowerInverse[k * d + col];
                }
                lowerInverse[row * d + col] = sum / chol[row * d + row];
            }
        }
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += lowerInverse[k * d + i] * lowerInverse[k * d + j];
                }
                out[i * d + j] = sum;
            }
        }
        GaussianMatrixOps.symmetrizeFlat(out, d);
    }

    private static void multiplyFlat(final double[] left, final double[] right,
                                     final double[] out, final int d) {
        GaussianMatrixOps.multiplyMatrixMatrixFlat(left, right, out, d);
    }

    private static void multiplyMatVecFlat(final double[] matrix, final double[] vector,
                                           final double[] out, final int d) {
        GaussianMatrixOps.multiplyMatrixVectorFlat(matrix, vector, out, d);
    }

    private static void multiplyTransposedRightFlat(final double[] left, final double[] right,
                                                    final double[] out, final int d) {
        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += left[iOff + k] * right[j * d + k];
                }
                out[iOff + j] = sum;
            }
        }
    }

    private static void outerProductFlat(final double[] left, final double[] right,
                                         final double[] out, final int d) {
        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                out[iOff + j] = left[i] * right[j];
            }
        }
    }

    private static void transposeIntoFlat(final double[] source, final double[] out, final int d) {
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                out[j * d + i] = source[i * d + j];
            }
        }
    }

    private static void addInPlaceFlat(final double[] target, final double[] delta, final int d2) {
        for (int k = 0; k < d2; ++k) {
            target[k] += delta[k];
        }
    }

    private static void subtractInPlaceFlat(final double[] target, final double[] delta, final int d2) {
        for (int k = 0; k < d2; ++k) {
            target[k] -= delta[k];
        }
    }

    private static void scaleInPlaceFlat(final double[] matrix, final double factor, final int d2) {
        for (int k = 0; k < d2; ++k) {
            matrix[k] *= factor;
        }
    }
}
