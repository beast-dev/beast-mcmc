package dr.inference.timeseries.engine.gaussian;

import dr.inference.timeseries.representation.CanonicalGaussianTransition;

/**
 * Direct pullback from canonical transition blocks to moment-form transition
 * adjoints.
 */
public final class CanonicalTransitionAdjointUtils {

    public static final class Workspace {
        final double[][] omega;
        final double[][] transitionMatrix;
        final double[] transitionOffset;
        final double[][] gXxPlusTranspose;
        final double[][] temp1;
        final double[][] temp2;
        final double[][] gP;
        final double[] tempv;

        public Workspace(final int dimension) {
            this.omega = new double[dimension][dimension];
            this.transitionMatrix = new double[dimension][dimension];
            this.transitionOffset = new double[dimension];
            this.gXxPlusTranspose = new double[dimension][dimension];
            this.temp1 = new double[dimension][dimension];
            this.temp2 = new double[dimension][dimension];
            this.gP = new double[dimension][dimension];
            this.tempv = new double[dimension];
        }
    }

    private CanonicalTransitionAdjointUtils() {
        // no instances
    }

    public static void fillFromCanonicalTransition(final CanonicalGaussianTransition transition,
                                                   final CanonicalBranchMessageContribution contribution,
                                                   final CanonicalLocalTransitionAdjoints out) {
        fillFromCanonicalTransition(transition, contribution, new Workspace(transition.getDimension()), out);
    }

    public static void fillFromCanonicalTransition(final CanonicalGaussianTransition transition,
                                                   final CanonicalBranchMessageContribution contribution,
                                                   final Workspace workspace,
                                                   final CanonicalLocalTransitionAdjoints out) {
        final int d = transition.getDimension();
        ensureWorkspaceDimension(workspace, d);
        invertPositiveDefinite(transition.precisionYY, workspace.omega);

        multiply(workspace.omega, transition.precisionYX, workspace.transitionMatrix);
        scaleInPlace(workspace.transitionMatrix, -1.0);

        multiply(workspace.omega, transition.informationY, workspace.transitionOffset);

        fillFromMoments(
                transition.precisionYY,
                workspace.omega,
                workspace.transitionMatrix,
                workspace.transitionOffset,
                contribution,
                workspace,
                out);
    }

    public static void fillFromMoments(final double[][] precision,
                                       final double[][] omega,
                                       final double[][] transitionMatrix,
                                       final double[] transitionOffset,
                                       final CanonicalBranchMessageContribution contribution,
                                       final CanonicalLocalTransitionAdjoints out) {
        fillFromMoments(
                precision,
                omega,
                transitionMatrix,
                transitionOffset,
                contribution,
                new Workspace(transitionOffset.length),
                out);
    }

    public static void fillFromMoments(final double[][] precision,
                                       final double[][] omega,
                                       final double[][] transitionMatrix,
                                       final double[] transitionOffset,
                                       final CanonicalBranchMessageContribution contribution,
                                       final Workspace workspace,
                                       final CanonicalLocalTransitionAdjoints out) {
        final int d = transitionOffset.length;
        ensureWorkspaceDimension(workspace, d);
        final double[][] gXx = contribution.dLogL_dPrecisionXX;
        final double[][] gXy = contribution.dLogL_dPrecisionXY;
        final double[][] gYx = contribution.dLogL_dPrecisionYX;
        final double[][] gYy = contribution.dLogL_dPrecisionYY;
        final double[] gX = contribution.dLogL_dInformationX;
        final double[] gY = contribution.dLogL_dInformationY;
        final double g0 = contribution.dLogL_dLogNormalizer;

        final double[][] gXxPlusTranspose = workspace.gXxPlusTranspose;
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                gXxPlusTranspose[i][j] = gXx[i][j] + gXx[j][i];
            }
        }

        final double[][] temp1 = workspace.temp1;
        final double[][] temp2 = workspace.temp2;
        multiply(transitionMatrix, gXxPlusTranspose, temp1);
        multiply(precision, temp1, out.dLogL_dF);

        transposeInto(gXy, temp1);
        multiply(precision, temp1, temp2);
        subtractInPlace(out.dLogL_dF, temp2);

        multiply(precision, gYx, temp2);
        subtractInPlace(out.dLogL_dF, temp2);

        outerProduct(transitionOffset, gX, temp1);
        multiply(precision, temp1, temp2);
        subtractInPlace(out.dLogL_dF, temp2);

        final double[] tempv = workspace.tempv;
        multiply(transitionMatrix, gX, tempv);
        for (int i = 0; i < d; ++i) {
            tempv[i] = -tempv[i] + gY[i] + g0 * transitionOffset[i];
        }
        multiply(precision, tempv, out.dLogL_df);

        final double[][] gP = workspace.gP;
        multiply(transitionMatrix, gXx, temp1);
        multiplyTransposedRight(temp1, transitionMatrix, gP);

        multiply(transitionMatrix, gXy, temp1);
        subtractInPlace(gP, temp1);

        multiplyTransposedRight(gYx, transitionMatrix, temp1);
        subtractInPlace(gP, temp1);

        addInPlace(gP, gYy);

        outerProduct(gX, transitionOffset, temp1);
        multiply(transitionMatrix, temp1, temp2);
        subtractInPlace(gP, temp2);

        outerProduct(gY, transitionOffset, temp1);
        addInPlace(gP, temp1);

        outerProduct(transitionOffset, transitionOffset, temp1);
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                gP[i][j] += 0.5 * g0 * temp1[i][j];
                gP[i][j] -= 0.5 * g0 * omega[i][j];
            }
        }

        multiply(precision, gP, temp1);
        multiply(temp1, precision, out.dLogL_dOmega);
        scaleInPlace(out.dLogL_dOmega, -1.0);
        symmetrize(out.dLogL_dOmega);
    }

    private static void ensureWorkspaceDimension(final Workspace workspace, final int dimension) {
        if (workspace.omega.length != dimension) {
            throw new IllegalArgumentException("Workspace dimension mismatch");
        }
    }

    private static void invertPositiveDefinite(final double[][] matrix, final double[][] out) {
        final int d = matrix.length;
        final double[][] work = new double[d][d];
        for (int i = 0; i < d; ++i) {
            System.arraycopy(matrix[i], 0, work[i], 0, d);
        }
        final double[][] chol = new double[d][d];
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = work[i][j];
                for (int k = 0; k < j; ++k) {
                    sum -= chol[i][k] * chol[j][k];
                }
                if (i == j) {
                    chol[i][j] = Math.sqrt(sum);
                } else {
                    chol[i][j] = sum / chol[j][j];
                }
            }
        }
        final double[][] lowerInverse = new double[d][d];
        for (int col = 0; col < d; ++col) {
            for (int row = 0; row < d; ++row) {
                double sum = (row == col) ? 1.0 : 0.0;
                for (int k = 0; k < row; ++k) {
                    sum -= chol[row][k] * lowerInverse[k][col];
                }
                lowerInverse[row][col] = sum / chol[row][row];
            }
        }
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += lowerInverse[k][i] * lowerInverse[k][j];
                }
                out[i][j] = sum;
            }
        }
        symmetrize(out);
    }

    private static void multiply(final double[][] left, final double[][] right, final double[][] out) {
        for (int i = 0; i < left.length; ++i) {
            for (int j = 0; j < right[0].length; ++j) {
                double sum = 0.0;
                for (int k = 0; k < right.length; ++k) {
                    sum += left[i][k] * right[k][j];
                }
                out[i][j] = sum;
            }
        }
    }

    private static void multiply(final double[][] matrix, final double[] vector, final double[] out) {
        for (int i = 0; i < matrix.length; ++i) {
            double sum = 0.0;
            for (int j = 0; j < vector.length; ++j) {
                sum += matrix[i][j] * vector[j];
            }
            out[i] = sum;
        }
    }

    private static void multiplyTransposedRight(final double[][] left, final double[][] right, final double[][] out) {
        for (int i = 0; i < left.length; ++i) {
            for (int j = 0; j < right.length; ++j) {
                double sum = 0.0;
                for (int k = 0; k < left[0].length; ++k) {
                    sum += left[i][k] * right[j][k];
                }
                out[i][j] = sum;
            }
        }
    }

    private static void outerProduct(final double[] left, final double[] right, final double[][] out) {
        for (int i = 0; i < left.length; ++i) {
            for (int j = 0; j < right.length; ++j) {
                out[i][j] = left[i] * right[j];
            }
        }
    }

    private static void transposeInto(final double[][] source, final double[][] out) {
        for (int i = 0; i < source.length; ++i) {
            for (int j = 0; j < source[i].length; ++j) {
                out[j][i] = source[i][j];
            }
        }
    }

    private static void addInPlace(final double[][] target, final double[][] delta) {
        for (int i = 0; i < target.length; ++i) {
            for (int j = 0; j < target[i].length; ++j) {
                target[i][j] += delta[i][j];
            }
        }
    }

    private static void subtractInPlace(final double[][] target, final double[][] delta) {
        for (int i = 0; i < target.length; ++i) {
            for (int j = 0; j < target[i].length; ++j) {
                target[i][j] -= delta[i][j];
            }
        }
    }

    private static void scaleInPlace(final double[][] matrix, final double factor) {
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < matrix[i].length; ++j) {
                matrix[i][j] *= factor;
            }
        }
    }

    private static void symmetrize(final double[][] matrix) {
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = i + 1; j < matrix.length; ++j) {
                final double average = 0.5 * (matrix[i][j] + matrix[j][i]);
                matrix[i][j] = average;
                matrix[j][i] = average;
            }
        }
    }
}
