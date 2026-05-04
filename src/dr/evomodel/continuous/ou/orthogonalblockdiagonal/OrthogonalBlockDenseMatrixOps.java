package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import org.ejml.data.DenseMatrix64F;

final class OrthogonalBlockDenseMatrixOps {

    private OrthogonalBlockDenseMatrixOps() {
    }

    static void mult(final DenseMatrix64F left,
                     final DenseMatrix64F right,
                     final DenseMatrix64F out) {
        MatrixOps.matMul(left.data, right.data, out.data, left.numRows);
    }

    static void multSymmetricLeft(final DenseMatrix64F symmetricLeft,
                                  final DenseMatrix64F right,
                                  final DenseMatrix64F out) {
        MatrixOps.multiplySymmetricLeft(symmetricLeft.data, right.data, out.data, symmetricLeft.numRows);
    }

    static void multSymmetricRight(final DenseMatrix64F left,
                                   final DenseMatrix64F symmetricRight,
                                   final DenseMatrix64F out) {
        MatrixOps.multiplySymmetricRight(left.data, symmetricRight.data, out.data, left.numRows);
    }

    static void multTransA(final DenseMatrix64F left,
                           final DenseMatrix64F right,
                           final DenseMatrix64F out) {
        final int dimension = left.numRows;
        final double[] leftData = left.data;
        final double[] rightData = right.data;
        final double[] outData = out.data;
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += leftData[k * dimension + i] * rightData[k * dimension + j];
                }
                outData[rowOffset + j] = sum;
            }
        }
    }

    static void subtract(final DenseMatrix64F left,
                         final DenseMatrix64F right,
                         final DenseMatrix64F out) {
        final int length = left.numRows * left.numCols;
        final double[] leftData = left.data;
        final double[] rightData = right.data;
        final double[] outData = out.data;
        for (int i = 0; i < length; ++i) {
            outData[i] = leftData[i] - rightData[i];
        }
    }

    static void addEquals(final DenseMatrix64F target,
                          final DenseMatrix64F increment) {
        final int length = target.numRows * target.numCols;
        final double[] targetData = target.data;
        final double[] incrementData = increment.data;
        for (int i = 0; i < length; ++i) {
            targetData[i] += incrementData[i];
        }
    }

    static void addEquals(final DenseMatrix64F target,
                          final double scale,
                          final DenseMatrix64F increment) {
        final int length = target.numRows * target.numCols;
        final double[] targetData = target.data;
        final double[] incrementData = increment.data;
        for (int i = 0; i < length; ++i) {
            targetData[i] += scale * incrementData[i];
        }
    }
}
