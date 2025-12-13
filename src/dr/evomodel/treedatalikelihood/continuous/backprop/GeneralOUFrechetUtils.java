package dr.evomodel.treedatalikelihood.continuous.backprop;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import dr.evomodel.treedatalikelihood.continuous.MatrixExponential;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

public class GeneralOUFrechetUtils {

    private GeneralOUFrechetUtils() {
        // no instantiation
    }
    public static void frechetAdjExp(DenseMatrix64F A, DenseMatrix64F X,
                                          DenseMatrix64F out) {
        /**
         * Compute the adjoint of the Fréchet derivative of matrix exponential.
         *
         * This is the adjoint operator that satisfies:
         *   <X, frechetExpBlock(A, E)> = <frechetAdjExpBlock(A, X), E>
         *
         * where <·,·> is the Frobenius inner product: <A,B> = tr(A^T·B) = sum_ij A[i,j]·B[i,j]
         *
         * CORRECT FORMULA:
         *   frechetAdjExpBlock(A, X) = frechetExpBlock(A^T, X)
         *
         * NOT: (frechetExpBlock(A^T, X^T))^T
         *
         * This difference is CRITICAL and causes systematic errors if wrong.
         */
        int k = A.numRows;

        // Transpose A
        DenseMatrix64F A_T = new DenseMatrix64F(k, k);
        CommonOps.transpose(A, A_T);

        // Compute frechetExpBlock(A^T, X) directly
        // NO TRANSPOSE on X! NO TRANSPOSE on output!
        frechetExp(A_T, X, out);
    }

    public static void frechetExp(final DenseMatrix64F A,
                                       final DenseMatrix64F E,
                                       final DenseMatrix64F out) {
        final int d = A.numRows;

        // Build 2d×2d block matrix [[A, E], [0, A]]
        DenseMatrix64F B = new DenseMatrix64F(2 * d, 2 * d);
        CommonOps.insert(A, B, 0, 0);
        CommonOps.insert(E, B, 0, d);
        CommonOps.insert(A, B, d, d);

        // Compute exp(B) using Padé
        double[][] MM = new double[2 * d][2 * d];
        for (int i = 0; i < 2 * d; i++) {
            for (int j = 0; j < 2 * d; j++) {
                MM[i][j] = B.get(i, j);
            }
        }
        DoubleMatrix2D Mc = new DenseDoubleMatrix2D(MM);
        DoubleMatrix2D EB = MatrixExponential.expmPade13(Mc);

        // Extract upper-right block
        DoubleMatrix2D Kc = EB.viewPart(0, d, d, d).copy();
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                out.set(i, j, Kc.getQuick(i, j));
            }
        }
    }

}
