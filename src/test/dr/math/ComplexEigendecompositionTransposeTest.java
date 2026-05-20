package test.dr.math;

import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.EigenSystem;
import dr.evomodel.substmodel.LapackEigenSystem;

import static dr.evomodel.substmodel.EigenDecomposition.*;

public class ComplexEigendecompositionTransposeTest extends MathTestCase {

    public void testTranspose() {

        double[][] matrix = {
                {-2, 0, 2},
                {1, -1, 0},
                {0, 1, -1}
        };

        double[][] matrixT = {
                {-2, 1, 0},
                {0, -1, 1},
                {2, 0, -1}
        };

        EigenSystem es = new LapackEigenSystem(3);

        EigenDecomposition edGold = es.decomposeMatrix(matrixT);
        double[] resultGold = new double[9];
        es.computeExponential(edGold, 1.0, resultGold);

        EigenDecomposition ed = es.decomposeMatrix(matrix);
        double[] resultRescale = new double[9];
        double[] resultConjugate = new double[9];

        EigenDecomposition tedRescale = ed.transpose();
        es.computeExponential(tedRescale, 1.0, resultRescale);

        EigenDecomposition tedConjugate = transposeViaConjugate(ed);
        es.computeExponential(tedConjugate, 1.0, resultConjugate);

        assertEquals(resultGold, resultRescale, 1E-12);
        assertEquals(resultGold, resultConjugate, 1E-12);
    }

    private EigenDecomposition transposeViaConjugate(EigenDecomposition in) {
        // note: exchange e/ivec
        int dim = (int) Math.sqrt(in.getInverseEigenVectors().length);
        double[] evec = in.getInverseEigenVectors().clone();
        transposeInPlace(evec, dim);

        double[] ievc = in.getEigenVectors().clone();
        transposeInPlace(ievc, dim);

        double[] eval = in.getEigenValues().clone();
        conjugate(eval, dim);
        return new EigenDecomposition(evec, ievc, eval);
    }

    public static void conjugate(double[] eval, int dim) {
        if (eval.length != 2 * dim) {
            return;
        }
        for (int i = 0; i < dim; ++i) {
            if (eval[dim + i] != 0.0) {
                eval[dim + i] = -eval[dim + i];
            }
        }
    }
}
