package test.dr.math;

import junit.framework.TestCase;
import dr.math.matrixAlgebra.Vector;
import dr.math.KroneckerOperation;

/**
 * @author Marc A. Suchard
 */
public class KroneckerOperationTest extends TestCase {

    private void testKroneckerSum(double[] A, int m, double[] B, int n, double[] trueC) {

        double[] C = KroneckerOperation.sum(A, m, B, n);
        printMatrix(A, m, m, "A = ");
        printMatrix(B, n, n, "B = ");
        printMatrix(C, m * n, m * n, "A %s% B = ");
        didPass(C, trueC);
    }

    private void testKroneckerProduct(double[] A, int m, int n, double[] B, int p, int q, double[] trueC) {

        double[] C = KroneckerOperation.product(A, m, n, B, p, q);
        printMatrix(A, m, n, "A = ");
        printMatrix(B, p, q, "B = ");
        printMatrix(C, m * p, n * q, "A %x% B = ");
        didPass(C, trueC);
    }

    private void printMatrix(double[] A, int m, int n, String text) {
        System.out.println(text + new Vector(A) + " (" + m + " x " + n + ")");
    }

    private void didPass(double[] C, double[] trueC) {
        boolean passed = true;
        for (int i = 0; i < C.length; i++) {
            assertEquals(C[i], trueC[i], 1E-10);
            if (C[i] != trueC[i]) {
                passed = false;
                break;
            }
        }
        System.out.println("Passed = " + passed + "\n");
    }

    public void testKroneckerSum() {
        double[] A = {
                1, 2,
                3, 4};
        double[] B = {
                0, 5,
                6, 7};
        double[] trueC = {
                1, 5, 2, 0,
                6, 8, 0, 2,
                3, 0, 4, 5,
                0, 3, 6, 11};
        testKroneckerSum(A, 2, B, 2, trueC);

        double[] D = {
                -4, 2,
                5, 6
        };
        double[] E = {
                1, 2, 3,
                4, 5, 6,
                7, 8, 9
        };
        double[] trueF = {
                -3, 2, 3, 2, 0, 0,
                4, 1, 6, 0, 2, 0,
                7, 8, 5, 0, 0, 2,
                5, 0, 0, 7, 2, 3,
                0, 5, 0, 4, 11, 6,
                0, 0, 5, 7, 8, 15
        };
        testKroneckerSum(D, 2, E, 3, trueF);
    }

    public void testKroneckerProduct() {

        double[] A = {1, 2, 3, 4};
        double[] B = {0, 5, 6, 7};
        double[] trueC = {0, 5, 0, 10, 6, 7, 12, 14, 0, 15, 0, 20, 18, 21, 24, 28};
        testKroneckerProduct(A, 2, 2, B, 2, 2, trueC);

        double[] D = {
                1, 2, 3,
                3, 2, 1};
        double[] E = {
                2, 1,
                2, 3};
        double[] trueF = {
                2, 1, 4, 2, 6, 3,
                2, 3, 4, 6, 6, 9,
                6, 3, 4, 2, 2, 1,
                6, 9, 4, 6, 2, 3};
        testKroneckerProduct(D, 2, 3, E, 2, 2, trueF);
    }
}
