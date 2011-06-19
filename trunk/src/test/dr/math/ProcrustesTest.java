package test.dr.math;

import dr.math.Procrustes;
import dr.math.matrixAlgebra.Matrix;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;

/* R test code

library(MCMCpack)
set.seed(666)
Xstar = matrix(runif(9,-10,10), nrow=3)
X = matrix(runif(9,-2,2), nrow=3)
out.0 = procrustes(X, Xstar, F, F)
out.1 = procrustes(X, Xstar, T, F)
out.2 = procrustes(X, Xstar, F, T)

*/

/**
 * @author Marc A. Suchard
 * @author Andrew Rambaut
 */
public class ProcrustesTest extends MathTestCase {

    public void testProcrustes() {
        double tolerance = 1E-4;
        runTest(X, Xstar, gldStd[0], false, false, tolerance, "No translation, no dilation");
        runTest(X, Xstar, gldStd[1], true, false, tolerance,  "Yes translation, no dilation");
        runTest(X, Xstar, gldStd[2], false, true, tolerance,  "No translation, yes dilation");

    }

    private void runTest(double[][] X, double[][] Xstar, double[][] goldStandard,
                         boolean allowTranslation, boolean allowDilation, double tolerance, String header) {

        RealMatrix rmX = new Array2DRowRealMatrix(X);
        RealMatrix rmXstar = new Array2DRowRealMatrix(Xstar);
        Procrustes procrustes = new Procrustes(rmX, rmXstar, allowTranslation, allowDilation);
        RealMatrix rmXnew = procrustes.procrustinate(rmX);

        System.out.println(header);
        System.out.print("Translation:\n" + new Matrix(procrustes.getTranslation().getData()));
        System.out.println("Dilation = " + procrustes.getDilation());
        System.out.println("Xnew:");
        System.out.println(new Matrix(rmXnew.getData()));
        assertEquals(rmXnew.getData(), goldStandard, tolerance);
    }

    private static double[][] X = {
            { -0.9602155, -1.617021,  1.245026 },
            {  1.1035723, -1.431346, -1.853811 },
            { -1.9344838, -1.155495,  1.566550 }
    };

    
    private static double[][] Xstar = {
            {  5.487370, -5.973453,  9.57456879 },
            { -6.055516, -2.775111, -0.03772582 },
            {  9.560277,  4.852239, -9.73368329 }
    };

    private static double[][][] gldStd = {
            {
                    {  2.062008, -0.8829124,  0.2354720 },  // No translation, no dilation
                    { -1.456369, -2.1056073,  0.3855449 },
                    {  2.669785, -0.3375384, -0.5383219 },
            },
            {
                    { 3.9949788, -1.2385190,  0.1258216 },
                    { 0.3349586, -1.9233277,  0.3040794 },
                    { 4.6621931, -0.7344788, -0.6267413 }
            },
            {

                    {  6.682553, -2.861341,  0.7631174 },
                    { -4.719800, -6.823849,  1.2494734 },
                    {  8.652235, -1.093894, -1.7445929 },
            }
    };
}
