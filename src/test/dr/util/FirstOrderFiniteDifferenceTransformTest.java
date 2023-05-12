package test.dr.util;

import dr.inference.model.Parameter;
import dr.math.MachineAccuracy;
import dr.math.MathUtils;
import dr.util.FirstOrderFiniteDifferenceTransform;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

public class FirstOrderFiniteDifferenceTransformTest extends TestCase {
    private final int d = 6;
    private final double acceptableProportionateError = 1e-8;
    private double[] unconstrained;
    private double[] sumScale;
    private double[] exponentialScale;
    private double[] logisticScale;


    private final FirstOrderFiniteDifferenceTransform noneFOFDT = new FirstOrderFiniteDifferenceTransform(d, FirstOrderFiniteDifferenceTransform.IncrementTransform.NONE, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
    private final FirstOrderFiniteDifferenceTransform logFOFDT = new FirstOrderFiniteDifferenceTransform(d, FirstOrderFiniteDifferenceTransform.IncrementTransform.LOG, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
    private final FirstOrderFiniteDifferenceTransform logitFOFDT = new FirstOrderFiniteDifferenceTransform(d, FirstOrderFiniteDifferenceTransform.IncrementTransform.LOGIT, 0.0, 3.14159);

    public FirstOrderFiniteDifferenceTransformTest(String name) { super(name); }

    public void setUp() throws Exception {
        super.setUp();

        unconstrained = new double[d];
        sumScale = new double[d];
        exponentialScale = new double[d];
        logisticScale = new double[d];

        // Values that don't get too far from 0 so as not to hit Infinities with the log-scale version
        for (int i = 0; i < d; i++) {
            unconstrained[i] = (i - 2.0)/12.0;
        }
        sumScale = noneFOFDT.inverse(unconstrained, 0, d);
        exponentialScale = logFOFDT.inverse(unconstrained, 0, d);
//        logisticScale = logitFOFDT.inverse(unconstrained, 0, d);

    }

    private double getH(double x) {
        // As in NumericalDerivative.java
        return MachineAccuracy.SQRT_EPSILON*(Math.abs(x) + 1.0);
    }

    private double getRelativeTolerance(double x) {
        return Math.abs(x) * acceptableProportionateError;
    }

    // jacobian[j][i] = d x_i / d y_j
    private double[][] computeNumericalJacobianInverse(double[] values, FirstOrderFiniteDifferenceTransform transform) {
        int dim = transform.getDimension();
        double[] y = values.clone();
        double[] x = transform.inverse(values, 0, dim);

        double[][] jacobian = new double[dim][dim];

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                double h = getH(values[i]);
                y[j] = values[j] + h;
                x = transform.inverse(y, 0, dim);
                double tmp = x[i];
                y[j] = values[j] - h;
                x = transform.inverse(y, 0, dim);
                tmp -= x[i];
                y[j] = values[j];
                jacobian[j][i] = tmp/(2.0 * h);
            }
        }

        return jacobian;
    }

    public double[] getNumericalGradientLogJacobianInverse(double[] values, FirstOrderFiniteDifferenceTransform transform) {
        int dim = transform.getDimension();
        double[] tmpv = values.clone();
        double[] numGrad = new double[dim];
        for (int i = 0; i < dim; i++) {
            double h = getH(values[i]);
            tmpv[i] = values[i] + h;
            double fx_plus_h = transform.getLogJacobianInverse(tmpv);
            tmpv[i] = values[i] - h;
            double fx_minus_h = transform.getLogJacobianInverse(tmpv);
            tmpv[i] = values[i];
            numGrad[i] = (fx_plus_h - fx_minus_h)/(2 * h);
        }
//        tmpv = transform.inverse(tmpv, 0, dim);
        return numGrad;
    }

    public void testForward() {
        double[] transformed = noneFOFDT.transform(sumScale, 0, d);
        for (int i = 0; i < d; i++) {
            assertEquals(unconstrained[i], transformed[i], 1e-8);
        }

        transformed = logFOFDT.transform(exponentialScale, 0, d);
        for (int i = 0; i < d; i++) {
            assertEquals(unconstrained[i], transformed[i], 1e-8);
        }
    }

    public void testReverse() {
        assertEquals(-0.1666666666666667, sumScale[0], 1e-8);
        assertEquals(-0.2500000000000000, sumScale[1], 1e-8);
        assertEquals(-0.2500000000000000, sumScale[2], 1e-8);
        assertEquals(-0.1666666666666667, sumScale[3], 1e-8);
        assertEquals( 0.0000000000000000, sumScale[4], 1e-8);
        assertEquals( 0.2500000000000000, sumScale[5], 1e-8);

        assertEquals(0.8464817248906141, exponentialScale[0], 1e-8);
        assertEquals(0.7788007830714049, exponentialScale[1], 1e-8);
        assertEquals(0.7788007830714049, exponentialScale[2], 1e-8);
        assertEquals(0.8464817248906140, exponentialScale[3], 1e-8);
        assertEquals(1.0000000000000000, exponentialScale[4], 1e-8);
        assertEquals(1.2840254166877414, exponentialScale[5], 1e-8);
    }

    public void testJacobian() {
        double[][] jacobianInverse = new double[d][d];
        double[][] numericJacobianInverse = new double[d][d];

        jacobianInverse = noneFOFDT.computeJacobianMatrixInverse(sumScale);
        numericJacobianInverse = computeNumericalJacobianInverse(sumScale, noneFOFDT);
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                assertEquals(numericJacobianInverse[i][j], jacobianInverse[i][j], getRelativeTolerance(jacobianInverse[i][j]));
            }
        }

        jacobianInverse = logFOFDT.computeJacobianMatrixInverse(exponentialScale);
        numericJacobianInverse = computeNumericalJacobianInverse(exponentialScale, logFOFDT);
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                assertEquals(numericJacobianInverse[i][j], jacobianInverse[i][j], getRelativeTolerance(jacobianInverse[i][j]));
            }
        }

    }

    public void testGradientLogJacobian() {
        double[] gradLogJacobianInv = new double[d];
        double[] numGradLogJacobianInv = new double[d];

        gradLogJacobianInv = noneFOFDT.getGradientLogJacobianInverse(sumScale);
        numGradLogJacobianInv = getNumericalGradientLogJacobianInverse(sumScale, noneFOFDT);
        for (int i = 0; i < d; i++) {
            assertEquals(numGradLogJacobianInv[i], gradLogJacobianInv[i], getRelativeTolerance(gradLogJacobianInv[i]));
        }

        gradLogJacobianInv = logFOFDT.getGradientLogJacobianInverse(exponentialScale);
        numGradLogJacobianInv = getNumericalGradientLogJacobianInverse(exponentialScale, logFOFDT);
        for (int i = 0; i < d; i++) {
            assertEquals(numGradLogJacobianInv[i], gradLogJacobianInv[i], getRelativeTolerance(gradLogJacobianInv[i]));
        }

    }

    public void testFirstOrderFiniteDifferenceTransform() {


//        System.err.println("unconstrained (\"transformed\") values: " + new dr.math.matrixAlgebra.Vector(unconstrained));

        // TODO: add test of logit-scale to all these components!
        // TODO: add test of updateGradientLogDensity once O(n) implementation is complete
        testForward();
        testReverse();
        testJacobian();
        testGradientLogJacobian();

    }

    public static Test suite() {
        return new TestSuite(FirstOrderFiniteDifferenceTransformTest.class);
    }
}
