package test.dr.util;

import dr.inference.model.Parameter;
import dr.math.MachineAccuracy;
import dr.math.MathUtils;
import dr.util.FirstOrderFiniteDifferenceTransform;
import dr.util.Transform;
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


    private final FirstOrderFiniteDifferenceTransform noneFOFDT = new FirstOrderFiniteDifferenceTransform(d, new Transform.NoTransform());
    private final FirstOrderFiniteDifferenceTransform logFOFDT = new FirstOrderFiniteDifferenceTransform(d, new Transform.LogTransform());
    private final FirstOrderFiniteDifferenceTransform logitFOFDT = new FirstOrderFiniteDifferenceTransform(d, new Transform.ScaledLogitTransform(Math.PI, Math.E));

    public FirstOrderFiniteDifferenceTransformTest(String name) { super(name); }

    public void setUp() throws Exception {
        super.setUp();

        // Set up values on all scales
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

    public static double getH(double x) {
        // As in NumericalDerivative.java
        return MachineAccuracy.SQRT_EPSILON*(Math.abs(x) + 1.0);
    }

    public static double getTolerance(double x, double acceptableProportionateError) {
        return Math.abs(x) * acceptableProportionateError;
    }

    // jacobian[j][i] = d x_i / d y_j
    public static double[][] computeNumericalJacobianInverse(double[] values, Transform.MultivariateTransform transform) {
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

    public static double getLogJacobianInverse(double[] values, Transform.MultivariateTransform transform, int dim) {
        return -transform.getLogJacobian(transform.inverse(values, 0, dim), 0, dim);
    }

    public static double[] getNumericalGradientLogJacobianInverse(double[] values, Transform.MultivariateTransform transform) {
        int dim = transform.getDimension();
        double[] tmpv = values.clone();
        double[] numGrad = new double[dim];
        for (int i = 0; i < dim; i++) {
            double h = getH(values[i]);
            tmpv[i] = values[i] + h;
            double fx_plus_h = getLogJacobianInverse(tmpv, transform, dim);
            tmpv[i] = values[i] - h;
            double fx_minus_h = getLogJacobianInverse(tmpv, transform, dim);
            tmpv[i] = values[i];
            numGrad[i] = (fx_plus_h - fx_minus_h)/(2 * h);
        }
//        tmpv = transform.inverse(tmpv, 0, dim);
        return numGrad;
    }

    public static double[] updateGradLogDens(double[] gradient, double[] value, double[] gradientLogJacobianInverse, Transform.MultivariateTransform transform) {
        // values = untransformed (R)
        double[] transformedValues = transform.transform(value, 0, value.length);
        double[][] jacobianInverse = transform.computeJacobianMatrixInverse(transformedValues);
        double[] updatedGradient = new double[gradient.length];
        for (int i = 0; i < gradient.length; i++) {
            for (int j = 0; j < gradient.length; j++) {
                updatedGradient[i] += jacobianInverse[i][j] * gradient[j];
            }
        }
//        double[] gradientLogJacobianInverse = transform.getGradientLogJacobianInverse(transformedValues);
        // Add gradient log jacobian
        for (int i = 0; i < gradient.length; i++) {
            updatedGradient[i] += gradientLogJacobianInverse[i];
        }
        return updatedGradient;
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

    // computeJacobianMatrixInverse is called on the transformed values
    public void testJacobian() {
        double[][] jacobianInverse = new double[d][d];
        double[][] numericJacobianInverse = new double[d][d];

        jacobianInverse = noneFOFDT.computeJacobianMatrixInverse(unconstrained);
        numericJacobianInverse = computeNumericalJacobianInverse(unconstrained, noneFOFDT);
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                assertEquals(numericJacobianInverse[i][j], jacobianInverse[i][j], getTolerance(jacobianInverse[i][j], acceptableProportionateError));
            }
        }

        jacobianInverse = logFOFDT.computeJacobianMatrixInverse(unconstrained);
        numericJacobianInverse = computeNumericalJacobianInverse(unconstrained, logFOFDT);
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                assertEquals(numericJacobianInverse[i][j], jacobianInverse[i][j], getTolerance(jacobianInverse[i][j], acceptableProportionateError));
            }
        }

    }

    // getGradientLogJacobianInverse is called on the transformed values
    public void testGradientLogJacobian() {
        double[] gradLogJacobianInv = new double[d];
        double[] numGradLogJacobianInv = new double[d];

        gradLogJacobianInv = noneFOFDT.getGradientLogJacobianInverse(unconstrained);
        numGradLogJacobianInv = getNumericalGradientLogJacobianInverse(unconstrained, noneFOFDT);
        for (int i = 0; i < d; i++) {
            assertEquals(numGradLogJacobianInv[i], gradLogJacobianInv[i], getTolerance(gradLogJacobianInv[i], acceptableProportionateError));
        }

        gradLogJacobianInv = logFOFDT.getGradientLogJacobianInverse(unconstrained);
        numGradLogJacobianInv = getNumericalGradientLogJacobianInverse(unconstrained, logFOFDT);
        for (int i = 0; i < d; i++) {
            assertEquals(numGradLogJacobianInv[i], gradLogJacobianInv[i], getTolerance(gradLogJacobianInv[i], acceptableProportionateError));
        }
    }

    private void test1UpdateGradient(double[] gradient, double[] untransformed, FirstOrderFiniteDifferenceTransform transform) {
        double[] transformed = transform.transform(untransformed, 0, d);
        double[] gradLogJacobian = transform.getGradientLogJacobianInverse(transformed);
        double[] updatedFullMatrix = updateGradLogDens(gradient, untransformed, gradLogJacobian, transform);
        double[] updated = transform.updateGradientLogDensity(gradient, untransformed, 0, d);
        for (int i = 0; i < d; i++) {
            assertEquals(updatedFullMatrix[i], updated[i], getTolerance(updated[i], acceptableProportionateError));
        }
    }

    // updateGradientLogDensity is called on the untransformed scale
    public void testUpdateGradient() {
        // An arbitrary value to take as the gradient
        double[] gradient = new double[d];
        for (int i = 0; i < d; i++) {
            gradient[i] = (((double)i) - 3.0) * 42.0;
        }

        test1UpdateGradient(gradient, unconstrained, noneFOFDT);
        test1UpdateGradient(gradient, unconstrained, logFOFDT);
    }

    public void testFirstOrderFiniteDifferenceTransform() {
        // TODO: add test of logit-scale to all these components!
        testForward();
        testReverse();
        testJacobian();
        testGradientLogJacobian();
        testUpdateGradient();
    }

    public static Test suite() {
        return new TestSuite(FirstOrderFiniteDifferenceTransformTest.class);
    }
}
