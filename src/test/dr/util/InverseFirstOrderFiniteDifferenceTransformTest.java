package test.dr.util;

import dr.inference.model.Parameter;
import dr.math.MachineAccuracy;
import dr.math.MathUtils;
import dr.util.FirstOrderFiniteDifferenceTransform;
import dr.util.InverseFirstOrderFiniteDifferenceTransform;
import dr.util.Transform;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import test.dr.util.FirstOrderFiniteDifferenceTransformTest;

public class InverseFirstOrderFiniteDifferenceTransformTest extends TestCase {
    private final int d = 6;
    private final double acceptableProportionateError = 1e-8;
    private double[] unconstrained;
    private double[] sumScale;
    private double[] exponentialScale;
    private double[] logisticScale;


    private final InverseFirstOrderFiniteDifferenceTransform noneIFOFDT = new InverseFirstOrderFiniteDifferenceTransform(d, new Transform.NoTransform());
    private final InverseFirstOrderFiniteDifferenceTransform logIFOFDT = new InverseFirstOrderFiniteDifferenceTransform(d, new Transform.LogTransform());
    private final InverseFirstOrderFiniteDifferenceTransform logitIFOFDT = new InverseFirstOrderFiniteDifferenceTransform(d, new Transform.ScaledLogitTransform(Math.PI, Math.E));

    public InverseFirstOrderFiniteDifferenceTransformTest(String name) { super(name); }

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
        sumScale = noneIFOFDT.transform(unconstrained, 0, d);
        exponentialScale = logIFOFDT.transform(unconstrained, 0, d);
//        logisticScale = logitIFOFDT.transform(unconstrained, 0, d);

    }

    public void testForward() {
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

    public void testReverse() {
        double[] inverse = noneIFOFDT.inverse(sumScale, 0, d);
        for (int i = 0; i < d; i++) {
            assertEquals(unconstrained[i], inverse[i], 1e-8);
        }

        inverse = logIFOFDT.inverse(exponentialScale, 0, d);
        for (int i = 0; i < d; i++) {
            assertEquals(unconstrained[i], inverse[i], 1e-8);
        }
    }

    // computeJacobianMatrixInverse is called on the transformed values
    public void testJacobian() {
        double[][] jacobianInverse = new double[d][d];
        double[][] numericJacobianInverse = new double[d][d];

        jacobianInverse = noneIFOFDT.computeJacobianMatrixInverse(sumScale);
        numericJacobianInverse = FirstOrderFiniteDifferenceTransformTest.computeNumericalJacobianInverse(sumScale, noneIFOFDT);
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                assertEquals(numericJacobianInverse[i][j], jacobianInverse[i][j], FirstOrderFiniteDifferenceTransformTest.getTolerance(jacobianInverse[i][j], acceptableProportionateError));
            }
        }

        jacobianInverse = logIFOFDT.computeJacobianMatrixInverse(exponentialScale);
        numericJacobianInverse = FirstOrderFiniteDifferenceTransformTest.computeNumericalJacobianInverse(exponentialScale, logIFOFDT);

        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                assertEquals(numericJacobianInverse[i][j], jacobianInverse[i][j], FirstOrderFiniteDifferenceTransformTest.getTolerance(jacobianInverse[i][j], acceptableProportionateError));
            }
        }

    }

    // getGradientLogJacobianInverse is called on the transformed values
    public void testGradientLogJacobian() {
        double[] gradLogJacobianInv = new double[d];
        double[] numGradLogJacobianInv = new double[d];

        gradLogJacobianInv = noneIFOFDT.getGradientLogJacobianInverse(sumScale);
        numGradLogJacobianInv = FirstOrderFiniteDifferenceTransformTest.getNumericalGradientLogJacobianInverse(sumScale, noneIFOFDT);
        for (int i = 0; i < d; i++) {
            assertEquals(numGradLogJacobianInv[i], gradLogJacobianInv[i], FirstOrderFiniteDifferenceTransformTest.getTolerance(gradLogJacobianInv[i], acceptableProportionateError));
        }

        gradLogJacobianInv = logIFOFDT.getGradientLogJacobianInverse(exponentialScale);
        numGradLogJacobianInv = FirstOrderFiniteDifferenceTransformTest.getNumericalGradientLogJacobianInverse(exponentialScale, logIFOFDT);
        for (int i = 0; i < d; i++) {
            assertEquals(numGradLogJacobianInv[i], gradLogJacobianInv[i], FirstOrderFiniteDifferenceTransformTest.getTolerance(gradLogJacobianInv[i], acceptableProportionateError));
        }

    }

    private void test1UpdateGradient(double[] gradient, double[] untransformed, InverseFirstOrderFiniteDifferenceTransform transform) {
        System.err.println(transform);
        double[] transformed = transform.transform(untransformed, 0, d);
        double[] gradLogJacobian = transform.getGradientLogJacobianInverse(transformed);
        double[] updatedFullMatrix = FirstOrderFiniteDifferenceTransformTest.updateGradLogDens(gradient, untransformed, gradLogJacobian, transform);
        double[] updated = transform.updateGradientLogDensity(gradient, untransformed, 0, d);
        for (int i = 0; i < d; i++) {
            assertEquals(updatedFullMatrix[i], updated[i], FirstOrderFiniteDifferenceTransformTest.getTolerance(updated[i], acceptableProportionateError));
        }
    }

    // updateGradientLogDensity is called on the untransformed scale
    public void testUpdateGradient() {
        // An arbitrary value to take as the gradient
        double[] gradient = new double[d];
        for (int i = 0; i < d; i++) {
            gradient[i] = (((double)i) - 3.0) * 42.0;
        }

        test1UpdateGradient(gradient, unconstrained, noneIFOFDT);
        test1UpdateGradient(gradient, unconstrained, logIFOFDT);
    }


//    public void testInverseFirstOrderFiniteDifferenceTransform() {
//        // TODO: add test of logit-scale to all these components!
//        testForward();
//        testReverse();
//        testJacobian();
//        testGradientLogJacobian();
////        testUpdateGradient();
//    }

    public static Test suite() {
        return new TestSuite(InverseFirstOrderFiniteDifferenceTransformTest.class);
    }
}
