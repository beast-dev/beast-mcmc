package test.dr.inference;

import dr.inference.model.CorrelationSymmetricMatrix;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.Matrix;
import test.dr.math.MathTestCase;

public class CorrelationSymmetricMatrixTest extends MathTestCase {

    public void testGetAttributeValue() {

        CorrelationSymmetricMatrix matrix = getMatrix(CorrelationSymmetricMatrix.Type.AS_CORRELATION);
        double[] values = matrix.getAttributeValue();

        assertEquals(values.length, matrix.getRowDimension() * matrix.getColumnDimension());
    }

    public void testGetParameterValue() {

        CorrelationSymmetricMatrix matrix = getMatrix(CorrelationSymmetricMatrix.Type.AS_CORRELATION);
        double value = matrix.getParameterValue(0, 0);
        assertEquals(value, 1.0);

        double top = matrix.getParameterValue(2, 1);
        double bottom = matrix.getParameterValue(1, 2);
        assertEquals(top, bottom);

        double test = 4.0 * Math.sqrt(2 * 3);
        assertEquals(top, test);
    }

    public void testGetParameterAsMatrix() {

        CorrelationSymmetricMatrix matrix = getMatrix(CorrelationSymmetricMatrix.Type.AS_IS);
        double[][] M = matrix.getParameterAsMatrix();
        System.out.println(new Matrix(M));

        double sum = 0.0;
        for (double x : M[0]) {
            sum += x;
        }

        assertEquals(sum, 7.0);

    }

    private CorrelationSymmetricMatrix getMatrix(CorrelationSymmetricMatrix.Type type) {
            Parameter diagonals = new Parameter.Default(new double[] { 1.0, 2.0, 3.0, 4.0 });
            Parameter offDiagonals = new Parameter.Default(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0});

            return new CorrelationSymmetricMatrix(diagonals, offDiagonals, type);
    }
}