package dr.math;

import dr.math.matrixAlgebra.*;

import static dr.math.matrixAlgebra.ReadableMatrix.Utils.product;
import static dr.math.matrixAlgebra.ReadableVector.Utils.innerProduct;
import static dr.math.matrixAlgebra.ReadableVector.Utils.norm;

/**
 * @author Marc A. Suchard
 */
public interface MaximumEigenvalue {

    double find(double[][] matrix);

    public class PowerMethod implements  MaximumEigenvalue {

        private final int numIterations;
        private final double err;

        public PowerMethod(int numIterations, double err) {
            this.numIterations = numIterations;
            this.err = err;
        }

        @Override
        public double find(double[][] matrix) {

            double[] y0 = new double[matrix.length];
            ReadableVector diff;
            double maxEigenvalue = 10.0; // TODO Bad magic number

            for (int i = 0; i < matrix.length; ++i) {
                y0[i] = MathUtils.nextDouble();
            }

            WrappedVector y = new WrappedVector.Raw(y0);
            final ReadableMatrix mat = new WrappedMatrix.ArrayOfArray(matrix);

            for (int i = 0; i < numIterations; ++i) {

                ReadableVector v = new ReadableVector.Scale(1 / norm(y), y);
                y = product(mat, v);
                maxEigenvalue = innerProduct(v, y);
                diff = new ReadableVector.Sum(y,
                        new ReadableVector.Scale(-maxEigenvalue, v));

                if (ReadableVector.Utils.norm(diff) < err) {
                    break;
                }
            }

            return maxEigenvalue;
        }
    }
}
