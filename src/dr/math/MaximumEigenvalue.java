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

    class PowerMethod implements  MaximumEigenvalue {

        private final int numIterations;
        private final double err;

        public PowerMethod(int numIterations, double err) {
            this.numIterations = numIterations;
            this.err = err;
        }

        @Override
        public double find(double[][] matrix) {

            final ReadableMatrix mat = new WrappedMatrix.ArrayOfArray(matrix);

            WrappedVector y = getInitialGuess(matrix.length);
            double maxEigenvalue = Double.NEGATIVE_INFINITY;

            for (int i = 0; i < numIterations; ++i) {

                ReadableVector v = new ReadableVector.Scale(1 / norm(y), y);
                y = product(mat, v);
                maxEigenvalue = innerProduct(v, y);

                ReadableVector diff = new ReadableVector.Sum(y,
                        new ReadableVector.Scale(-maxEigenvalue, v));

                if (ReadableVector.Utils.norm(diff) < err) {
                    break;
                }
            }

            return maxEigenvalue;
        }

        private static WrappedVector getInitialGuess(int dim) {

            double[] y = new double[dim];
            for (int i = 0; i < dim; ++i) {
                y[i] = MathUtils.nextDouble();
            }

            return new WrappedVector.Raw(y);
        }
    }
}
