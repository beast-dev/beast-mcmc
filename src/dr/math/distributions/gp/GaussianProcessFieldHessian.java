package dr.math.distributions.gp;

import dr.inference.distribution.RandomField;
import dr.inference.distribution.RandomFieldHessian;
import dr.inference.model.Parameter;

/**
 * Hessian provider for Gaussian process random fields.
 */
public class GaussianProcessFieldHessian extends RandomFieldHessian {

    public GaussianProcessFieldHessian(RandomField randomField,
                                       Parameter parameter) {
        super(randomField, parameter);
        GaussianProcessFieldGradient.getGaussianProcessField(randomField);
    }

    static double[] getDiagonalHessianLogDensity(GaussianProcessField distribution) {
        final int dim = distribution.getDimension();
        final double[] precision = distribution.getPrecision();
        final double[] diagonal = new double[dim];
        for (int i = 0; i < dim; ++i) {
            diagonal[i] = -precision[i * dim + i];
        }
        return diagonal;
    }

    static double[][] getHessianLogDensity(GaussianProcessField distribution) {
        final int dim = distribution.getDimension();
        final double[] precision = distribution.getPrecision();
        final double[][] hessian = new double[dim][dim];
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                hessian[i][j] = -precision[i * dim + j];
            }
        }
        return hessian;
    }
}
