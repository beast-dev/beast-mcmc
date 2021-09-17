package dr.math.distributions;


import dr.inference.model.GradientProvider;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */

public class ConstrainedDeterminantDistributionModel implements MultivariateDistribution, GradientProvider {

    private final int dim;
    private final double shape;

    public ConstrainedDeterminantDistributionModel(double shape, int dim) {
        this.dim = dim;
        this.shape = shape;
    }


    @Override
    public int getDimension() {
        return dim * dim;
    }

    @Override
    public double[] getGradientLogDensity(Object x) {
        return gradLogPdf((double[]) x);
    }

    private double[] gradLogPdf(double[] x) {
        DenseMatrix64F X = DenseMatrix64F.wrap(dim, dim, x);

        DenseMatrix64F Xinv = new DenseMatrix64F(dim, dim);
        CommonOps.invert(X, Xinv);

        CommonOps.scale(shape, Xinv);
        CommonOps.transpose(Xinv);
        return Xinv.getData();
    }

    @Override
    public double logPdf(double[] x) {
        DenseMatrix64F X = DenseMatrix64F.wrap(dim, dim, x);
        double det = CommonOps.det(X);

        return shape * Math.log(Math.abs(det)); //TODO: normalizing constant
    }

    @Override
    public double[][] getScaleMatrix() {
        throw new RuntimeException("not implemented");
    }

    @Override
    public double[] getMean() {
        throw new RuntimeException("not implemented");
    }

    @Override
    public String getType() {
        return "ConstrainedDeterminant";
    }
}
