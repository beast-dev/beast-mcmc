package dr.util;

public class LKJCholeskyTransformConstrainedWithDiag extends LKJCholeskyTransformConstrained {

    public LKJCholeskyTransformConstrainedWithDiag(int dim) {
        super(dim);
    }

    @Override
    public double[] transform(double[] values, int from, int to) { //todo: to ensure Cholesky factor comes first in "values"

        assert from == 0 && to == values.length : "The transform function can only be applied to the whole array of values.";

        double[] choleskyFactor = new double[dim];
        double[] appendedTransformedValue = new double[2 * dim];

        System.arraycopy(values, 0, choleskyFactor, 0, dim);

        double[] CRCs = super.transform(choleskyFactor, 0, dim);

        System.arraycopy(CRCs, 0, appendedTransformedValue, 0, dim);
        System.arraycopy(values, dim, appendedTransformedValue, dim, dim);

        return appendedTransformedValue;
    }
}
