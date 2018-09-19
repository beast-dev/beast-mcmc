package dr.util;

public class LKJCholeskyTransformConstrainedWithDiag extends LKJCholeskyTransformConstrained {

    protected int CPCdimension;

    public LKJCholeskyTransformConstrainedWithDiag(int dim) {
        super(dim);
        this.CPCdimension = dim * (dim - 1) / 2;
    }

    @Override
    public double[] transform(double[] values, int from, int to) {

        assert from == 0 && to == values.length : "The transform function can only be applied to the whole array of " +
                "values.";

        double[] choleskyFactor = subsetCholeskyOrCPCs(values);
        double[] diagonals = subsetDiagonals(values);
        double[] CPCs = super.transform(choleskyFactor, 0, dim);

        return pasteTogether(CPCs, diagonals);
    }

    @Override
    //values: CPCs appended with log-transformed diagonals
    public double[] inverse(double[] values, int from, int to) {

        assert from == 0 && to == values.length : "The transform function can only be applied to the whole array of " +
                "values.";
        assert values.length == dim * (dim + 1) / 2 : "The transform function can only be applied to CPCs appended " +
                "with diagonals";
        for (int k = 0; k < dim * (dim - 1) / 2; k++) {
            assert values[k] <= 1.0 && values[k] >= -1.0 : "CPCs must be between -1.0 and 1.0";
        }

        double[] CPCs = subsetCholeskyOrCPCs(values);
        double[] diagonals = subsetDiagonals(values);

        return pasteTogether(CPCs, diagonals);
    }

    @Override
    public double[][] computeJacobianMatrixInverse(double[] values) {

        double[] choleskyFactor = subsetCholeskyOrCPCs(values);

        double[][] jacobian = super.computeJacobianMatrixInverse(choleskyFactor);

        return appendZeros(jacobian);
    }

    private double[] subsetCholeskyOrCPCs(double[] values) { //todo: to ensure Cholesky factor comes first in "values"

        assert values.length == dim * (dim + 1) / 2;
        double[] choleskyFactor = new double[CPCdimension];
        System.arraycopy(values, 0, choleskyFactor, 0, CPCdimension);
        return choleskyFactor;
    }

    private double[] subsetDiagonals(double[] values) { //todo: to ensure Cholesky factor comes first in "values"

        assert values.length == dim * (dim + 1) / 2;
        double[] diagonals = new double[dim];
        System.arraycopy(values, CPCdimension, diagonals, 0, dim);
        return diagonals;
    }


    private double[][] appendZeros(double[][] jacobian) {

        assert jacobian.length == 3;
        int length = dim * (dim - 1) / 2 + dim;
        double[][] appendedJacobian = new double[length][length];

        for (int i = 0; i < length; i++) {
            for (int j = 0; j < length; j++) {
                if (i >= dim || j >= dim) {
                    appendedJacobian[i][j] = 0;
                } else {
                    appendedJacobian[i][j] = jacobian[i][j];
                }
            }
        }

        return appendedJacobian;
    }

    private double[] pasteTogether(double[] choleskyOrCPCs, double[] diagonals) {

        double[] concatenatedArray = new double[CPCdimension + dim];
        System.arraycopy(choleskyOrCPCs, 0, concatenatedArray, 0, CPCdimension);
        System.arraycopy(diagonals, 0, concatenatedArray, CPCdimension, dim);
        return concatenatedArray;
    }
}
