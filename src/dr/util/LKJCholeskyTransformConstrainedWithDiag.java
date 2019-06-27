package dr.util;


/**
 * @author Zhenyu Zhang
 */
public class LKJCholeskyTransformConstrainedWithDiag extends LKJCholeskyTransformConstrained {

    private int dimCPC;
    private final int totalDimension;

    LKJCholeskyTransformConstrainedWithDiag(int dim) {
        super(dim);
        this.dimCPC = dim * (dim - 1) / 2;
        this.totalDimension = dim * (dim + 1) / 2;
    }

    @Override
    public double[] transform(double[] values) {

        double[] choleskyFactor = subsetCholeskyOrCPCs(values);
        double[] diagonals = subsetDiagonals(values);
        double[] CPCs = super.transform(choleskyFactor, 0, dimCPC);

        return pasteTogether(CPCs, diagonals);
    }

    @Override
    //values: CPCs appended with log-transformed diagonals
    protected double[] inverse(double[] values) {

        assert values.length == totalDimension : "The transform function can only be applied to CPCs appended " +
                "with diagonals";
        for (int k = 0; k < dimCPC; k++) {
            assert values[k] <= 1.0 && values[k] >= -1.0 : "CPCs must be between -1.0 and 1.0";
        }

        double[] CPCs = subsetCholeskyOrCPCs(values);
        double[] diagonals = subsetDiagonals(values);
        double[] choleskyFactor = super.inverse(CPCs, 0, CPCs.length);
        return pasteTogether(choleskyFactor, diagonals);
    }

    @Override
    public double[][] computeJacobianMatrixInverse(double[] values) {

        double[] CPCs = subsetCholeskyOrCPCs(values);
        double[][] jacobian = super.computeJacobianMatrixInverse(CPCs);
        return appendIdentityMatrix(jacobian);
    }

    public double[] getGradientLogJacobianInverse(double[] values) {

        double[] CPCs = subsetCholeskyOrCPCs(values);
        double[] gradientLogJacobianInverse = super.getGradientLogJacobianInverse(CPCs);
        return pasteTogether(gradientLogJacobianInverse, new double[dimVector]);
    }

    private double[] subsetCholeskyOrCPCs(double[] values) { //todo: to ensure Cholesky factor comes first in "values"

        assert values.length == totalDimension;
        double[] choleskyOrCPC = new double[dimCPC];
        System.arraycopy(values, 0, choleskyOrCPC, 0, dimCPC);
        return choleskyOrCPC;
    }

    private double[] subsetDiagonals(double[] values) { //todo: to ensure Cholesky factor comes first in "values"

        assert values.length == totalDimension;
        double[] diagonals = new double[dimVector];
        System.arraycopy(values, dimCPC, diagonals, 0, dimVector);
        return diagonals;
    }

    private double[][] appendIdentityMatrix(double[][] jacobian) {

        assert jacobian.length == dimCPC;
        int length = dimCPC + dimVector;
        double[][] appendedJacobian = new double[length][length];

        for (int i = 0; i < length; i++) {
            for (int j = 0; j < length; j++) {
                if (i >= dimCPC || j >= dimCPC) {
                    if (i == j) {
                        appendedJacobian[i][j] = 1;
                    }
                } else {
                    appendedJacobian[i][j] = jacobian[i][j];
                }
            }
        }

        return appendedJacobian;
    }

    private double[] pasteTogether(double[] choleskyOrCPCs, double[] diagonals) {

        double[] concatenatedArray = new double[dimCPC + dimVector];
        System.arraycopy(choleskyOrCPCs, 0, concatenatedArray, 0, dimCPC);
        System.arraycopy(diagonals, 0, concatenatedArray, dimCPC, dimVector);
        return concatenatedArray;
    }

    public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value, int from, int to) {
        throw new RuntimeException("Not yet implemented");
    }
}
