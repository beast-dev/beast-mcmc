package dr.math.matrixAlgebra;

public class GramSchmidtProjection {

    public static void projectOrthonormal(WrappedMatrix matrix) {

        double[] sumSquares = projectOrthogonal(matrix);

        for (int j = 0; j < matrix.getMinorDim(); j++) {

            double norm = Math.sqrt(sumSquares[j]);

            for (int i = 0; i < matrix.getMajorDim(); i++) {
                matrix.set(i, j, matrix.get(i, j) / norm);
            }
        }
    }

    public static double[] projectOrthogonal(WrappedMatrix matrix) {
        int rowDim = matrix.getMajorDim();
        int colDim = matrix.getMinorDim();
        double[] sumSquares = new double[colDim];
        for (int col = 0; col < colDim; col++) {
            for (int col2 = 0; col2 < col; col2++) {

                double vdotu = 0;
                for (int row = 0; row < rowDim; row++) {
                    double xij = matrix.get(row, col);
                    double uij = matrix.get(row, col2);
                    vdotu += xij * uij;
                }

                vdotu = vdotu / sumSquares[col2];
                for (int row = 0; row < rowDim; row++) {
                    matrix.set(row, col, matrix.get(row, col) - vdotu * matrix.get(row, col2));
                }

            }

            double sumSquaresCol = 0;

            for (int row = 0; row < rowDim; row++) {
                double x = matrix.get(row, col);
                sumSquaresCol += x * x;
            }

            sumSquares[col] = sumSquaresCol;
        }
        return sumSquares;
    }


}
