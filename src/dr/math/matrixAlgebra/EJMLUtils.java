package dr.math.matrixAlgebra;

import org.ejml.data.DenseMatrix64F;

import java.util.Arrays;

public class EJMLUtils {

    public static void addWithTransposed(DenseMatrix64F X) {

        checkSquare(X);

        for (int i = 0; i < X.numCols; i++) {
            double xii = X.get(i, i);
            X.set(i, i, xii * 2);
            for (int j = i + 1; j < X.numCols; j++) {
                double xij = X.get(i, j);
                double xji = X.get(j, i);
                double x = xij + xji;
                X.set(i, j, x);
                X.set(j, i, x);
            }
        }
    }

    public static void setIdentity(DenseMatrix64F X) {
        setScaledIdentity(X, 1);
    }

    public static void setScaledIdentity(DenseMatrix64F X, double scale) {
        checkSquare(X);

        Arrays.fill(X.data, 0);
        for (int i = 0; i < X.numCols; i++) {
            X.set(i, i, scale);
        }
    }


    private static void checkSquare(DenseMatrix64F X) {
        if (X.numCols != X.numRows) {
            throw new IllegalArgumentException("matrix must be square.");
        }
    }
}
