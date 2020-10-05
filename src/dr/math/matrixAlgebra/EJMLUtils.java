package dr.math.matrixAlgebra;

import org.ejml.data.DenseMatrix64F;

public class EJMLUtils {

    public static void addWithTransposed(DenseMatrix64F X) {

        if (X.numCols != X.numRows) {
            throw new IllegalArgumentException("matrix must be square.");
        }

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
}
