package dr.math.matrixAlgebra;

import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.SingularValueDecomposition;
import org.ejml.ops.CommonOps;

import java.util.Arrays;

public class EJMLUtils {

    private static final Boolean DEBUG = false;

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

    public static DenseMatrix64F computeRobustAdjugate(DenseMatrix64F X) throws IllegalDimension {
        // algorithm taken from: Stewart, G. W. "On the adjugate matrix." Linear Algebra and its Applications 283.1-3 (1998): 151-164.

        if (DEBUG) {
            System.out.println(X);
        }


        checkSquare(X);

        int dim = X.numRows;

        SingularValueDecomposition svd = DecompositionFactory.svd(dim, dim,
                true, true, true);

        svd.decompose(X);
        DenseMatrix64F Ut = new DenseMatrix64F(dim, dim);
        DenseMatrix64F V = new DenseMatrix64F(dim, dim);

        svd.getU(Ut, true);
        svd.getV(V, false);
        double[] svs = svd.getSingularValues();

        double du = CommonOps.det(Ut);
        double dv = CommonOps.det(V); //TODO: should be same as du if matrix is symmetric

        // gamma_i is just the product of the singular values excluding sv_i
        double[] gamma = new double[dim];
        double[] forwardProds = new double[dim];
        double[] backwardProds = new double[dim];
        double fp = 1;
        double bp = 1;
        for (int i = 0; i < dim; i++) {
            fp *= svs[i];
            forwardProds[i] = fp;

            int backInd = dim - i - 1;
            bp *= svs[backInd];
            backwardProds[backInd] = bp;
        }

        gamma[0] = backwardProds[1];
        gamma[dim - 1] = forwardProds[dim - 2];

        for (int i = 1; i < dim - 1; i++) {
            gamma[i] = forwardProds[i - 1] * backwardProds[i + 1];
        }


        // Ut = Diag(gamma) * Ut
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                Ut.set(i, j, Ut.get(i, j) * gamma[i]);
            }
        }

        DenseMatrix64F adjugate = new DenseMatrix64F(dim, dim);
        CommonOps.mult(V, Ut, adjugate);
        CommonOps.scale(du * dv, adjugate);

        if (DEBUG) {
            System.out.println(adjugate);
        }

        return adjugate;
    }

}
