package dr.math.geodesics;

import dr.math.matrixAlgebra.SkewSymmetricMatrixExponential;
import org.ejml.alg.dense.decomposition.TriangularSolver;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.CholeskyDecomposition;
import org.ejml.ops.CommonOps;


public class StiefelManifold implements Manifold {

    private final int rowDim;
    private final int colDim;
    
    private final DenseMatrix64F innerProduct1;
    private final DenseMatrix64F innerProduct2;
    public StiefelManifold(int rowDim, int colDim) {
        this.rowDim = rowDim;
        this.colDim = colDim;
        
        this.innerProduct1 = new DenseMatrix64F(colDim, colDim);
        this.innerProduct2 = new DenseMatrix64F(colDim, colDim);


    }

    @Override
    public void projectTangent(double[] tangent, double[] point) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void projectPoint(double[] point) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void geodesic(double[] point, double[] velocity, double t) {
        DenseMatrix64F positionMatrix = DenseMatrix64F.wrap(colDim, rowDim, point);
        DenseMatrix64F velocityMatrix = DenseMatrix64F.wrap(colDim, rowDim, velocity);

        CommonOps.multTransB(positionMatrix, velocityMatrix, innerProduct1);
        CommonOps.multTransB(velocityMatrix, velocityMatrix, innerProduct2);

        double[][] VtV = new double[2 * colDim][2 * colDim];

        for (int i = 0; i < colDim; i++) {
            VtV[i + colDim][i] = 1;
            for (int j = 0; j < colDim; j++) {
                VtV[i][j] = innerProduct1.get(i, j);
                VtV[i + colDim][j + colDim] = innerProduct1.get(i, j);
                VtV[i][j + colDim] = -innerProduct2.get(j, i);
            }
        }

        double[] expBuffer = new double[colDim * colDim];
        CommonOps.scale(-t, innerProduct1);
        SkewSymmetricMatrixExponential matExp1 = new SkewSymmetricMatrixExponential(colDim);
        matExp1.exponentiate(innerProduct1.data, expBuffer);

        double[] expBuffer2 = new double[colDim * colDim * 4];
        SkewSymmetricMatrixExponential matExp2 = new SkewSymmetricMatrixExponential(colDim * 2); //TODO: better matrix exponential
        DenseMatrix64F VtVmat = new DenseMatrix64F(VtV);
        CommonOps.scale(t, VtVmat);
        matExp2.exponentiate(VtVmat.data, expBuffer2);

        DenseMatrix64F X = new DenseMatrix64F(colDim * 2, colDim * 2);
        DenseMatrix64F Y = new DenseMatrix64F(colDim * 2, colDim * 2);

        for (int i = 0; i < colDim; i++) {
            for (int j = 0; j < colDim; j++) {
                X.set(i, j, expBuffer[i * colDim + j]);
                X.set(i + colDim, j + colDim, expBuffer[i * colDim + j]);
            }
        }
        Y.setData(expBuffer2);

        DenseMatrix64F Z = new DenseMatrix64F(colDim * 2, colDim * 2);

        CommonOps.mult(Y, X, Z);

        DenseMatrix64F PM = new DenseMatrix64F(colDim * 2, rowDim);
        for (int i = 0; i < rowDim; i++) {
            for (int j = 0; j < colDim; j++) {
                PM.set(j, i, positionMatrix.get(j, i));
                PM.set(j + colDim, i, velocityMatrix.get(j, i));
            }
        }

        DenseMatrix64F W = new DenseMatrix64F(2 * colDim, rowDim);
        CommonOps.transpose(Z);
        CommonOps.mult(Z, PM, W);

        for (int i = 0; i < rowDim; i++) {
            for (int j = 0; j < colDim; j++) {
                positionMatrix.set(j, i, W.get(j, i));
                velocityMatrix.set(j, i, W.get(j + colDim, i));
            }
        }

        //TODO: only run chunk below occasionally
        CommonOps.multTransB(positionMatrix, positionMatrix, innerProduct1);
        CholeskyDecomposition cholesky = DecompositionFactory.chol(colDim, true);
        cholesky.decompose(innerProduct1);
        TriangularSolver.invertLower(innerProduct1.data, colDim);

        DenseMatrix64F projection = new DenseMatrix64F(colDim, rowDim);

        CommonOps.mult(innerProduct1, positionMatrix, projection);

        double sse = 0;
        for (int i = 0; i < positionMatrix.data.length; i++) {
            double diff = projection.data[i] - positionMatrix.data[i];
            sse += diff * diff;
        }

        if (sse / point.length > 1e-2) { //TODO: actually figure out if I want this
            System.err.println("unstable"); //TODO: REMOVE
            throw new RuntimeException(); //TODO: need to handle this in HMC, don't throw runtime
        }


        System.arraycopy(projection.data, 0, point, 0, point.length);
        System.arraycopy(velocityMatrix.data, 0, velocity, 0, velocity.length);
        
    }
}
