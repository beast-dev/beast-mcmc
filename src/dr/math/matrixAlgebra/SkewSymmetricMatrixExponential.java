package dr.math.matrixAlgebra;

import dr.math.MathUtils;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.SingularValueDecomposition;
import org.ejml.ops.CommonOps;

// Derived from:
//  Cardoso, João R., and F. Silva Leite. "Exponentials of skew-symmetric matrices and logarithms of orthogonal matrices."
//      Journal of computational and applied mathematics 233.11 (2010): 2867-2875.

public class SkewSymmetricMatrixExponential {

    private final int dim;
    private final DenseMatrix64F S;
    private final DenseMatrix64F B1;
    private final DenseMatrix64F B2;
    private final DenseMatrix64F B3;
    private final DenseMatrix64F P1;
    private final DenseMatrix64F P2;
    private final DenseMatrix64F buffer;
    //TODO: do I need so many buffers?


    public SkewSymmetricMatrixExponential(int dim) {
        this.dim = dim;
        this.S = new DenseMatrix64F(dim, dim);
        this.B1 = new DenseMatrix64F(dim, dim);
        this.B2 = new DenseMatrix64F(dim, dim);
        this.B3 = new DenseMatrix64F(dim, dim);
        this.P1 = new DenseMatrix64F(dim, dim);
        this.P2 = new DenseMatrix64F(dim, dim);
        this.buffer = new DenseMatrix64F(dim, dim);


    }

    public void exponentiate(double[] src, double[] dest) {
        int dimSquared = dim * dim;
        if (src.length != dimSquared || dest.length != dimSquared) {
            throw new IllegalArgumentException("At least one matrix is of wrong dimension.");
        }

        System.arraycopy(src, 0, S.data, 0, dimSquared);


        SingularValueDecomposition svd = DecompositionFactory.svd(dim, dim, false, false, true);
        svd.decompose(S);

        double[] singularValues = svd.getSingularValues();
        double norm = MathUtils.maximum(singularValues);
        double scale = 1.0;
        int k = 0;
        while (norm > 1) {
            k++;
            norm *= 0.5;
            scale *= 0.5;
        }

        CommonOps.scale(scale, S);
        CommonOps.mult(S, S, B1);
        CommonOps.mult(B1, B1, B2);
        CommonOps.mult(B1, B2, B3);

        EJMLUtils.setScaledIdentity(P1, 17297280);


        CommonOps.addEquals(P1, 1995840, B1);
        CommonOps.addEquals(P1, 25200, B2);
        CommonOps.addEquals(P1, 56, B3);

        EJMLUtils.setScaledIdentity(buffer, 8648640);
        CommonOps.addEquals(buffer, 277200, B1);
        CommonOps.addEquals(buffer, 1512, B2);
        CommonOps.addEquals(buffer, B3);
        CommonOps.mult(S, buffer, P2);

        CommonOps.subtract(P1, P2, B1); //B1 just a buffer
        CommonOps.add(P1, P2, B2); //B2 just a buffer

        CommonOps.solve(B1, B2, S);

        for (int i = 0; i < k; i++) {
            System.arraycopy(S.data, 0, B1.data, 0, S.data.length);
            CommonOps.mult(B1, B1, S);
        }

        System.arraycopy(S.data, 0, dest, 0, S.data.length);
    }


}