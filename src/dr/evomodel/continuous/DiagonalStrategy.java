package dr.evomodel.continuous;

import dr.inference.model.DiagonalMatrix;
import dr.inference.model.MatrixParameterInterface;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * Strategy for diagonal strength-of-selection matrices.
 * A = diag(λ_1, ..., λ_d)
 *
 * Basis is simply the identity matrix.
 */
class DiagonalStrategy implements MultivariateElasticModel.SelectionMatrixStrategy {

    @Override
    public BasisRepresentation computeBasis(MatrixParameterInterface param, int dim) {
        if (!(param instanceof DiagonalMatrix)) {
            throw new IllegalArgumentException(
                    "DiagonalStrategy requires DiagonalMatrix parameter");
        }

        DiagonalMatrix diagMatrix = (DiagonalMatrix) param;
        double[] eigenvalues = diagMatrix.getDiagonalParameter().getParameterValues();

        // Basis is identity: R = I, R^{-1} = I
//        DenseMatrix64F identity = CommonOps.identity(dim);
//        double[] identityVector = createIdentityVector(dim);

        return new BasisRepresentation(
                dim,
                eigenvalues.length,
                eigenvalues,
                null,  // avoids unnecessary storage
                null // avoids unnecessary storage
        );
    }

    @Override
    public boolean isDiagonal() {
        return true;
    }

    @Override
    public boolean isSymmetric() {
        return true;  // Diagonal matrices are symmetric
    }

    @Override
    public boolean isBlockDiagonal() {
        return false;
    }

//    private static double[] createIdentityVector(int dim) {
//        double[] identity = new double[dim * dim];
//        for (int i = 0; i < dim; i++) {
//            identity[i * dim + i] = 1.0;
//        }
//        return identity;
//    }
}