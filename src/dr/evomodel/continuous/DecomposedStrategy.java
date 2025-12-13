package dr.evomodel.continuous;

import dr.inference.model.MatrixParameterInterface;
import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import dr.inference.model.CompoundEigenMatrix;

/**
 * Strategy for pre-decomposed strength-of-selection matrices.
 * Uses CompoundEigenMatrix which already provides eigenvalues and eigenvectors.
 *
 * A = R Λ R^{-1}, where R and Λ are provided directly.
 */
class DecomposedStrategy implements MultivariateElasticModel.SelectionMatrixStrategy {

    private final boolean isSymmetric;

    public DecomposedStrategy(boolean isSymmetric) {
        this.isSymmetric = isSymmetric;
    }

    @Override
    public BasisRepresentation computeBasis(MatrixParameterInterface param, int dim) {
        if (!(param instanceof CompoundEigenMatrix)) {
            throw new IllegalArgumentException(
                    "DecomposedStrategy requires CompoundEigenMatrix parameter");
        }

        CompoundEigenMatrix eigenMatrix = (CompoundEigenMatrix) param;

        double[] eigenvalues = eigenMatrix.getEigenValues();
        double[] eigenvectors = eigenMatrix.getEigenVectors();

        DenseMatrix64F R = MissingOps.wrap(eigenvectors, 0, dim, dim);
        DenseMatrix64F Rinv = computeInverse(R, isSymmetric);

        int nParameters = eigenvalues.length + eigenvectors.length;

        return new BasisRepresentation(dim, nParameters, eigenvalues, eigenvectors, Rinv.getData());
    }

    @Override
    public boolean isDiagonal() {
        return false;
    }

    @Override
    public boolean isSymmetric() {
        return isSymmetric;
    }

    @Override
    public boolean isBlockDiagonal() {
        return false;
    }

    /**
     * Compute R^{-1}.
     * For symmetric matrices, R is orthonormal so R^{-1} = R^T.
     */
    private static DenseMatrix64F computeInverse(DenseMatrix64F R, boolean isSymmetric) {
        if (isSymmetric) {
            DenseMatrix64F Rinv = new DenseMatrix64F(R.numRows, R.numCols);
            CommonOps.transpose(R, Rinv);
            return Rinv;
        } else {
            DenseMatrix64F Rinv = new DenseMatrix64F(R.numRows, R.numCols);
            CommonOps.invert(R, Rinv);
            return Rinv;
        }
    }
}