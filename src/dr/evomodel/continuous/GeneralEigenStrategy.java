package dr.evomodel.continuous;

import dr.inference.model.MatrixParameterInterface;
import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.Complex64F;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.EigenDecomposition;
import org.ejml.ops.CommonOps;
import org.ejml.ops.EigenOps;

/**
 * Strategy for general strength-of-selection matrices.
 * Performs eigen-decomposition using EJML to obtain A = R Λ R^{-1}.
 *
 * Validates that all eigenvalues are positive and real.
 */
class GeneralEigenStrategy implements MultivariateElasticModel.SelectionMatrixStrategy {

    private final boolean isSymmetric;

    public GeneralEigenStrategy(boolean isSymmetric) {
        this.isSymmetric = isSymmetric;
    }

    @Override
    public BasisRepresentation computeBasis(MatrixParameterInterface param, int dim) {

        DenseMatrix64F A = MissingOps.wrap(param);

        EigenDecomposition<DenseMatrix64F> eigenDecomp =
                DecompositionFactory.eig(dim, true, isSymmetric);

        if (!eigenDecomp.decompose(A)) {
            throw new RuntimeException("Eigen decomposition of selection matrix A failed.");
        }

        double[] eigenvalues = extractEigenvalues(eigenDecomp, dim);
        double[] eigenvectors = extractEigenvectors(eigenDecomp);

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
     * Extract eigenvalues and validate they are positive and real.
     */
    private static double[] extractEigenvalues(EigenDecomposition<DenseMatrix64F> eigenDecomp,
                                               int dim) {
        double[] eigenvalues = new double[dim];

        for (int i = 0; i < dim; i++) {
            Complex64F eigenvalue = eigenDecomp.getEigenvalue(i);

            if (!eigenvalue.isReal()) {
                throw new IllegalStateException(
                        "Selection strength matrix A must have only real eigenvalues. " +
                                "Found complex eigenvalue at index " + i + ": " + eigenvalue);
            }

            if (eigenvalue.real <= 0.0) {
                throw new IllegalStateException(
                        "Selection strength matrix A must have only positive eigenvalues. " +
                                "Found eigenvalue " + eigenvalue.real + " at index " + i);
            }

            eigenvalues[i] = eigenvalue.real;
        }

        return eigenvalues;
    }

    /**
     * Extract eigenvector matrix as flattened row-major array.
     */
    private static double[] extractEigenvectors(EigenDecomposition<DenseMatrix64F> eigenDecomp) {
        DenseMatrix64F V = EigenOps.createMatrixV(eigenDecomp);
        return V.getData();  // EJML uses row-major order
    }

    /**
     * Compute R^{-1}.
     * For symmetric matrices, R is orthonormal so R^{-1} = R^T.
     */
    private static DenseMatrix64F computeInverse(DenseMatrix64F R, boolean isSymmetric) {
        DenseMatrix64F Rinv = new DenseMatrix64F(R.numRows, R.numCols);

        if (isSymmetric) {
            CommonOps.transpose(R, Rinv);
        } else {
            CommonOps.invert(R, Rinv);
        }

        return Rinv;
    }
}