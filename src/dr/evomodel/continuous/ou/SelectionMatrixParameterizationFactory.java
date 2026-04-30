package dr.evomodel.continuous.ou;

import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.OrthogonalMatrixProvider;

/**
 * Factory for selection-matrix parametrizations.
 */
public final class SelectionMatrixParameterizationFactory {

    private SelectionMatrixParameterizationFactory() {
    }

    public static SelectionMatrixParameterization create(final MatrixParameterInterface matrixParameter) {
        if (matrixParameter instanceof AbstractBlockDiagonalTwoByTwoMatrixParameter) {
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter =
                    (AbstractBlockDiagonalTwoByTwoMatrixParameter) matrixParameter;
            if (blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider) {
                return new OrthogonalBlockDiagonalSelectionMatrixParameterization(
                        blockParameter,
                        (OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter());
            }
        }
        return new DenseSelectionMatrixParameterization(matrixParameter);
    }
}
