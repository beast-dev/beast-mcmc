package dr.evomodel.continuous.ou;

import dr.inference.model.MatrixParameterInterface;

/**
 * Shared validation and accessors for OU selection-matrix parametrizations.
 */
public abstract class AbstractSelectionMatrixParameterization implements SelectionMatrixParameterization {

    private final MatrixParameterInterface matrixParameter;
    private final int dimension;

    protected AbstractSelectionMatrixParameterization(final MatrixParameterInterface matrixParameter) {
        if (matrixParameter == null) {
            throw new IllegalArgumentException("matrixParameter must not be null");
        }
        if (matrixParameter.getRowDimension() != matrixParameter.getColumnDimension()) {
            throw new IllegalArgumentException("selection matrix must be square");
        }
        this.matrixParameter = matrixParameter;
        this.dimension = matrixParameter.getRowDimension();
    }

    @Override
    public MatrixParameterInterface getMatrixParameter() {
        return matrixParameter;
    }

    @Override
    public int getDimension() {
        return dimension;
    }

    protected final void checkFlatSquare(final double[] matrix, final String label) {
        if (matrix == null || matrix.length != dimension * dimension) {
            throw new IllegalArgumentException(
                    label + " must have length " + (dimension * dimension)
                            + " for a " + dimension + "x" + dimension + " matrix");
        }
    }

    protected final void checkVector(final double[] vector, final String label) {
        if (vector == null || vector.length < dimension) {
            throw new IllegalArgumentException(
                    label + " must have length at least " + dimension);
        }
    }
}
