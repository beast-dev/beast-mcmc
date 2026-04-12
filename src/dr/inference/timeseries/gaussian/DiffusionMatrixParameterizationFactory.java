package dr.inference.timeseries.gaussian;

import dr.inference.model.CompoundSymmetricMatrix;
import dr.inference.model.MatrixParameterInterface;

/**
 * Factory for diffusion-matrix parametrizations.
 */
public final class DiffusionMatrixParameterizationFactory {

    private DiffusionMatrixParameterizationFactory() {
    }

    public static DiffusionMatrixParameterization create(final MatrixParameterInterface matrixParameter) {
        if (matrixParameter instanceof CompoundSymmetricMatrix) {
            return new CompoundSymmetricDiffusionMatrixParameterization(
                    (CompoundSymmetricMatrix) matrixParameter);
        }
        return new DenseDiffusionMatrixParameterization(matrixParameter);
    }
}
