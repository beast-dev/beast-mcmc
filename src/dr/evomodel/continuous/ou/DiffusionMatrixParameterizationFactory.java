package dr.evomodel.continuous.ou;

import dr.inference.model.CompoundSymmetricMatrix;
import dr.inference.model.MatrixParameterInterface;

/**
 * Factory for OU diffusion-matrix parametrizations.
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
