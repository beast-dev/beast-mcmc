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
        return create(matrixParameter, 0.0);
    }

    /**
     * @param diffusionFloor constant added to every diagonal entry of the diffusion matrix the
     *                       parameterization reports ({@code Q + floor * I}); the pull-back of
     *                       gradients to the underlying parameters is unaffected because the
     *                       added term is constant.
     */
    public static DiffusionMatrixParameterization create(final MatrixParameterInterface matrixParameter,
                                                         final double diffusionFloor) {
        if (matrixParameter instanceof CompoundSymmetricMatrix) {
            return new CompoundSymmetricDiffusionMatrixParameterization(
                    (CompoundSymmetricMatrix) matrixParameter, diffusionFloor);
        }
        return new DenseDiffusionMatrixParameterization(matrixParameter, diffusionFloor);
    }
}
