package dr.evomodel.continuous.ou;

import dr.inference.model.CompoundSymmetricMatrix;
import dr.inference.model.Parameter;

/**
 * OU diffusion parametrization for covariance matrices represented by a
 * diagonal variance vector together with an off-diagonal correlation (or
 * Cholesky-correlation) vector through {@link CompoundSymmetricMatrix}.
 */
public final class CompoundSymmetricDiffusionMatrixParameterization implements DiffusionMatrixParameterization {

    private final CompoundSymmetricMatrix matrixParameter;
    private final int dimension;
    private final double diffusionFloor;

    public CompoundSymmetricDiffusionMatrixParameterization(final CompoundSymmetricMatrix matrixParameter) {
        this(matrixParameter, 0.0);
    }

    public CompoundSymmetricDiffusionMatrixParameterization(final CompoundSymmetricMatrix matrixParameter,
                                                            final double diffusionFloor) {
        if (matrixParameter == null) {
            throw new IllegalArgumentException("matrixParameter must not be null");
        }
        if (!Double.isFinite(diffusionFloor) || diffusionFloor < 0.0) {
            throw new IllegalArgumentException("diffusionFloor must be finite and non-negative");
        }
        this.matrixParameter = matrixParameter;
        this.dimension = matrixParameter.getRowDimension();
        this.diffusionFloor = diffusionFloor;
    }

    public double getDiffusionFloor() {
        return diffusionFloor;
    }

    @Override
    public CompoundSymmetricMatrix getMatrixParameter() {
        return matrixParameter;
    }

    @Override
    public int getDimension() {
        return dimension;
    }

    @Override
    public void fillDiffusionMatrixFlat(final double[] out) {
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                out[rowOffset + j] = matrixParameter.getParameterValue(i, j);
            }
            out[rowOffset + i] += diffusionFloor;
        }
    }

    @Override
    public boolean supportsParameter(final Parameter parameter) {
        return parameter == matrixParameter
                || parameter == matrixParameter.getDiagonalParameter()
                || parameter == matrixParameter.getOffDiagonalParameter()
                || parameter == matrixParameter.getUntransformedOffDiagonalParameter()
                || parameter == matrixParameter.getUntransformedCompoundParameter();
    }

    @Override
    public double[] pullBackGradient(final Parameter parameter, final double[] denseGradient) {
        final double[] diagonalGradient = matrixParameter.updateGradientDiagonal(denseGradient);
        final double[] offDiagonalGradient = matrixParameter.updateGradientOffDiagonal(denseGradient);

        if (parameter == matrixParameter) {
            return join(diagonalGradient, offDiagonalGradient);
        }
        if (parameter == matrixParameter.getDiagonalParameter()) {
            return diagonalGradient;
        }
        if (parameter == matrixParameter.getOffDiagonalParameter()
                || parameter == matrixParameter.getUntransformedOffDiagonalParameter()) {
            return offDiagonalGradient;
        }
        if (parameter == matrixParameter.getUntransformedCompoundParameter()) {
            return join(diagonalGradient, offDiagonalGradient);
        }
        throw new IllegalArgumentException("Unsupported compound-symmetric diffusion parameter");
    }

    private static double[] join(final double[] diagonalGradient,
                                 final double[] offDiagonalGradient) {
        final double[] out = new double[diagonalGradient.length + offDiagonalGradient.length];
        System.arraycopy(diagonalGradient, 0, out, 0, diagonalGradient.length);
        System.arraycopy(offDiagonalGradient, 0, out, diagonalGradient.length, offDiagonalGradient.length);
        return out;
    }
}
