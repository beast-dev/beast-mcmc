package dr.evomodel.treedatalikelihood.hmc;

import dr.inference.model.*;
import dr.util.MatrixInnerProductTransform;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */

public class FullCorrelationPrecisionGradient extends CorrelationPrecisionGradient {

    private final MatrixParameterInterface decomposedMatrix;
    private static final RuntimeException PARAMETER_EXCEPTION = new RuntimeException("off-diagonal parameter must be a mask of a inner product transform.");

    public FullCorrelationPrecisionGradient(GradientWrtPrecisionProvider gradientWrtPrecisionProvider, Likelihood likelihood, MatrixParameterInterface parameter) {
        super(gradientWrtPrecisionProvider, likelihood, parameter);

        // TODO: this is super messy but get's the job done (for now). Maybe create a class just to wrap these different transforms/masks
        Parameter correlationParameter = compoundSymmetricMatrix.getOffDiagonalParameter();
        if (correlationParameter instanceof MaskedParameter) {
            MaskedParameter maskedCorrelation = (MaskedParameter) correlationParameter;
            if (maskedCorrelation.getUnmaskedParameter() instanceof TransformedMultivariateParameter) {
                TransformedMultivariateParameter transformedParameter =
                        (TransformedMultivariateParameter) maskedCorrelation.getUnmaskedParameter();
                if (transformedParameter.getTransform() instanceof MatrixInnerProductTransform) { //TODO: chain rul in transforms
                    decomposedMatrix = (MatrixParameterInterface) transformedParameter.getUntransformedParameter();
                } else {
                    throw PARAMETER_EXCEPTION;
                }
            } else {
                throw PARAMETER_EXCEPTION;
            }
        } else {
            throw PARAMETER_EXCEPTION;
        }

    }


    @Override
    public int getDimension() {
        return getDimensionFull();
    }

    @Override
    public Parameter getParameter() {
        return decomposedMatrix;
    }

    @Override
    double[] getGradientParameter(double[] gradient) {
        int dim = decomposedMatrix.getRowDimension();

        double[] correlationGradient = compoundSymmetricMatrix.updateGradientFullOffDiagonal(gradient);
        double[] decomposedGradient = new double[gradient.length];

        DenseMatrix64F corGradMat = DenseMatrix64F.wrap(dim, dim, correlationGradient);
        DenseMatrix64F decompMat = DenseMatrix64F.wrap(dim, dim, decomposedMatrix.getParameterValues());
        DenseMatrix64F decompGradMat = DenseMatrix64F.wrap(dim, dim, decomposedGradient);

        CommonOps.multTransA(corGradMat, decompMat, decompGradMat);

        CommonOps.scale(2.0, decompGradMat);

        return decomposedGradient;
    }


}
