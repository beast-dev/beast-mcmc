package dr.inference.operators.hmc;

import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.WrappedVector;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;


public class GeodesicLeapFrogEngine extends HamiltonianMonteCarloOperator.LeapFrogEngine.Default {

    private final MatrixParameterInterface matrixParameter;
    private final DenseMatrix64F positionMatrix;
    private final DenseMatrix64F innerProduct;
    private final DenseMatrix64F projection;
    private final DenseMatrix64F momentumMatrix;


    GeodesicLeapFrogEngine(Parameter parameter, HamiltonianMonteCarloOperator.InstabilityHandler instabilityHandler,
                           MassPreconditioner preconditioning, double[] mask) {
        super(parameter, instabilityHandler, preconditioning, mask);
        this.matrixParameter = (MatrixParameterInterface) parameter;
        this.positionMatrix = new DenseMatrix64F(matrixParameter.getRowDimension(), matrixParameter.getColumnDimension());
        this.innerProduct = new DenseMatrix64F(matrixParameter.getRowDimension(), matrixParameter.getRowDimension());
        this.projection = new DenseMatrix64F(matrixParameter.getRowDimension(), 1);
        this.momentumMatrix = new DenseMatrix64F(matrixParameter.getRowDimension(), 1);
    }

    @Override
    public void updateMomentum(double[] position, double[] momentum, double[] gradient,
                               double functionalStepSize) throws HamiltonianMonteCarloOperator.NumericInstabilityException {


    }

    @Override
    public void updatePosition(double[] position, WrappedVector momentum,
                               double functionalStepSize) throws HamiltonianMonteCarloOperator.NumericInstabilityException {

        //TODO:

    }

    @Override
    public void projectMomentum(double[] momentum, double[] position) {
        positionMatrix.setData(position);
        CommonOps.multTransB(positionMatrix, positionMatrix, innerProduct);
        for (int i = 0; i < matrixParameter.getRowDimension(); i++) {
            innerProduct.set(i, i, innerProduct.get(i, i) - 1);
        }

        momentumMatrix.setData(momentum);

        CommonOps.mult(positionMatrix, momentumMatrix, projection);
        for (int i = 0; i < projection.numRows; i++) {
            projection.data[i] *= -1;
        }
        System.arraycopy(projection.data, 0, momentum, 0, projection.data.length);
    }
}

