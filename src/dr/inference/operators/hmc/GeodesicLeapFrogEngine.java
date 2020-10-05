package dr.inference.operators.hmc;

import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.EJMLUtils;
import dr.math.matrixAlgebra.WrappedVector;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.EigenDecomposition;
import org.ejml.ops.CommonOps;


public class GeodesicLeapFrogEngine extends HamiltonianMonteCarloOperator.LeapFrogEngine.Default {

    private final MatrixParameterInterface matrixParameter;
    private final DenseMatrix64F positionMatrix;
    private final DenseMatrix64F innerProduct;
    private final DenseMatrix64F innerProduct2;
    private final DenseMatrix64F projection;
    private final DenseMatrix64F momentumMatrix;


    GeodesicLeapFrogEngine(Parameter parameter, HamiltonianMonteCarloOperator.InstabilityHandler instabilityHandler,
                           MassPreconditioner preconditioning, double[] mask) {
        super(parameter, instabilityHandler, preconditioning, mask);
        this.matrixParameter = (MatrixParameterInterface) parameter;
        int nRows = matrixParameter.getRowDimension();
        int nCols = matrixParameter.getColumnDimension();
        this.positionMatrix = new DenseMatrix64F(nRows, nCols);
        this.innerProduct = new DenseMatrix64F(nCols, nCols);
        this.innerProduct2 = new DenseMatrix64F(nCols, nCols);
        this.projection = new DenseMatrix64F(nRows, nCols);
        this.momentumMatrix = new DenseMatrix64F(nRows, nCols);
    }

    @Override
    public void updateMomentum(double[] position, double[] momentum, double[] gradient,
                               double functionalStepSize) throws HamiltonianMonteCarloOperator.NumericInstabilityException {
        super.updateMomentum(position, momentum, gradient, functionalStepSize);
        projectMomentum(momentum, position);

    }

    @Override
    public void updatePosition(double[] position, WrappedVector momentum,
                               double functionalStepSize) throws HamiltonianMonteCarloOperator.NumericInstabilityException {

        positionMatrix.setData(position);
        System.arraycopy(momentum.getBuffer(), momentum.getOffset(), momentumMatrix.data, 0, momentum.getDim());
        CommonOps.multTransA(positionMatrix, momentumMatrix, innerProduct);

        EigenDecomposition eig = DecompositionFactory.eig(innerProduct.numCols, true);
        eig.decompose(innerProduct);

        //TODO


    }

    @Override
    public void projectMomentum(double[] momentum, double[] position) {
        positionMatrix.setData(position);
        momentumMatrix.setData(momentum);

        CommonOps.multTransA(positionMatrix, momentumMatrix, innerProduct);
        EJMLUtils.addWithTransposed(innerProduct);

        CommonOps.mult(0.5, positionMatrix, innerProduct, projection);
        CommonOps.subtractEquals(momentumMatrix, projection);
    }
}

