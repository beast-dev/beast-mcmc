package dr.inference.operators.hmc;

import dr.evomodel.substmodel.ComplexColtEigenSystem;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.EJMLUtils;
import dr.math.matrixAlgebra.WrappedVector;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;


public class GeodesicLeapFrogEngine extends HamiltonianMonteCarloOperator.LeapFrogEngine.Default {

    private final MatrixParameterInterface matrixParameter;
    private final DenseMatrix64F positionMatrix;
    private final DenseMatrix64F innerProduct;
    private final DenseMatrix64F innerProduct2;
    private final DenseMatrix64F projection;
    private final DenseMatrix64F momentumMatrix;
    private final int nRows;
    private final int nCols;


    GeodesicLeapFrogEngine(Parameter parameter, HamiltonianMonteCarloOperator.InstabilityHandler instabilityHandler,
                           MassPreconditioner preconditioning, double[] mask) {
        super(parameter, instabilityHandler, preconditioning, mask);
        this.matrixParameter = (MatrixParameterInterface) parameter;
        this.nRows = matrixParameter.getRowDimension();
        this.nCols = matrixParameter.getColumnDimension();
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
        CommonOps.multTransA(momentumMatrix, momentumMatrix, innerProduct2);

        double[][] XtV = new double[nCols][nCols];
        double[][] VtV = new double[2 * nCols][2 * nCols];

        for (int i = 0; i < nCols; i++) {
            VtV[i + nCols][i] = 1;
            for (int j = 0; j < nCols; j++) {
                XtV[i][j] = innerProduct.get(i, j);
                VtV[i][j] = innerProduct.get(i, j);
                VtV[i + nCols][j + nCols] = innerProduct.get(i, j);
                VtV[i][j + nCols] = innerProduct2.get(i, j);
            }
        }

        double[] expBuffer = new double[nCols * nCols];


        ComplexColtEigenSystem eigSystem = new ComplexColtEigenSystem(nCols);
        EigenDecomposition eigDecomposition = eigSystem.decomposeMatrix(XtV);
        eigSystem.computeExponential(eigDecomposition, functionalStepSize, expBuffer);


        double[] expBuffer2 = new double[nCols * nCols * 4];

        ComplexColtEigenSystem eigSystem2 = new ComplexColtEigenSystem(nCols * 2);
        EigenDecomposition eigDecomposition2 = eigSystem2.decomposeMatrix(VtV);
        eigSystem2.computeExponential(eigDecomposition2, functionalStepSize, expBuffer2);

        DenseMatrix64F X = new DenseMatrix64F(nCols * 2, nCols * 2);
        DenseMatrix64F Y = new DenseMatrix64F(nCols * 2, nCols * 2);

        X.setData(expBuffer);
        Y.setData(expBuffer2);

        DenseMatrix64F Z = new DenseMatrix64F(nCols * 2, nCols * 2);

        System.out.println("A"); //This gets printed

        CommonOps.mult(Y, X, Z);

        System.out.println("B"); //This does not get printed

        DenseMatrix64F PM = new DenseMatrix64F(nRows, nCols * 2);
        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCols; j++) {
                PM.set(i, j, positionMatrix.get(i, j));
                PM.set(i, j + nCols, momentumMatrix.get(i, j));
            }
        }

        DenseMatrix64F W = new DenseMatrix64F(nRows, 2 * nCols);
        CommonOps.mult(PM, Z, W);

        DenseMatrix64F newPosition = new DenseMatrix64F(nRows, nCols);
        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCols; i++) {
                newPosition.set(i, j, W.get(i, j));
            }
        }

        System.arraycopy(newPosition.data, 0, position, 0, position.length);


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

