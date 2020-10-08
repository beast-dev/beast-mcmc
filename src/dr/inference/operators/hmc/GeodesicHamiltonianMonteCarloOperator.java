package dr.inference.operators.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.math.matrixAlgebra.EJMLUtils;
import dr.math.matrixAlgebra.SkewSymmetricMatrixExponential;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Transform;
import org.ejml.alg.dense.decomposition.TriangularSolver;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.CholeskyDecomposition;
import org.ejml.ops.CommonOps;

public class GeodesicHamiltonianMonteCarloOperator extends HamiltonianMonteCarloOperator {

    public GeodesicHamiltonianMonteCarloOperator(AdaptationMode mode, double weight, GradientWrtParameterProvider gradientProvider, Parameter parameter, Transform transform, Parameter maskParameter, Options runtimeOptions, MassPreconditioner.Type preconditioningType) {
        super(mode, weight, gradientProvider, parameter, transform, maskParameter, runtimeOptions, preconditioningType);
        this.leapFrogEngine = new GeodesicLeapFrogEngine(parameter, getDefaultInstabilityHandler(), preconditioning, mask);
    }

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
            this.positionMatrix = new DenseMatrix64F(nCols, nRows);
            this.innerProduct = new DenseMatrix64F(nCols, nCols);
            this.innerProduct2 = new DenseMatrix64F(nCols, nCols);
            this.projection = new DenseMatrix64F(nCols, nRows);
            this.momentumMatrix = new DenseMatrix64F(nCols, nRows);
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
            CommonOps.multTransB(positionMatrix, momentumMatrix, innerProduct);
            CommonOps.multTransB(momentumMatrix, momentumMatrix, innerProduct2);

            double[][] VtV = new double[2 * nCols][2 * nCols];

            for (int i = 0; i < nCols; i++) {
                VtV[i + nCols][i] = 1;
                for (int j = 0; j < nCols; j++) {
                    VtV[i][j] = innerProduct.get(i, j);
                    VtV[i + nCols][j + nCols] = innerProduct.get(i, j);
                    VtV[i][j + nCols] = -innerProduct2.get(j, i);
                }
            }

            double[] expBuffer = new double[nCols * nCols];
            CommonOps.scale(-functionalStepSize, innerProduct);
            SkewSymmetricMatrixExponential matExp1 = new SkewSymmetricMatrixExponential(nCols);
            matExp1.exponentiate(innerProduct.data, expBuffer);

            double[] expBuffer2 = new double[nCols * nCols * 4];
            SkewSymmetricMatrixExponential matExp2 = new SkewSymmetricMatrixExponential(nCols * 2); //TODO: better matrix exponential
            DenseMatrix64F VtVmat = new DenseMatrix64F(VtV);
            CommonOps.scale(functionalStepSize, VtVmat);
            matExp2.exponentiate(VtVmat.data, expBuffer2);

            DenseMatrix64F X = new DenseMatrix64F(nCols * 2, nCols * 2);
            DenseMatrix64F Y = new DenseMatrix64F(nCols * 2, nCols * 2);

            for (int i = 0; i < nCols; i++) {
                for (int j = 0; j < nCols; j++) {
                    X.set(i, j, expBuffer[i * nCols + j]);
                    X.set(i + nCols, j + nCols, expBuffer[i * nCols + j]);
                }
            }
            Y.setData(expBuffer2);

            DenseMatrix64F Z = new DenseMatrix64F(nCols * 2, nCols * 2);

            CommonOps.mult(Y, X, Z);

            DenseMatrix64F PM = new DenseMatrix64F(nCols * 2, nRows);
            for (int i = 0; i < nRows; i++) {
                for (int j = 0; j < nCols; j++) {
                    PM.set(j, i, positionMatrix.get(j, i));
                    PM.set(j + nCols, i, momentumMatrix.get(j, i));
                }
            }

            DenseMatrix64F W = new DenseMatrix64F(2 * nCols, nRows);
            CommonOps.transpose(Z);
            CommonOps.mult(Z, PM, W);

            for (int i = 0; i < nRows; i++) {
                for (int j = 0; j < nCols; j++) {
                    positionMatrix.set(j, i, W.get(j, i));
                    momentumMatrix.set(j, i, W.get(j + nCols, i));
                }
            }

            //TODO: only run chunk below occasionally
            CommonOps.multTransB(positionMatrix, positionMatrix, innerProduct);
            CholeskyDecomposition cholesky = DecompositionFactory.chol(nCols, true);
            cholesky.decompose(innerProduct);
            TriangularSolver.invertLower(innerProduct.data, nCols);
            CommonOps.mult(innerProduct, positionMatrix, projection);
            System.arraycopy(projection.data, 0, positionMatrix.data, 0, positionMatrix.data.length);

            System.arraycopy(positionMatrix.data, 0, position, 0, position.length);
            System.arraycopy(momentumMatrix.data, 0, momentum.getBuffer(), momentum.getOffset(), momentum.getDim());

            matrixParameter.setAllParameterValuesQuietly(position, 0);
            matrixParameter.fireParameterChangedEvent();


        }

        @Override
        public void projectMomentum(double[] momentum, double[] position) {
            positionMatrix.setData(position);
            momentumMatrix.setData(momentum);

            CommonOps.multTransB(positionMatrix, momentumMatrix, innerProduct);
            EJMLUtils.addWithTransposed(innerProduct);

            CommonOps.mult(0.5, innerProduct, positionMatrix, projection);
            CommonOps.subtractEquals(momentumMatrix, projection);
        }
    }
}
