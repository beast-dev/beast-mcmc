package dr.inference.operators.hmc;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.math.MathUtils;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.RobustEigenDecomposition;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Transform;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public interface MassPreconditioner {

    WrappedVector drawInitialMomentum();

    double getVelocity(int index, ReadableVector momentum);

    void storeSecant(ReadableVector gradient, ReadableVector position);

    abstract class Base implements MassPreconditioner {
        // TODO Remove
    }

    class NoPreconditioning extends  Base implements MassPreconditioner {

        final int dim;

        NoPreconditioning(int dim) { this.dim = dim; }

        @Override
        public WrappedVector drawInitialMomentum() {

            double[] momentum = new double[dim];

            for (int i = 0; i < dim; ++i) {
                momentum[i] = MathUtils.nextGaussian();
            }

            return new WrappedVector.Raw(momentum);
        }

        @Override
        public double getVelocity(int i, ReadableVector momentum) {
            return momentum.get(i);
        }

        @Override
        public void storeSecant(ReadableVector gradient, ReadableVector position) { }
    }

    abstract class HessianBased extends Base implements MassPreconditioner {
        final protected int dim;
        final protected HessianWrtParameterProvider hessian;
        final protected Transform transform;
        final double[] inverseMass;

        // TODO Should probably make a TransformedHessian so that this class does not need to know about transformations
        HessianBased(HessianWrtParameterProvider hessian,
                     Transform transform) {
            
            this.dim = hessian.getDimension();
            this.hessian = hessian;
            this.transform = transform;

            this.inverseMass = computeInverseMass();
        }

        public void storeSecant(ReadableVector gradient, ReadableVector position) {
            // Do nothing
        }

        abstract protected double[] computeInverseMass();
    }

    class DiagonalPreconditioning extends HessianBased {


        DiagonalPreconditioning(HessianWrtParameterProvider hessian,
                                Transform transform) {
            super(hessian, transform);
        }

        @Override
        public WrappedVector drawInitialMomentum() {

            double[] momentum = new double[dim];

            for (int i = 0; i < dim; i++) {
                momentum[i] = MathUtils.nextGaussian() * Math.sqrt(1.0 / inverseMass[i]);
            }

            return new WrappedVector.Raw(momentum);
        }

        @Override
        protected double[] computeInverseMass() {

            double[] diagonalHessian = hessian.getDiagonalHessianLogDensity();

            // TODO Check transformation
            if (transform != null) {

                double[] untransformedValues = hessian.getParameter().getParameterValues();

                double[] gradient = hessian.getGradientLogDensity();

                diagonalHessian = transform.updateDiagonalHessianLogDensity(
                        diagonalHessian, gradient, untransformedValues, 0, dim

                );
            }

            return boundMassInverse(diagonalHessian);
        }
        
        private double[] boundMassInverse(double[] diagonalHessian) {

            double sum = 0.0;
            final double lowerBound = 1E-2; //TODO bad magic numbers
            final double upperBound = 1E2;
            double[] boundedMassInverse = new double[dim];

            for (int i = 0; i < dim; i++) {
                boundedMassInverse[i] = -1.0 / diagonalHessian[i];
                if (boundedMassInverse[i] < lowerBound) {
                    boundedMassInverse[i] = lowerBound;
                } else if (boundedMassInverse[i] > upperBound) {
                    boundedMassInverse[i] = upperBound;
                }
                sum += 1.0 / boundedMassInverse[i];
            }
            final double mean = sum / dim;
            for (int i = 0; i < dim; i++) {
                boundedMassInverse[i] = boundedMassInverse[i] * mean;
            }
            return boundedMassInverse;
        }

        @Override
        public double getVelocity(int i, ReadableVector momentum) {
            return momentum.get(i) * inverseMass[i];
        }
    }

    // TODO Implement
    class FullPreconditioning extends HessianBased {

        FullPreconditioning(HessianWrtParameterProvider hessian,
                            Transform transform) {
            super(hessian, transform);
        }

        @Override
        protected double[] computeInverseMass() {
            double[][] hessianMatrix = hessian.getHessianLogDensity();
            if (transform != null) {
                //TODO: hessianMatrix transform
            }

            Algebra algebra = new Algebra();

            DoubleMatrix2D H = new DenseDoubleMatrix2D(hessianMatrix);
            RobustEigenDecomposition decomposition = new RobustEigenDecomposition(H);
            DoubleMatrix1D eigenvalues = decomposition.getRealEigenvalues();

            normalizeEigenvalues(eigenvalues);

            DoubleMatrix2D V = decomposition.getV();

            inverseNegateEigenvalues(eigenvalues);

            double[][] negativeHessianInverse = algebra.mult(
                    algebra.mult(V, DoubleFactory2D.dense.diagonal(eigenvalues)),
                    algebra.inverse(V)
            ).toArray();

            double[] inverseMassArray = new double[hessian.getDimension() * hessian.getDimension()];
            for (int i = 0; i < hessian.getDimension(); i++) {
                System.arraycopy(negativeHessianInverse[i], 0, inverseMassArray, i * hessian.getDimension(), hessian.getDimension());
            }
            return inverseMassArray;
        }

        private static final double MIN_EIGENVALUE = -10.0; // TODO Bad magic number
        private static final double MAX_EIGENVALUE = -0.1; // TODO Bad magic number

        private void boundEigenvalues(DoubleMatrix1D eigenvalues) {

            for (int i = 0; i < eigenvalues.cardinality(); ++i) {
                if (eigenvalues.get(i) > MAX_EIGENVALUE) {
                    eigenvalues.set(i, MAX_EIGENVALUE);
                } else if (eigenvalues.get(i) < MIN_EIGENVALUE) {
                    eigenvalues.set(i, MIN_EIGENVALUE);
                }
            }
        }

        private void scaleEigenvalues(DoubleMatrix1D eigenvalues) {
            double sum = 0.0;
            for (int i = 0; i < eigenvalues.cardinality(); ++i) {
                sum += eigenvalues.get(i);
            }

            double mean = -sum / eigenvalues.cardinality();

            for (int i = 0; i < eigenvalues.cardinality(); ++i) {
                eigenvalues.set(i, eigenvalues.get(i) / mean);
            }
        }

        private void normalizeEigenvalues(DoubleMatrix1D eigenvalues) {
            boundEigenvalues(eigenvalues);
            scaleEigenvalues(eigenvalues);
        }

        private void inverseNegateEigenvalues(DoubleMatrix1D eigenvalues) {
            for (int i = 0; i < eigenvalues.cardinality(); i++) {
                eigenvalues.set(i, -1.0 / eigenvalues.get(i));
            }
        }

        @Override
        public WrappedVector drawInitialMomentum() {

            MultivariateNormalDistribution mvn = new MultivariateNormalDistribution(
                    new double[dim], toArray(inverseMass, dim, dim)
            );

            return new WrappedVector.Raw(mvn.nextMultivariateNormal());
        }

        @Override
        public double getVelocity(int i, ReadableVector momentum) {
            double velocity = 0.0;

            for (int j = 0; j < dim; ++j) {
                velocity += inverseMass[i * dim + j] * momentum.get(j);
            }

            return velocity;
        }

        private static double[][] toArray(double[] vector, int rowDim, int colDim) {
            double[][] array = new double[rowDim][];

            for (int row = 0; row < rowDim; ++row) {
                array[row] = new double[colDim];
                System.arraycopy(vector, colDim * row,
                        array[row], 0, colDim);
            }

            return array;
        }
    }

    abstract // TODO Implement interface
    class SecantPreconditioing extends HessianBased {

        private final Secant[] queue;
        private int secantIndex;
        private int secantUpdateCount;

        SecantPreconditioing(HessianWrtParameterProvider hessian, int secantSize) {
            super(hessian, null);

            this.queue = new Secant[secantSize];
            this.secantIndex = 0;
            this.secantUpdateCount = 0;
        }

        @Override
        public void storeSecant(ReadableVector gradient, ReadableVector position) {
            queue[secantIndex] = new Secant(gradient, position);
            secantIndex = (secantIndex + 1) % queue.length;
            ++secantUpdateCount;
        }

        private class Secant { // TODO Inner class because we may change the storage depending on efficiency

            ReadableVector gradient;
            ReadableVector position;

            Secant(ReadableVector gradient, ReadableVector position) {
                this.gradient = gradient;
                this.position = position;
            }
        }
    }
}
