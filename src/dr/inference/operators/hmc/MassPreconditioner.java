package dr.inference.operators.hmc;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.math.AdaptableCovariance;
import dr.math.MathUtils;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.*;
import dr.util.Transform;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public interface MassPreconditioner {

    WrappedVector drawInitialMomentum();

    double getVelocity(int index, ReadableVector momentum);

    void storeSecant(ReadableVector gradient, ReadableVector position);

    void updateMass();

    enum Type {

        NONE("none") {
            @Override
            public MassPreconditioner factory(GradientWrtParameterProvider gradient, Transform transform) {
                return new NoPreconditioning(gradient.getDimension());
            }
        },
        DIAGONAL("diagonal") {
            @Override
            public MassPreconditioner factory(GradientWrtParameterProvider gradient, Transform transform) {
                return new DiagonalPreconditioning((HessianWrtParameterProvider) gradient, transform);
            }
        },
        FULL("full") {
            @Override
            public MassPreconditioner factory(GradientWrtParameterProvider gradient, Transform transform) {
                return new FullPreconditioning((HessianWrtParameterProvider) gradient, transform);
            }
        },
        SECANT("secant") {
            @Override
            public MassPreconditioner factory(GradientWrtParameterProvider gradient, Transform transform) {
                SecantHessian secantHessian = new SecantHessian(gradient, 3);  // TODO make size an option
                return new Secant(secantHessian, transform);
            }
        },
        ADAPTIVE("adaptive") {
            @Override
            public MassPreconditioner factory(GradientWrtParameterProvider gradient, Transform transform) {
                AdaptableCovariance adaptableCovariance = new AdaptableCovariance(gradient.getDimension());
                return new AdaptivePreconditioning(adaptableCovariance, transform, gradient.getDimension());
            }
        };

        private final String name;

        Type(String name) {
            this.name = name;
        }

        public abstract MassPreconditioner factory(GradientWrtParameterProvider gradient, Transform transform);

        public String getName() { return name; }

        public static Type parseFromString(String text) {
            for (Type type : Type.values()) {
                if (type.name.toLowerCase().compareToIgnoreCase(text) == 0) {
                    return type;
                }
            }
            return Type.NONE;
        }
    }


    class NoPreconditioning implements MassPreconditioner {

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

        @Override
        public void updateMass() {
            // Do nothing
        }
    }

    abstract class HessianBased implements MassPreconditioner {
        final protected int dim;
        final protected HessianWrtParameterProvider hessian;
        final protected Transform transform;
        double[] inverseMass;

        // TODO Should probably make a TransformedHessian so that this class does not need to know about transformations
        HessianBased(HessianWrtParameterProvider hessian,
                     Transform transform) {

            this(hessian, transform, hessian.getDimension());

        }

        HessianBased(HessianWrtParameterProvider hessian, Transform transform, int dim) {
            this.dim = dim;
            this.hessian = hessian;
            this.transform = transform;

            updateMass();
        }

        public void storeSecant(ReadableVector gradient, ReadableVector position) {
            // Do nothing
        }

        abstract protected double[] computeInverseMass();

        public void updateMass() {
            this.inverseMass = computeInverseMass();
        }
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

    class FullPreconditioning extends HessianBased {

        FullPreconditioning(HessianWrtParameterProvider hessian,
                            Transform transform) {
            super(hessian, transform);
        }

        FullPreconditioning(HessianWrtParameterProvider hessian, Transform transform, int dim) {
            super(hessian, transform, dim);
        }

        private double[] computeInverseMass(WrappedMatrix.ArrayOfArray hessianMatrix) {

            double[][] transformedHessian = hessianMatrix.getArrays();

            if (transform != null) {
                transformedHessian = transform.updateHessianLogDensity(transformedHessian, new double[dim][dim], hessian.getGradientLogDensity(), hessian.getParameter().getParameterValues(), 0, dim);
            }
            Algebra algebra = new Algebra();

            DoubleMatrix2D H = new DenseDoubleMatrix2D(transformedHessian);
            RobustEigenDecomposition decomposition = new RobustEigenDecomposition(H);
            DoubleMatrix1D eigenvalues = decomposition.getRealEigenvalues();

            normalizeEigenvalues(eigenvalues);

            DoubleMatrix2D V = decomposition.getV();

            inverseNegateEigenvalues(eigenvalues);

            double[][] negativeHessianInverse = algebra.mult(
                    algebra.mult(V, DoubleFactory2D.dense.diagonal(eigenvalues)),
                    algebra.inverse(V)
            ).toArray();

            double[] inverseMassArray = new double[dim * dim];
            for (int i = 0; i < dim; i++) {
                System.arraycopy(negativeHessianInverse[i], 0, inverseMassArray, i * dim, dim);
            }
            return inverseMassArray;
        }

        @Override
        protected double[] computeInverseMass() {
            //TODO: change to ReadableMatrix
            WrappedMatrix.ArrayOfArray hessianMatrix = new WrappedMatrix.ArrayOfArray(hessian.getHessianLogDensity());

            return computeInverseMass(hessianMatrix);
        }

        private static final double MIN_EIGENVALUE = -10.0; // TODO Bad magic number
        private static final double MAX_EIGENVALUE = -0.5; // TODO Bad magic number

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

    class Secant extends FullPreconditioning {

        private final SecantHessian secantHessian;

        Secant(SecantHessian secantHessian, Transform transform) {
            super(secantHessian, transform);
            this.secantHessian = secantHessian;
        }

        @Override
        public void storeSecant(ReadableVector gradient, ReadableVector position) {
            secantHessian.storeSecant(gradient, position);
        }
    }

    class AdaptivePreconditioning extends FullPreconditioning {

        private final AdaptableCovariance adaptableCovariance;

        AdaptivePreconditioning(AdaptableCovariance adaptableCovariance, Transform transform, int dim) {
            super(null, transform, dim);
            this.adaptableCovariance = adaptableCovariance;
        }

        @Override
        protected double[] computeInverseMass() {

            if (adaptableCovariance != null) {

                WrappedMatrix.ArrayOfArray covariance = (WrappedMatrix.ArrayOfArray) adaptableCovariance.getCovariance();

                return super.computeInverseMass(covariance);
            } else {
                double[] output = new double[dim * dim];
                for (int i = 0; i < dim; i++) {
                    output[i * dim + i] = 1.0;
                }
                return output;
            }

        }

        @Override
        public void storeSecant(ReadableVector gradient, ReadableVector position) {
            adaptableCovariance.update(gradient);
        }
    }
}
