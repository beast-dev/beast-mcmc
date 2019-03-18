package dr.inference.operators.hmc;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.math.*;
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
            public MassPreconditioner factory(GradientWrtParameterProvider gradient, Transform transform, HamiltonianMonteCarloOperator.Options options) {
                return new NoPreconditioning(gradient.getDimension());
            }
        },
        DIAGONAL("diagonal") {
            @Override
            public MassPreconditioner factory(GradientWrtParameterProvider gradient, Transform transform, HamiltonianMonteCarloOperator.Options options) {
                return new DiagonalPreconditioning((HessianWrtParameterProvider) gradient, transform, options.preconditioningMemory);
            }
        },
        FULL("full") {
            @Override
            public MassPreconditioner factory(GradientWrtParameterProvider gradient, Transform transform, HamiltonianMonteCarloOperator.Options options) {
                return new FullPreconditioning((HessianWrtParameterProvider) gradient, transform);
            }
        },
        SECANT("secant") {
            @Override
            public MassPreconditioner factory(GradientWrtParameterProvider gradient, Transform transform, HamiltonianMonteCarloOperator.Options options) {
                SecantHessian secantHessian = new SecantHessian(gradient, options.preconditioningMemory);
                return new Secant(secantHessian, transform);
            }
        },
        ADAPTIVE("adaptive") {
            @Override
            public MassPreconditioner factory(GradientWrtParameterProvider gradient, Transform transform, HamiltonianMonteCarloOperator.Options options) {
//                AdaptableCovariance adaptableCovariance = new AdaptableCovariance.WithSubsampling(gradient.getDimension(), 1000);
                AdaptableCovariance adaptableCovariance = new AdaptableCovariance(gradient.getDimension());
                return new AdaptivePreconditioning(gradient, adaptableCovariance, transform, gradient.getDimension());
            }
        };

        private final String name;

        Type(String name) {
            this.name = name;
        }

        public abstract MassPreconditioner factory(GradientWrtParameterProvider gradient, Transform transform, HamiltonianMonteCarloOperator.Options options);

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

            initializeMass();
        }

        public void storeSecant(ReadableVector gradient, ReadableVector position) {
            // Do nothing
        }

        abstract protected void initializeMass();

        abstract protected double[] computeInverseMass();

        public void updateMass() {
            this.inverseMass = computeInverseMass();
        }
    }

    class DiagonalPreconditioning extends HessianBased {

        private AdaptableVector adaptiveDiagonalHessian;

        DiagonalPreconditioning(HessianWrtParameterProvider hessian,
                                Transform transform) {
            this(hessian, transform, 0);
        }

        DiagonalPreconditioning(HessianWrtParameterProvider hessian,
                                Transform transform,
                                int memorySize) {
            super(hessian, transform);
            if (memorySize > 0) {
                this.adaptiveDiagonalHessian = new AdaptableVector.LimitedMemory(hessian.getDimension(), memorySize);
            } else {
                this.adaptiveDiagonalHessian = new AdaptableVector.Default(hessian.getDimension());
            }
        }



        @Override
        protected void initializeMass() {

            double[] result = new double[dim];
            for (int i = 0; i < dim; i++) {
                result[i] = 1.0;
            }
            inverseMass = result;
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

            double[] newDiagonalHessian = hessian.getDiagonalHessianLogDensity();

            if (transform != null) {

                double[] untransformedValues = hessian.getParameter().getParameterValues();

                double[] gradient = hessian.getGradientLogDensity();

                newDiagonalHessian = transform.updateDiagonalHessianLogDensity(
                        newDiagonalHessian, gradient, untransformedValues, 0, dim

                );
            }

            adaptiveDiagonalHessian.update(new WrappedVector.Raw(newDiagonalHessian));

            return boundMassInverse(((WrappedVector) adaptiveDiagonalHessian.getMean()).getBuffer());
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

        @Override
        protected void initializeMass() {

            double[] result = new double[dim * dim];
            for (int i = 0; i < dim; i++) {
                result[i * dim + i] = 1.0;
            }
            inverseMass = result;

        }

        enum PDTransformMatrix{
            Invert("Transform inverse matrix into a PD matrix") {
                @Override
                protected void transformEigenvalues(DoubleMatrix1D eigenvalues) {
                    inverseNegateEigenvalues(eigenvalues);
                }
            },
            Default("Transform matrix into a PD matrix") {
                @Override
                protected void transformEigenvalues(DoubleMatrix1D eigenvalues) {
                    negateEigenvalues(eigenvalues);
                }
            },
            Negate("Transform negative matrix into a PD matrix") {
                @Override
                protected void transformEigenvalues(DoubleMatrix1D eigenvalues) {
                    negateEigenvalues(eigenvalues);
                }
                @Override
                protected void normalizeEigenvalues(DoubleMatrix1D eigenvalues) {
                    negateEigenvalues(eigenvalues);
                    boundEigenvalues(eigenvalues);
                    scaleEigenvalues(eigenvalues);
                }
            },
            NegateInvert("Transform negative inverse matrix into a PD matrix") {
                @Override
                protected void transformEigenvalues(DoubleMatrix1D eigenvalues) {
                    inverseNegateEigenvalues(eigenvalues);
                }
                @Override
                protected void normalizeEigenvalues(DoubleMatrix1D eigenvalues) {
                    negateEigenvalues(eigenvalues);
                    boundEigenvalues(eigenvalues);
                    scaleEigenvalues(eigenvalues);
                }
            };

            String desc;

            PDTransformMatrix(String s) {
                desc = s;
            }

            public String toString() {
                return desc;
            }

            private static final double MIN_EIGENVALUE = -10.0; // TODO Bad magic number
            private static final double MAX_EIGENVALUE = -0.5; // TODO Bad magic number

            protected void boundEigenvalues(DoubleMatrix1D eigenvalues) {

                for (int i = 0; i < eigenvalues.cardinality(); ++i) {
                    if (eigenvalues.get(i) > MAX_EIGENVALUE) {
                        eigenvalues.set(i, MAX_EIGENVALUE);
                    } else if (eigenvalues.get(i) < MIN_EIGENVALUE) {
                        eigenvalues.set(i, MIN_EIGENVALUE);
                    }
                }
            }

            protected void scaleEigenvalues(DoubleMatrix1D eigenvalues) {
                double sum = 0.0;
                for (int i = 0; i < eigenvalues.cardinality(); ++i) {
                    sum += eigenvalues.get(i);
                }

                double mean = -sum / eigenvalues.cardinality();

                for (int i = 0; i < eigenvalues.cardinality(); ++i) {
                    eigenvalues.set(i, eigenvalues.get(i) / mean);
                }
            }

            protected void normalizeEigenvalues(DoubleMatrix1D eigenvalues) {
                boundEigenvalues(eigenvalues);
                scaleEigenvalues(eigenvalues);
            }

            protected void inverseNegateEigenvalues(DoubleMatrix1D eigenvalues) {
                for (int i = 0; i < eigenvalues.cardinality(); i++) {
                    eigenvalues.set(i, -1.0 / eigenvalues.get(i));
                }
            }

            protected void negateEigenvalues(DoubleMatrix1D eigenvalues) {
                for (int i = 0; i < eigenvalues.cardinality(); i++) {
                    eigenvalues.set(i, -eigenvalues.get(i));
                }
            }

            public double[] transformMatrix(double[][] inputMatrix, int dim) {
                Algebra algebra = new Algebra();

                DoubleMatrix2D H = new DenseDoubleMatrix2D(inputMatrix);
                RobustEigenDecomposition decomposition = new RobustEigenDecomposition(H);
                DoubleMatrix1D eigenvalues = decomposition.getRealEigenvalues();

                normalizeEigenvalues(eigenvalues);

                DoubleMatrix2D V = decomposition.getV();

                transformEigenvalues(eigenvalues);

                double[][] negativeHessianInverse = algebra.mult(
                        algebra.mult(V, DoubleFactory2D.dense.diagonal(eigenvalues)),
                        algebra.inverse(V)
                ).toArray();

                double[] massArray = new double[dim * dim];
                for (int i = 0; i < dim; i++) {
                    System.arraycopy(negativeHessianInverse[i], 0, massArray, i * dim, dim);
                }
                return massArray;
            }

            abstract protected  void transformEigenvalues(DoubleMatrix1D eigenvalues);
        }

        private double[] computeInverseMass(WrappedMatrix.ArrayOfArray hessianMatrix, GradientWrtParameterProvider gradientProvider, PDTransformMatrix pdTransformMatrix) {

            double[][] transformedHessian = hessianMatrix.getArrays();

            if (transform != null) { // TODO: change 0 matrix into general TransformationHessian
                transformedHessian = transform.updateHessianLogDensity(transformedHessian, new double[dim][dim], gradientProvider.getGradientLogDensity(), gradientProvider.getParameter().getParameterValues(), 0, dim);
            }

            return pdTransformMatrix.transformMatrix(transformedHessian, dim);
        }

        @Override
        protected double[] computeInverseMass() {
            //TODO: change to ReadableMatrix
            WrappedMatrix.ArrayOfArray hessianMatrix = new WrappedMatrix.ArrayOfArray(hessian.getHessianLogDensity());

            return computeInverseMass(hessianMatrix, hessian, PDTransformMatrix.Invert);
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
        private final GradientWrtParameterProvider gradientProvider;

        AdaptivePreconditioning(GradientWrtParameterProvider gradientProvider, AdaptableCovariance adaptableCovariance, Transform transform, int dim) {
            super(null, transform, dim);
            this.adaptableCovariance = adaptableCovariance;
            this.gradientProvider = gradientProvider;
        }

        @Override
        protected double[] computeInverseMass() {

            WrappedMatrix.ArrayOfArray covariance = (WrappedMatrix.ArrayOfArray) adaptableCovariance.getCovariance();

            double[][] numericHessian = NumericalDerivative.getNumericalHessian(numeric1, gradientProvider.getParameter().getParameterValues());

            return super.computeInverseMass(covariance, gradientProvider, PDTransformMatrix.Negate);

        }

        @Override
        public void storeSecant(ReadableVector gradient, ReadableVector position) {
//            adaptableCovariance.update(gradient);
            adaptableCovariance.update(position);
        }

        protected MultivariateFunction numeric1 = new MultivariateFunction() {
            @Override
            public double evaluate(double[] argument) {

                for (int i = 0; i < argument.length; ++i) {
                    gradientProvider.getParameter().setParameterValue(i, argument[i]);
                }

//            treeDataLikelihood.makeDirty();
                return gradientProvider.getLikelihood().getLogLikelihood();
            }

            @Override
            public int getNumArguments() {
                return gradientProvider.getParameter().getDimension();
            }

            @Override
            public double getLowerBound(int n) {
                return 0;
            }

            @Override
            public double getUpperBound(int n) {
                return Double.POSITIVE_INFINITY;
            }
        };

    }
}
