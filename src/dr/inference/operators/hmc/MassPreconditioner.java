package dr.inference.operators.hmc;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.Parameter;
import dr.inference.model.PriorPreconditioningProvider;
import dr.math.*;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.RobustEigenDecomposition;
import dr.math.matrixAlgebra.WrappedMatrix;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Transform;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public interface MassPreconditioner {

    WrappedVector drawInitialMomentum();

    double getVelocity(int index, ReadableVector momentum);

    void storeSecant(ReadableVector gradient, ReadableVector position);

    void updateMass();

    WrappedVector getMass();

    void updateVariance(WrappedVector position);

    ReadableVector doCollision(int[] indices, ReadableVector momentum);

    int getDimension();

    enum Type {

        NONE("none") {
            @Override
            public MassPreconditioner factory(GradientWrtParameterProvider gradient, Transform transform, MassPreconditioningOptions options) {
                final Parameter parameter = gradient.getParameter();
                int dim = parameter.getDimension();
                if (transform != null && transform instanceof Transform.MultivariableTransform) {
                    dim = ((Transform.MultivariableTransform) transform).getDimension();
                }
                return new NoPreconditioning(dim);
            }
        },
        DIAGONAL("diagonal") {
            @Override
            public MassPreconditioner factory(GradientWrtParameterProvider gradient, Transform transform, MassPreconditioningOptions options) {
                return new DiagonalHessianPreconditioning((HessianWrtParameterProvider) gradient, transform, options.preconditioningMemory());
            }
        },
        ADAPTIVE_DIAGONAL("adaptiveDiagonal") {
            @Override
            public MassPreconditioner factory(GradientWrtParameterProvider gradient, Transform transform, MassPreconditioningOptions options) {
                int dimension = transform instanceof Transform.MultivariableTransform ?
                        ((Transform.MultivariableTransform) transform).getDimension() : gradient.getDimension();
                return new AdaptiveDiagonalPreconditioning(dimension, gradient, transform, options.preconditioningDelay());
            }
        },
        PRIOR_DIAGONAL("priorDiagonal") {
            @Override
            public MassPreconditioner factory(GradientWrtParameterProvider gradient, Transform transform, MassPreconditioningOptions options) {

                if (!(gradient instanceof PriorPreconditioningProvider)) {
                    throw new RuntimeException("Gradient must be a PriorPreconditioningProvider for prior preconditioning!");
                }

                return new PriorPreconditioner((PriorPreconditioningProvider) gradient, transform);

            }
        },
        FULL("full") {
            @Override
            public MassPreconditioner factory(GradientWrtParameterProvider gradient, Transform transform, MassPreconditioningOptions options) {
                return new FullHessianPreconditioning((HessianWrtParameterProvider) gradient, transform);
            }
        },
        SECANT("secant") {
            @Override
            public MassPreconditioner factory(GradientWrtParameterProvider gradient, Transform transform, MassPreconditioningOptions options) {
                SecantHessian secantHessian = new SecantHessian(gradient, options.preconditioningMemory());
                return new Secant(secantHessian, transform);
            }
        },
        ADAPTIVE("adaptive") {
            @Override
            public MassPreconditioner factory(GradientWrtParameterProvider gradient, Transform transform, MassPreconditioningOptions options) {
//                AdaptableCovariance adaptableCovariance = new AdaptableCovariance.WithSubsampling(gradient.getDimension(), 1000);
                AdaptableCovariance adaptableCovariance = new AdaptableCovariance(gradient.getDimension());
                return new AdaptiveFullHessianPreconditioning(gradient, adaptableCovariance, transform, gradient.getDimension(), options.preconditioningDelay());
            }
        };

        private final String name;

        Type(String name) {
            this.name = name;
        }

        public abstract MassPreconditioner factory(GradientWrtParameterProvider gradient, Transform transform, MassPreconditioningOptions options);

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

    class CompoundPreconditioning implements MassPreconditioner {

        final int dim;
        final List<MassPreconditioner> preconditionerList;
        boolean velocityKnown = false;
        double[] velocity;

        CompoundPreconditioning(List<MassPreconditioner> preconditionerList) {

            int thisDim = 0;
            for (MassPreconditioner preconditioner : preconditionerList) {
                thisDim += preconditioner.getDimension();
            }
            this.dim = thisDim;
            this.preconditionerList = preconditionerList;
            this.velocity = new double[dim];
        }

        @Override
        public WrappedVector drawInitialMomentum() {
            WrappedVector initialMomentum = new WrappedVector.Raw(new double[dim]);
            int currentIndex = 0;
            for (MassPreconditioner preconditioner : preconditionerList) {
                WrappedVector currentMomentum = preconditioner.drawInitialMomentum();
                for (int i = 0; i < preconditioner.getDimension(); i++) {
                    initialMomentum.set(currentIndex + i, currentMomentum.get(i));
                }
                currentIndex += preconditioner.getDimension();
            }
            return initialMomentum;
        }

        @Override
        public double getVelocity(int index, ReadableVector momentum) {
            getVelocityVector(momentum);
            return velocity[index];
        }

        private void getVelocityVector(ReadableVector momentum) {
            if (!velocityKnown) {
                int currentIndex = 0;
                List<ReadableVector> separatedMomentum = separateVectors(momentum);
                for (int j = 0; j < preconditionerList.size(); j++) {
                    MassPreconditioner preconditioner = preconditionerList.get(j);
                    ReadableVector currentMomentum = separatedMomentum.get(j);
                    for (int i = 0; i < preconditioner.getDimension(); i++) {
                        velocity[currentIndex + i] = preconditioner.getVelocity(i, currentMomentum);
                    }
                    currentIndex += preconditioner.getDimension();
                }
                velocityKnown = true;
            }
        }

        private List<ReadableVector> separateVectors(ReadableVector rawVector) {
            List<ReadableVector> vectors = new ArrayList<>();
            int currentIndex = 0;
            for (MassPreconditioner preconditioner : preconditionerList) {
                WrappedVector thisVector = new WrappedVector.Raw(new double[preconditioner.getDimension()]);
                for (int i = 0; i < preconditioner.getDimension(); i++) {
                    thisVector.set(i, rawVector.get(currentIndex + i));
                }
                vectors.add(thisVector);
                currentIndex += preconditioner.getDimension();
            }
            return vectors;
        }

        private ReadableVector combineVectors(List<ReadableVector> vectors) {
            WrappedVector combinedVector = new WrappedVector.Raw(new double[dim]);
            int currentIndex = 0;
            for (ReadableVector readableVector : vectors) {
                for (int i = 0; i < readableVector.getDim(); i++) {
                    combinedVector.set(currentIndex + i, readableVector.get(i));
                }
            }
            return combinedVector;
        }

        @Override
        public void storeSecant(ReadableVector gradient, ReadableVector position) {
            List<ReadableVector> separatedGradient = separateVectors(gradient);
            List<ReadableVector> separatedPosition = separateVectors(position);
            for (int i = 0; i < preconditionerList.size(); i++) {
                preconditionerList.get(i).storeSecant(separatedGradient.get(i), separatedPosition.get(i));
            }
        }

        @Override
        public void updateMass() {
            for (MassPreconditioner preconditioner : preconditionerList) {
                preconditioner.updateMass();
            }
            velocityKnown = false;
        }

        @Override
        public ReadableVector doCollision(int[] indices, ReadableVector momentum) {
            throw new RuntimeException("Not yet implemented!");
        }

        @Override
        public WrappedVector getMass() {
            throw new RuntimeException("Not yet implemented!");
        }

        @Override
        public void updateVariance(WrappedVector position) {
            // Do nothing
        }

        @Override
        public int getDimension() {
            return dim;
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

        @Override
        public ReadableVector doCollision(int[] indices, ReadableVector momentum) {
            if (indices.length != 2) {
                throw new RuntimeException("Not implemented for more than two dimensions yet.");
            }
            WrappedVector updatedMomentum = new WrappedVector.Raw(new double[momentum.getDim()]);

            for (int i = 0; i < momentum.getDim(); i++) {
                updatedMomentum.set(i, momentum.get(i));
            }

            updatedMomentum.set(indices[0], momentum.get(indices[1]));
            updatedMomentum.set(indices[1], momentum.get(indices[0]));
            return updatedMomentum;
        }

        @Override
        public WrappedVector getMass() {
            double[] mass = new double[dim];
            Arrays.fill(mass, 1.0);
            return new WrappedVector.Raw(mass);
        }

        @Override
        public void updateVariance(WrappedVector position) {
            // Do nothing
        }

        @Override
        public int getDimension() {
            return dim;
        }
    }



    abstract class AbstractMassPreconditioning implements MassPreconditioner {
        final protected int dim;
        final protected Transform transform;
        double[] inverseMass;

        protected AbstractMassPreconditioning(int dim, Transform transform) {
            this.dim = dim;
            this.transform = transform;
        }

        abstract protected void initializeMass();

        abstract protected double[] computeInverseMass();

        public void updateMass() {
            this.inverseMass = computeInverseMass();
        }

        @Override
        public int getDimension() {
            return dim;
        }

        abstract public void storeSecant(ReadableVector gradient, ReadableVector position);

    }

    abstract class HessianBased extends AbstractMassPreconditioning {

        final protected HessianWrtParameterProvider hessian;

        // TODO Should probably make a TransformedHessian so that this class does not need to know about transformations
        HessianBased(HessianWrtParameterProvider hessian,
                     Transform transform) {

            this(hessian, transform, hessian.getDimension());

        }

        HessianBased(HessianWrtParameterProvider hessian, Transform transform, int dim) {

            super(dim, transform);
            this.hessian = hessian;

            initializeMass();
        }

        @Override
        public ReadableVector doCollision(int[] indices, ReadableVector momentum) {
            throw new RuntimeException("Not yet implemented.");
        }
    }

    abstract class DiagonalPreconditioning extends AbstractMassPreconditioning {

        protected AdaptableVector adaptiveDiagonal;

        protected DiagonalPreconditioning(int dim, Transform transform) {
            super(dim, transform);
            this.adaptiveDiagonal = new AdaptableVector.Default(dim);
            initializeMass();
        }

        @Override
        protected void initializeMass() {

            double[] result = new double[dim];
            Arrays.fill(result, 1.0);
            inverseMass = normalizeVector(new WrappedVector.Raw(result), dim);
        }

        protected double[] normalizeVector(ReadableVector values, double targetSum) {
            double sum = 0.0;
            for (int i = 0; i < values.getDim(); i++) {
                sum += values.get(i);
            }
            final double multiplier = targetSum / sum;
            double[] normalizedValues = new double[values.getDim()];
            for (int i = 0; i < values.getDim(); i++) {
                normalizedValues[i] = values.get(i) * multiplier;
            }
            return normalizedValues;
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
        public double getVelocity(int i, ReadableVector momentum) {
            return momentum.get(i) * inverseMass[i];
        }

        @Override
        public ReadableVector doCollision(int[] indices, ReadableVector momentum) {
            if (indices.length != 2) {
                throw new RuntimeException("Not implemented for more than two dimensions yet.");
            }

            WrappedVector updatedMomentum = new WrappedVector.Raw(new double[momentum.getDim()]);

            for (int i = 0; i < momentum.getDim(); i++) {
                updatedMomentum.set(i, momentum.get(i));
            }

            final int index1 = indices[0];
            final int index2 = indices[1];
            final double updatedMomentum1 = ((inverseMass[index2] - inverseMass[index1]) * momentum.get(index1)
                    + 2.0 * inverseMass[index2] * momentum.get(index2)) / (inverseMass[index1] + inverseMass[index2]);
            final double updatedMomentum2 = ((inverseMass[index1] - inverseMass[index2]) * momentum.get(index2)
                    + 2.0 * inverseMass[index1] * momentum.get(index1)) / (inverseMass[index1] + inverseMass[index2]);

            updatedMomentum.set(index1, updatedMomentum1);
            updatedMomentum.set(index2, updatedMomentum2);
            return updatedMomentum;
        }
    }

    class PriorPreconditioner extends DiagonalPreconditioning {

        PriorPreconditioningProvider priorDistribution;

        public PriorPreconditioner(PriorPreconditioningProvider priorDistribution, Transform transform){
            super(priorDistribution.getDimension(), transform);
            this.priorDistribution = priorDistribution;
            inverseMass = computeInverseMass();
        }

        public MassPreconditioner factory(PriorPreconditioningProvider priorDistribution, Transform transform) {
            return new PriorPreconditioner(priorDistribution, transform);
        }

        @Override
        protected double[] computeInverseMass() {
            double diagonal[] = new double[(priorDistribution.getDimension())];

            for (int i = 0; i < priorDistribution.getDimension(); i++){
                double stDev = priorDistribution.getStandardDeviation(i);
                diagonal[i] = stDev * stDev;
            }
            return diagonal;
        }

        @Override
        public void storeSecant(ReadableVector gradient, ReadableVector position) {
            // Do nothing
        }

        @Override
        public void updateVariance(WrappedVector position) {
            // Do nothing
        }

        @Override
        public WrappedVector getMass() {
            throw new RuntimeException("Not yet implemented!");
        }
    }

    class DiagonalHessianPreconditioning extends DiagonalPreconditioning {

        final protected HessianWrtParameterProvider hessian;

        DiagonalHessianPreconditioning(HessianWrtParameterProvider hessian,
                                       Transform transform,
                                       int memorySize) {
            super(hessian.getDimension(), transform);
            this.hessian = hessian;
            if (memorySize > 0) {
                this.adaptiveDiagonal = new AdaptableVector.LimitedMemory(hessian.getDimension(), memorySize);
            } else {
                this.adaptiveDiagonal = new AdaptableVector.Default(hessian.getDimension());
            }
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

            adaptiveDiagonal.update(new WrappedVector.Raw(newDiagonalHessian));

            return boundMassInverse(((WrappedVector) adaptiveDiagonal.getMean()).getBuffer());
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
        public void storeSecant(ReadableVector gradient, ReadableVector position) {
            // Do nothing
        }

        @Override
        public void updateVariance(WrappedVector position) {

        }

        @Override
        public WrappedVector getMass() {
            throw new RuntimeException("Not yet implemented!");
        }

    }


    class AdaptiveDiagonalPreconditioning extends DiagonalPreconditioning {

        private AdaptableVector.AdaptableVariance variance;
        private final int minimumUpdates;
        private final GradientWrtParameterProvider gradient;

        AdaptiveDiagonalPreconditioning(int dim,
                                        GradientWrtParameterProvider gradient,
                                        Transform transform, int preconditioningDelay) {
            super(dim, transform);
            this.variance = new AdaptableVector.AdaptableVariance(dim);
            this.minimumUpdates = preconditioningDelay;
            this.gradient = gradient;
        }

        private void setInitialMass() {
            double[] values = gradient.getParameter().getParameterValues();
            double[] storedValues = values.clone();
            for (int i = 0; i < dim; i++) {
                gradient.getParameter().setParameterValueQuietly(i, values[i] + MachineAccuracy.SQRT_SQRT_EPSILON);
            }
            gradient.getParameter().fireParameterChangedEvent();
            double[] gradientPlus = gradient.getGradientLogDensity();

            for (int i = 0; i < dim; i++) {
                gradient.getParameter().setParameterValueQuietly(i, values[i] - MachineAccuracy.SQRT_SQRT_EPSILON);
            }
            gradient.getParameter().fireParameterChangedEvent();
            double[] gradientMinus = gradient.getGradientLogDensity();

            for (int i = 0; i < dim; i++) {
                gradient.getParameter().setParameterValueQuietly(i, values[i]);
            }
            gradient.getParameter().fireParameterChangedEvent();


            for (int i = 0; i < dim; i++) {
                values[i] = Math.abs((gradientPlus[i] - gradientMinus[i]) / (2.0 * MachineAccuracy.SQRT_SQRT_EPSILON));
                gradient.getParameter().setParameterValueQuietly(i, storedValues[i]);
            }
            gradient.getParameter().fireParameterChangedEvent();
            inverseMass = normalizeVector(new WrappedVector.Raw(values), dim);
        }

        @Override
        protected double[] computeInverseMass() {

            if (variance.getUpdateCount() > minimumUpdates) {
                double[] newVariance = variance.getVariance();
//                adaptiveDiagonal.update(new WrappedVector.Raw(newVariance));
//                return normalizeVector(adaptiveDiagonal.getMean(), dim);
                return normalizeVector(new WrappedVector.Raw(newVariance), dim);
            } else {
                return inverseMass;
            }

        }

        @Override
        public void storeSecant(ReadableVector gradient, ReadableVector position) {
             variance.update(position);
        }

        @Override
        public void updateVariance(WrappedVector position) {
            variance.update(position);
        }

        @Override
        public WrappedVector getMass() {
            double[] mass = new double[dim];
            for (int i = 0; i < dim; i++) {
                mass[i] = 1 / inverseMass[i];
            }
            return new WrappedVector.Raw(mass);
        }
    }

    class FullHessianPreconditioning extends HessianBased {

        FullHessianPreconditioning(HessianWrtParameterProvider hessian,
                                   Transform transform) {
            super(hessian, transform);
        }

        FullHessianPreconditioning(HessianWrtParameterProvider hessian, Transform transform, int dim) {
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

            private static final double MIN_EIGENVALUE = -20.0; // TODO Bad magic number
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
        public void storeSecant(ReadableVector gradient, ReadableVector position) {
            // Do nothing
        }

        @Override
        public void updateVariance(WrappedVector position) {

        }

        @Override
        public WrappedVector getMass() {
            throw new RuntimeException("Not yet implemented!");
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

    class Secant extends FullHessianPreconditioning {

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

    class AdaptiveFullHessianPreconditioning extends FullHessianPreconditioning {

        private final AdaptableCovariance adaptableCovariance;
        private final GradientWrtParameterProvider gradientProvider;
        private final AdaptableVector averageCovariance;
        private final double[] inverseMassBuffer;
        private final int minimumUpdates;

        AdaptiveFullHessianPreconditioning(GradientWrtParameterProvider gradientProvider,
                                           AdaptableCovariance adaptableCovariance,
                                           Transform transform,
                                           int dim,
                                           int preconditioningDelay) {
            super(null, transform, dim);
            this.adaptableCovariance = adaptableCovariance;
            this.gradientProvider = gradientProvider;
            this.averageCovariance = new AdaptableVector.Default(dim * dim);
            this.inverseMassBuffer = new double[dim * dim];
            this.minimumUpdates = preconditioningDelay;
        }

        @Override
        protected double[] computeInverseMass() {

            if (adaptableCovariance.getUpdateCount() > minimumUpdates) {
                WrappedMatrix.ArrayOfArray covariance = (WrappedMatrix.ArrayOfArray) adaptableCovariance.getCovariance();

                double[] flatCovariance = new double[dim * dim];
                for (int i = 0; i < dim; i++) {
                    System.arraycopy(covariance.getArrays()[i], 0, flatCovariance, i * dim, dim);
                }

                averageCovariance.update(new WrappedVector.Raw(flatCovariance));

//            double[][] numericHessian = NumericalDerivative.getNumericalHessian(numeric1, gradientProvider.getParameter().getParameterValues());

                cacheAverageCovariance(normalizeCovariance((WrappedVector.Raw) averageCovariance.getMean()));

                return inverseMassBuffer;
            } else {
                return inverseMass;
            }
        }

        private ReadableVector normalizeCovariance(WrappedVector flatCovariance) {
            double sum = 0.0;
            for (int i = 0; i < dim; i++) {
                sum += flatCovariance.get(i * dim + i);
            }

            final double multiplier = dim / sum;
            for (int i = 0; i < dim * dim; i++) {
                flatCovariance.set(i, flatCovariance.get(i) * multiplier);
            }
            return flatCovariance;
        }

        private void cacheAverageCovariance(ReadableVector mean) {

            double[][] tempVariance = new double[dim][dim];
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    tempVariance[i][j] = -mean.get(i * dim + j);
                }
            }

            double[] transformedVariance = PDTransformMatrix.Default.transformMatrix(tempVariance, dim);


            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    inverseMassBuffer[i * dim + j] = transformedVariance[i * dim + j];
                }
            }
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
