package dr.inference.operators.hmc;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.math.MachineAccuracy;
import dr.math.MultivariateFunction;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.RobustEigenDecomposition;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.util.Transform;

import java.util.Arrays;

/**
 * @author Marc A. Suchard
 */
interface MassProvider {

    double[][] getMass();

    double[][] getMassInverse();

    void updateMass();

    class Default implements MassProvider {
        final double drawVariance;
        final int dim;
        final double[][] mass;
        final double[][] massInverse;
        MultivariateNormalDistribution drawDistribution;

        Default(int dim, double drawVariance) {
            this.dim = dim;
            this.drawVariance = drawVariance;
            this.drawDistribution = setDrawDistribution(drawVariance);

            this.mass = new double[dim][dim];
            this.massInverse = new double[dim][dim];
            for (int i = 0; i < dim; i++) {
                Arrays.fill(mass[i], 0.0);
                Arrays.fill(massInverse[i], 0.0);
                massInverse[i][i] = drawVariance;
                mass[i][i] = 1.0 / massInverse[i][i];
            }
        }

        private MultivariateNormalDistribution setDrawDistribution(double drawVariance) {
            double[] mean = new double[dim];
            Arrays.fill(mean, 0.0);
            return new MultivariateNormalDistribution(mean, 1.0 / drawVariance);
        }

        public double[][] getMass() {
            return mass;
        }

        @Override
        public double[][] getMassInverse() {
            return massInverse;
        }

        @Override
        public void updateMass() {
            // Do nothing;
        }
    }

    class PreConditioning extends Default {

        final HessianWrtParameterProvider hessianWrtParameterProvider;

        PreConditioning(double drawVariance, HessianWrtParameterProvider hessianWrtParameterProvider) {
            super(hessianWrtParameterProvider.getDimension(), drawVariance);
            if (!(hessianWrtParameterProvider instanceof HessianWrtParameterProvider)) {
                throw new IllegalArgumentException("Must provide a HessianProvider for preConditioning.");
            }
            this.hessianWrtParameterProvider = hessianWrtParameterProvider;
//                setMassMatrices(hessianWrtParameterProvider.getDiagonalHessianLogDensity());
        }

        public void setMassMatrices(double[] diagonalHessian) {
            for (int i = 0; i < dim; i++) {
                Arrays.fill(massInverse[i], 0.0);
            }
            double[] boundedMassInverse = boundMassInverse(diagonalHessian);
            for (int i = 0; i < dim; i++) {
                massInverse[i][i] = boundedMassInverse[i];
                mass[i][i] = 1.0 / massInverse[i][i];
            }
        }

        @Override
        public void updateMass() {
            setMassMatrices(hessianWrtParameterProvider.getDiagonalHessianLogDensity());
        }

        private double[] boundMassInverse(double[] diagonalHessian) {

            double sum = 0.0;
            final double lowerBound = 1E-2;
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

    }

    class PreConditioningWithTransform extends PreConditioning {

        final Transform transform;

        PreConditioningWithTransform(double drawVariance, HessianWrtParameterProvider hessianWrtParameterProvider,
                                     Transform transform) {
            super(drawVariance, hessianWrtParameterProvider);
            this.transform = transform;
        }

        private void setArbitraryMatrix(double[][] matrix) {
            double[] diagonal = drawDistribution.nextMultivariateNormal();
            double multiplier = 1.2;
            for (int i = 0; i < dim; i++) {
                matrix[i][i] = Math.abs(diagonal[i]) * multiplier;
            }
        }

        public void setMassMatrices(double[] diagonalHessian) {
//                double[] gradient = hessianWrtParameterProvider.getGradientLogDensity();
//                double[] unTransformedPosition = hessianWrtParameterProvider.getParameter().getParameterValues();
//                diagonalHessian = transform.updateDiagonalHessianLogDensity(diagonalHessian, gradient, unTransformedPosition,
//                        0, diagonalHessian.length);
//
//                double[] testHessian = NumericalDerivative.diagonalHessian(numeric1, transform.transform(hessianWrtParameterProvider.getParameter().getParameterValues(), 0, dim));
            SymmetricMatrix hessian = getNumericalHessianCentral();
            double[][] hessianInverse = hessian.inverse().toComponents();


//                super.setMassMatrices(diagonalHessian);
//                setArbitraryMatrix(massMatrixInverse);
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    massInverse[i][j] = -hessianInverse[i][j];
                    mass[i][j] = -hessian.component(i, j);
                }
            }
        }

        private static final double MIN_EIGENVALUE = -0.5; // TODO Bad magic number

        private void boundEigenvalues(DoubleMatrix1D eigenvalues) {

            for (int i = 0; i < eigenvalues.cardinality(); ++i) {
                if (eigenvalues.get(i) > MIN_EIGENVALUE) {
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

        private SymmetricMatrix getNumericalHessianCentral() {
            double[][] hessian = new double[dim][dim];

            double[] oldUntransformedPosition = hessianWrtParameterProvider.getParameter().getParameterValues();
            double[] oldTransformedPosition = transform.transform(oldUntransformedPosition, 0, dim);

            double[][] gradientPlus = new double[dim][dim];
            double[][] gradientMinus = new double[dim][dim];
            
            double[] h = new double[dim];
            for (int i = 0; i < dim; i++) {
                h[i] = MachineAccuracy.SQRT_SQRT_EPSILON * (Math.abs(oldTransformedPosition[i]) + 1.0);
                hessianWrtParameterProvider.getParameter().setParameterValue(i, Math.exp(oldTransformedPosition[i] + h[i]));
                double[] tempGradient = hessianWrtParameterProvider.getGradientLogDensity();
                gradientPlus[i] = transform.updateGradientLogDensity(tempGradient, hessianWrtParameterProvider.getParameter().getParameterValues(),
                        0, dim);

                hessianWrtParameterProvider.getParameter().setParameterValue(i, Math.exp(oldTransformedPosition[i] - h[i]));
                tempGradient = hessianWrtParameterProvider.getGradientLogDensity();
                gradientMinus[i] = transform.updateGradientLogDensity(tempGradient, hessianWrtParameterProvider.getParameter().getParameterValues(),
                        0, dim);
                hessianWrtParameterProvider.getParameter().setParameterValue(i, Math.exp(oldTransformedPosition[i]));
            }
            
            for (int i = 0; i < dim; i++) {
                for (int j = i; j < dim; j++) {
                    hessian[j][i] = hessian[i][j] = (gradientPlus[i][j] - gradientMinus[i][j]) / (4.0 * h[j]) + (gradientPlus[j][i] - gradientMinus[j][i]) / (4.0 * h[i]);
                }
            }
//                double[] gradient = hessianWrtParameterProvider.getGradientLogDensity();
//                double[] unTransformedPosition = hessianWrtParameterProvider.getParameter().getParameterValues();
//                double[] diagonalHessian = hessianWrtParameterProvider.getDiagonalHessianLogDensity();
//                diagonalHessian = transform.updateDiagonalHessianLogDensity(diagonalHessian, gradient, unTransformedPosition,
//                        0, diagonalHessian.length);
//                double[][] hessianCheck = getNumericalHessianCheck(numeric1,
//                        transform.transform(hessianWrtParameterProvider.getParameter().getParameterValues(), 0, dim));
//                double[][] hessianCheck2 = getNumericalHessianCheckForward(dim,
//                        transform.transform(hessianWrtParameterProvider.getParameter().getParameterValues(), 0, dim));
            Algebra algebra = new Algebra();

            DoubleMatrix2D H = new DenseDoubleMatrix2D(hessian);
            RobustEigenDecomposition decomposition = new RobustEigenDecomposition(H);
            DoubleMatrix1D eigenvalues = decomposition.getRealEigenvalues();

            normalizeEigenvalues(eigenvalues);

            DoubleMatrix2D V = decomposition.getV();
            DoubleMatrix2D newHessian = algebra.mult(
                    algebra.mult(V, DoubleFactory2D.dense.diagonal(eigenvalues)),
                    algebra.inverse(V)
            );

            return new SymmetricMatrix(newHessian.toArray());
        }

        private double[][] getNumericalHessianCheckForward(int dim, double[] x) {
            double[][] hessian = new double[dim][dim];
            double[][] gradientMatrix = new double[dim][dim];
            final double[] gradient = transform.updateGradientLogDensity(hessianWrtParameterProvider.getGradientLogDensity(),
                    hessianWrtParameterProvider.getParameter().getParameterValues(), 0, dim);
            double[] oldUntransformedPosition = hessianWrtParameterProvider.getParameter().getParameterValues();
            double[] oldTransformedPosition = transform.transform(oldUntransformedPosition, 0, dim);
            double[] h = new double[dim];
            for (int i = 0; i < dim; i++) {
                h[i] = MachineAccuracy.SQRT_EPSILON * (Math.abs(x[i]) + 1.0);
                hessianWrtParameterProvider.getParameter().setParameterValue(i, Math.exp(oldTransformedPosition[i] + h[i]));
                double[] gradientTmp = hessianWrtParameterProvider.getGradientLogDensity();
                gradientMatrix[i] = transform.updateGradientLogDensity(gradientTmp, hessianWrtParameterProvider.getParameter().getParameterValues(),
                        0, dim);
                hessianWrtParameterProvider.getParameter().setParameterValue(i, Math.exp(oldTransformedPosition[i]));
            }
            for (int i = 0; i < dim; i++) {
                for (int j = i; j < dim; j++) {
                    hessian[i][j] = hessian[j][i] = (gradientMatrix[j][i] - gradient[i]) / (2.0 * h[j]) + (gradientMatrix[i][j] - gradient[j]) / (2.0 * h[i]);
                }
            }
            return hessian;
        }

        private double[][] getNumericalHessianCheck(MultivariateFunction f, double[] x) {
            double[][] hessian = new double[f.getNumArguments()][f.getNumArguments()];
            for (int i = 0; i < f.getNumArguments(); i++) {
                double hi = MachineAccuracy.SQRT_SQRT_EPSILON * (Math.abs(x[i]) + 1.0);
                double oldXi = x[i];
                double f__ = f.evaluate(x);
                x[i] = oldXi + hi;
                double fp_ = f.evaluate(x);
                x[i] = oldXi - hi;
                double fm_ = f.evaluate(x);
                x[i] = oldXi + 2.0 * hi;
                double fpp = f.evaluate(x);
                x[i] = oldXi - 2.0 * hi;
                double fmm = f.evaluate(x);
                hessian[i][i] = (-fpp + 16.0 * fp_ - 30.0 * f__ + 16.0 * fm_ - fmm) / (12.0 * hi * hi);
                for (int j = i + 1; j < f.getNumArguments(); j++) { //forward difference approximation
                    double hj = MachineAccuracy.SQRT_SQRT_EPSILON * (Math.abs(x[j]) + 1.0);
                    double oldXj = x[j];
                    x[i] = oldXi + hi;
                    x[j] = oldXj + hj;
                    fpp = f.evaluate(x);
                    x[i] = oldXi + hi;
                    x[j] = oldXj - hj;
                    double fpm = f.evaluate(x);
                    x[i] = oldXi - hi;
                    x[j] = oldXj + hj;
                    double fmp = f.evaluate(x);
                    x[i] = oldXi - hi;
                    x[j] = oldXj - hj;
                    fmm = f.evaluate(x);
                    x[i] = oldXi;
                    x[j] = oldXj;
                    hessian[i][j] = hessian[j][i] = (fpp - fpm - fmp + fmm) / (4.0 * hi * hj);
                }
            }
            return hessian;
        }

        private MultivariateFunction numeric1 = new MultivariateFunction() {
            @Override
            public double evaluate(double[] argument) {

                for (int i = 0; i < argument.length; ++i) {
                    hessianWrtParameterProvider.getParameter().setParameterValue(i, Math.exp(argument[i]));
                }

                return hessianWrtParameterProvider.getLikelihood().getLogLikelihood();
            }

            @Override
            public int getNumArguments() {
                return hessianWrtParameterProvider.getDimension();
            }

            @Override
            public double getLowerBound(int n) {
                return Double.NEGATIVE_INFINITY;
            }

            @Override
            public double getUpperBound(int n) {
                return Double.POSITIVE_INFINITY;
            }
        };
    }

}
