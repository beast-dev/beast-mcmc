/*
 * AbstractPrecisionGradient.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.treedatalikelihood.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.math.matrixAlgebra.Vector;
import dr.math.matrixAlgebra.WrappedVector;
import dr.xml.Reportable;

import static dr.math.matrixAlgebra.SymmetricMatrix.extractUpperTriangular;

/**
 * @author Paul Bastide
 * @author Marc A. Suchard
 */
public abstract class AbstractPrecisionGradient implements GradientWrtParameterProvider, Reportable {

    private final ConjugateWishartStatisticsProvider wishartStatistics;
    final Likelihood likelihood;
    final CompoundSymmetricMatrix compoundSymmetricMatrix;
    private final int dim;
    private final Parametrization parametrization;

    private final MatrixParameterInterface precision;
    private final MatrixParameterInterface variance;

    AbstractPrecisionGradient(ConjugateWishartStatisticsProvider wishartStatistics,
                              Likelihood likelihood,
                              MatrixParameterInterface parameter) {

        this.precision = parameter;

        if (parameter instanceof CachedMatrixInverse) {

            this.compoundSymmetricMatrix = (CompoundSymmetricMatrix)
                    ((CachedMatrixInverse) parameter).getBaseParameter();
            this.variance = compoundSymmetricMatrix;
            this.parametrization = Parametrization.AS_VARIANCE;

        } else if (parameter instanceof CompoundSymmetricMatrix) {

            this.compoundSymmetricMatrix = (CompoundSymmetricMatrix) parameter;
            this.variance = new CachedMatrixInverse("", parameter);
            this.parametrization = Parametrization.AS_PRECISION;

        } else {
            throw new IllegalArgumentException("Unimplemented type");
        }

        assert compoundSymmetricMatrix.asCorrelation()
                : "PrecisionGradient can only be applied to a CompoundSymmetricMatrix with off-diagonal as correlation.";

        this.wishartStatistics = wishartStatistics;
        this.likelihood = likelihood;
        this.dim = parameter.getColumnDimension();
    }

    enum Parametrization {
        AS_PRECISION {
            @Override
            double[] chainRule(double[] x) {
                return x; // Do nothing
            }
        },
        AS_VARIANCE {
            @Override
            double[] chainRule(double[] x) {
                // TODO Handle matrix-inverse
                return x;
            }
        };

        abstract double[] chainRule(double[] x);
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return compoundSymmetricMatrix;
    }

    @Override
    public int getDimension() {
        return getParameter().getDimension();
    }

    int getDimensionCorrelation() {
        return dim * (dim - 1) / 2;
    }

    int getDimensionDiagonal() {
        return dim;
    }

//    static int getVechDimension(int dim) {
//        return dim * (dim + 1) / 2;
//    }

    static int getVechuDimension(int dim) {
        return dim * (dim - 1) / 2;
    }

    @Override
    public double[] getGradientLogDensity() {

        // Statistics
        WishartSufficientStatistics statistics = wishartStatistics.getWishartStatistics();
        SymmetricMatrix weightedSumOfSquares = new SymmetricMatrix(statistics.getScaleMatrix(), dim);
        int numberTips = statistics.getDf();

        // TODO For non-matrix-normal models, use sum of gradients w.r.t. branch-specific precisions

        // parameters
        SymmetricMatrix correlation = new SymmetricMatrix(compoundSymmetricMatrix.getCorrelationMatrix());
        double[] diagonal = compoundSymmetricMatrix.getDiagonal();

        // TODO Chain-rule w.r.t. to parametrization

        if (CHECK_GRADIENT) {
            System.err.println("Analytic at: \n" + new Vector(compoundSymmetricMatrix.getOffDiagonalParameter().getParameterValues())
                    + " " + new Vector(compoundSymmetricMatrix.getDiagonal()));
        }

        double[] gradient = getGradientParameter(weightedSumOfSquares, numberTips,
                correlation, diagonal);

        if (CHECK_GRADIENT) {
            System.err.println(checkNumeric(gradient));
        }

        return gradient;
    }

    String getReportString(double[] analytic, double[] numeric) {

        return getClass().getCanonicalName() + "\n" +
                "analytic: " + new Vector(analytic) +
                "\n" +
                "numeric: " + new Vector(numeric) +
                "\n";
    }

    String getReportString(double[] analytic, double[] numeric, double[] numericTrans) {

        return getClass().getCanonicalName() + "\n" +
                "analytic: " + new Vector(analytic) +
                "\n" +
                "numeric (no Cholesky): " + new Vector(numeric) +
                "\n" +
                "numeric (with Cholesky): " + new Vector(numericTrans) +
                "\n";
    }

    abstract String checkNumeric(double[] analytic);

    @Override
    public String getReport() {
        return checkNumeric(getGradientLogDensity());
    }

    abstract double[] getGradientParameter(SymmetricMatrix weightedSumOfSquares,
                                           int numberTips,
                                           SymmetricMatrix correlationPrecision,
                                           double[] precisionDiagonal);

    // Gradient w.r.t. correlation
    double[] getGradientCorrelation(SymmetricMatrix weightedSumOfSquares,
                                    int numberTips,
                                    SymmetricMatrix correlation,
                                    double[] diagonal) {

        // Gradient w.r.t. the correlation matrix (strictly upper diagonal)
        double[] gradientCorrelation = extractUpperTriangular((SymmetricMatrix) correlation.inverse());

        int k = 0;
        for (int i = 0; i < dim - 1; i++) {
            for (int j = i + 1; j < dim; j++) {
                gradientCorrelation[k] = numberTips * gradientCorrelation[k]
                        - weightedSumOfSquares.component(i, j) * Math.sqrt(diagonal[i] * diagonal[j]);
                k++;
            }
        }

        if (USE_CHAIN_RULE) {


            double[] diagQ = compoundSymmetricMatrix.getDiagonal();

            double[] vecV = flatten(variance.getParameterAsMatrix());
            double[] vecP = flatten(precision.getParameterAsMatrix());

            WishartSufficientStatistics wss = wishartStatistics.getWishartStatistics();
            double[] vecS = wss.getScaleMatrix();
            int n = wss.getDf();

            double[] gradient = getGradientWrtPrecision(vecV, n, vecS);
            gradient = chainCorrelation(gradient, vecP, diagQ);

            System.err.println("Correlation 2: " + new WrappedVector.Raw(gradient));
            System.err.println("Correlation 1: " + new WrappedVector.Raw(gradientCorrelation));
            System.err.println();
        }

        // TODO Handle chain-rule for parametrization

        // If necessary, apply chain rule to get the gradient w.r.t. cholesky of correlation matrix
        gradientCorrelation = compoundSymmetricMatrix.updateGradientCorrelation(gradientCorrelation);

        return gradientCorrelation;
    }

    private double[] chainDiagonal(double[] gradient, double[] vecP, double[] vecV,
                                   double[] vecC, double[] diagQ) {

        if (parametrization == Parametrization.AS_VARIANCE) { // TODO Delegate
            MultivariateChainRule ruleI = new MultivariateChainRule.Inverse(vecP, vecV);
            gradient = ruleI.chainGradient(gradient);
        }

//        MultivariateChainRule ruleD = new MultivariateChainRule.DecomposedDiagonals(diagD, vecC);
//        gradient = ruleD.chainGradient(gradient);
        MultivariateChainRule ruleQ = new MultivariateChainRule.DecomposedSqrtDiagonals(diagQ, vecC);
        return ruleQ.chainGradient(gradient);
    }

    private double[] chainCorrelation(double[] gradient, double[] vecP, double[] diagQ) {
        if (parametrization == Parametrization.AS_VARIANCE) { // TODO Delegate
            MultivariateChainRule ruleI = new MultivariateChainRule.InverseGeneral(vecP);
            gradient = ruleI.chainGradient(gradient);
        }

        MultivariateChainRule rule = new MultivariateChainRule.DecomposedCorrelation(diagQ);
        return rule.chainGradient(gradient);
    }

    private double[] getGradientWrtPrecision(double[] vecV, int n, double[] vecS) {

        assert vecV.length == dim * dim;
        assert vecS.length == dim * dim;
        assert n > 0;

        double[] gradient = new double[dim * dim];

        for (int i = 0; i < dim * dim; ++i) {
            gradient[i] = 0.5 * (n * vecV[i] - vecS[i]);
        }

        return gradient;
    }

    private double[] flatten(double[][] matrix) {
        int dim = 0;
        for (double[] vector : matrix) {
            dim += vector.length;
        }

        double[] result = new double[dim];

        int offset = 0;
        for (double[] vector : matrix) {
            System.arraycopy(vector, 0, result, offset, vector.length);
            offset += vector.length;
        }

        return result;
    }

//    private double[] sqrt(double[] vector) {
//
//        double[] result = new double[vector.length];
//
//        for (int i = 0; i < result.length; ++i) {
//            result[i] = Math.sqrt(vector[i]);
//        }
//
//        return result;
//    }

    // Gradient w.r.t. diagonal
    double[] getGradientDiagonal(SymmetricMatrix weightedSumOfSquares,
                                 int numberTips,
                                 SymmetricMatrix correlationPrecision,
                                 double[] precisionDiagonal) {

        // Gradient w.r.t. to the diagonal values of the precision
        double[] gradientDiagonal = new double[dim];

        for (int i = 0; i < dim; i++) {
            // Product
            double innerProduct = 0.0;
            for (int j = 0; j < dim; j++) {
                innerProduct += correlationPrecision.component(i, j) * weightedSumOfSquares.component(i, j)
                        * Math.sqrt(precisionDiagonal[j] / precisionDiagonal[i]);
            }
            // diagonal
            gradientDiagonal[i] = numberTips * 0.5 / precisionDiagonal[i] - 0.5 * innerProduct;
        }

        if (USE_CHAIN_RULE) {

            double[] vecC = flatten(compoundSymmetricMatrix.getCorrelationMatrix());
            double[] diagQ = compoundSymmetricMatrix.getDiagonal();
//            double[] diagD = sqrt(diagQ); //DONE //TODO Reparameterize s.t. P = D^{1/2} C D^{1/2}, unless there is a more general case I am missing

            double[] vecV = flatten(variance.getParameterAsMatrix());
            double[] vecP = flatten(precision.getParameterAsMatrix());

            WishartSufficientStatistics wss = wishartStatistics.getWishartStatistics();
            double[] vecS = wss.getScaleMatrix();
            int n = wss.getDf();

            double[] gradient = getGradientWrtPrecision(vecV, n, vecS);
            gradient = chainDiagonal(gradient, vecP, vecV, vecC, diagQ);

            System.err.println("vecC: " + new WrappedVector.Raw(vecC));
            System.err.println("Diagonal 2: " + new WrappedVector.Raw(gradient));
            System.err.println("Diagonal 1: " + new WrappedVector.Raw(gradientDiagonal));
            System.err.println();
        }
        return gradientDiagonal;
    }

    private static final boolean CHECK_GRADIENT = true;
    private static final boolean USE_CHAIN_RULE = true;

    interface MultivariateChainRule {

        double[] chainGradient(double[] lhs);

        class DecomposedCorrelation implements MultivariateChainRule {

            private final double[] diagQ;
            private final int dim;

            DecomposedCorrelation(double[] diagQ) {
                this.diagQ = diagQ;
                this.dim = diagQ.length;
            }

            @Override
            public double[] chainGradient(double[] vecX) {

                assert vecX.length == dim * dim;

                double[] vechuGradient = new double[getVechuDimension(dim)];

                int k = 0;
                for (int i = 0; i < dim - 1; ++i) {
                    for (int j = i + 1; j < dim; ++j) {
                        vechuGradient[k] = 2.0 * vecX[i * dim + j] * Math.sqrt(diagQ[i] * diagQ[j]);
                        ++k;
                    }
                }

                return vechuGradient;
            }
        }

//        class DecomposedDiagonals implements  MultivariateChainRule {
//
//            private final double[] diagD;
//            private final double[] vecC;
//            private final int dim;
//
//            DecomposedDiagonals(double[] diagD, double[] vecC) {
//                this.diagD = diagD;
//                this.vecC = vecC;
//                this.dim = diagD.length;
//
//                assert vecC.length == dim * dim;
//            }
//
//            @Override
//            public double[] chainGradient(double[] vecX) {
//
//                assert vecX.length == dim * dim;
//
//                double[] diagGradient = new double[dim];
//
//                for (int i = 0; i < dim; ++i) {
//                    double sum = 0.0;
//                    for (int k = 0; k < dim; ++k) {
//                        sum += 2.0 * vecX[i * dim + k] * diagD[k] * vecC[i * dim + k];
//                    }
//
//                    diagGradient[i] = sum;
//                }
//
//                return diagGradient;
//            }
//        }
//
//        class SqrtDiagonals implements MultivariateChainRule {
//
//            private final double[] diagQ;
//            private final int dim;
//
//            SqrtDiagonals(double[] diagQ) {
//                this.diagQ = diagQ;
//                this.dim = diagQ.length;
//            }
//
//
//            @Override
//            public double[] chainGradient(double[] diagX) {
//
//                assert diagX.length == dim;
//
//                double[] diagGradient = new double[dim];
//
//                for (int i = 0; i < dim; ++i) {
//                    diagGradient[i] = 0.5 * diagX[i] / Math.sqrt(diagQ[i]);
//                }
//
//                return diagGradient;
//            }
//        }

        class DecomposedSqrtDiagonals implements MultivariateChainRule {

            private final double[] diagQ;
            private final double[] vecC;
            private final int dim;

            DecomposedSqrtDiagonals(double[] diagQ, double[] vecC) {
                this.diagQ = diagQ;
                this.vecC = vecC;
                this.dim = diagQ.length;
            }


            @Override
            public double[] chainGradient(double[] vecX) {

                assert vecX.length == dim * dim;

                double[] diagGradient = new double[dim];

                for (int i = 0; i < dim; ++i) {
                    double sum = 0.0;
                    for (int k = 0; k < dim; ++k) {
                        sum += vecX[i * dim + k] * Math.sqrt(diagQ[k] / diagQ[i]) * vecC[i * dim + k];
                    }

                    diagGradient[i] = sum;
                }

                return diagGradient;
            }
        }

        class Inverse implements MultivariateChainRule {

            private final double[] vecP;
            private final double[] vecV;
            private final int dim;

            Inverse(double[] vecP, double[] vecV) {
                this.vecP = vecP;
                this.vecV = vecV;
                this.dim = (int) Math.sqrt(vecP.length);
            }

            @Override
            public double[] chainGradient(double[] lhs) {

                assert lhs.length == dim * dim;

                double[] gradient = new double[dim * dim];

                for (int i = 0; i < dim * dim; ++i) {
                    gradient[i] = -lhs[i] * vecP[i] / vecV[i];
                }

                return gradient;
            }
        }

        class InverseGeneral implements MultivariateChainRule {

            private final double[] vecP;
            private final int dim;

            InverseGeneral(double[] vecP) {
                this.vecP = vecP;
                this.dim = (int) Math.sqrt(vecP.length);
            }

            @Override
            public double[] chainGradient(double[] lhs) {

                assert lhs.length == dim * dim;

                double[] gradient = new double[dim * dim];

                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++) {
                        double sum = 0.0;
                        for (int k = 0; k < dim; k++) {
                            for (int l = 0; l < dim; l++) {
                                sum += vecP[i * dim + l] * lhs[l * dim + k] * vecP[k * dim + j];
                            }
                        }
                        gradient[i * dim + j] = -sum;
                    }
                }

                return gradient;
            }
        }
    }
}
