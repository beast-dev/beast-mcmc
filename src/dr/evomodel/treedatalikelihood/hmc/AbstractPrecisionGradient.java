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
import dr.inference.model.CompoundSymmetricMatrix;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
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
    final CompoundSymmetricMatrix parameter;
    private final int dim;
    private final Parameterization parameterization;

    AbstractPrecisionGradient(ConjugateWishartStatisticsProvider wishartStatistics,
                              Likelihood likelihood,
                              CompoundSymmetricMatrix parameter) {
        this(wishartStatistics, likelihood, parameter, Parameterization.AS_PRECISION);
    }

    AbstractPrecisionGradient(ConjugateWishartStatisticsProvider wishartStatistics,
                              Likelihood likelihood,
                              CompoundSymmetricMatrix parameter, Parameterization parameterization) {
        assert parameter.asCorrelation()
                : "PrecisionGradient can only be applied to a CompoundSymmetricMatrix with off-diagonal as correlation.";

        this.wishartStatistics = wishartStatistics;
        this.likelihood = likelihood;
        this.parameter = parameter;
        this.dim = parameter.getColumnDimension();
        this.parameterization = parameterization;
    }

    enum Parameterization {
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
        return parameter;
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

    static int getVechDimension(int dim) {
        return dim * (dim + 1) / 2;
    }

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
        SymmetricMatrix correlationPrecision = new SymmetricMatrix(parameter.getCorrelationMatrix());
        double[] precisionDiagonal = parameter.getDiagonal();

        // TODO Chain-rule w.r.t. to parametrization

        if (CHECK_GRADIENT) {
            System.err.println("Analytic at: \n" + new Vector(parameter.getOffDiagonalParameter().getParameterValues())
                    + " " + new Vector(parameter.getDiagonal()));
        }

        double[] gradient = getGradientParameter(weightedSumOfSquares, numberTips,
                correlationPrecision, precisionDiagonal);

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
                                    SymmetricMatrix correlationPrecision,
                                    double[] precisionDiagonal) {

        // Gradient w.r.t. the correlation matrix (strictly upper diagonal)
        double[] gradientCorrelation = extractUpperTriangular((SymmetricMatrix) correlationPrecision.inverse());

        int k = 0;
        for (int i = 0; i < dim - 1; i++) {
            for (int j = i + 1; j < dim; j++) {
                gradientCorrelation[k] = numberTips * gradientCorrelation[k]
                        - weightedSumOfSquares.component(i, j) * Math.sqrt(precisionDiagonal[i] * precisionDiagonal[j]);
                k++;
            }
        }

        if (USE_CHAIN_RULE) {

            
            double[] diagD = sqrt(parameter.getDiagonal()); // TODO Reparameterize s.t. P = D^{1/2} C D^{1/2}, unless there is a more general case I am missing

            double[] vecV = new SymmetricMatrix(parameter.getParameterAsMatrix()).inverse().toArrayComponents(); // TODO Make inverse accessible from `parameter` to avoid recomputation

            WishartSufficientStatistics wss = wishartStatistics.getWishartStatistics();
            double[] vecS = wss.getScaleMatrix();
            int n = wss.getDf();

            double[] gradient = getGradientWrtPrecision(vecV, n, vecS);
            gradient = chainCorrelation(gradient, diagD);

            System.err.println("Correlation 2: " + new WrappedVector.Raw(gradient));
            System.err.println("Correlation 1: " + new WrappedVector.Raw(gradientCorrelation));
            System.err.println();
        }

        // TODO Handle chain-rule for parameterization

        // If necessary, apply chain rule to get the gradient w.r.t. cholesky of correlation matrix
        gradientCorrelation = parameter.updateGradientCorrelation(gradientCorrelation);

        return gradientCorrelation;
    }

    double[] chainDiagonal(double[] gradient, double[] diagD, double[] vecC, double[] diagQ) {
        MultivariateChainRule ruleD = new MultivariateChainRule.DecomposedDiagonals(diagD, vecC);
        gradient = ruleD.chainGradient(gradient);
        MultivariateChainRule ruleQ = new MultivariateChainRule.SqrtDiagonals(diagQ);
        return ruleQ.chainGradient(gradient);
    }

    double[] chainCorrelation(double[] gradient, double[] diagD) {
        MultivariateChainRule rule = new MultivariateChainRule.DecomposedCorrelation(diagD);
        return rule.chainGradient(gradient);
    }

    double[] getGradientWrtPrecision(double[] vecV, int n, double[] vecS) {

        assert vecV.length == dim * dim;
        assert vecS.length == dim * dim;
        assert n > 0;

        double[] gradient = new double[dim * dim];

        for (int i = 0; i < dim * dim; ++i) {
            gradient[i] = 0.5 * (n * vecV[i] - vecS[i]);
        }

        return gradient;
    }

    double[] flatten(double[][] matrix) {
        int dim = 0;
        for (int i = 0; i < matrix.length; ++i) {
            dim += matrix[i].length;
        }

        double[] result = new double[dim];

        int offset = 0;
        for (int i = 0; i < matrix.length; ++i) {
            System.arraycopy(matrix[i], 0, result, offset, matrix[i].length);
            offset += matrix[i].length;
        }

        return result;
    }

    double[] sqrt(double[] vector) {

        double[] result = new double[vector.length];

        for (int i = 0; i < result.length; ++i) {
            result[i] = Math.sqrt(vector[i]);
        }

        return result;
    }

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

            double[] vecC = flatten(parameter.getCorrelationMatrix());
            double[] diagQ = parameter.getDiagonal();
            double[] diagD = sqrt(diagQ); // TODO Reparameterize s.t. P = D^{1/2} C D^{1/2}, unless there is a more general case I am missing

            double[] vecV = new SymmetricMatrix(parameter.getParameterAsMatrix()).inverse().toArrayComponents(); // TODO Make inverse accessible from `parameter` to avoid recomputation
            
            WishartSufficientStatistics wss = wishartStatistics.getWishartStatistics();
            double[] vecS = wss.getScaleMatrix();
            int n = wss.getDf();

            double[] gradient = getGradientWrtPrecision(vecV, n, vecS);
            gradient = chainDiagonal(gradient, diagD, vecC, diagQ);

            System.err.println("vecC: " + new WrappedVector.Raw(vecC));
            System.err.println("Diagonal 2: " + new WrappedVector.Raw(gradient));
            System.err.println("Diagonal 1: " + new WrappedVector.Raw(gradientDiagonal));
            System.err.println();
        }
        return gradientDiagonal;
    }

    private static final boolean CHECK_GRADIENT = false;
    private static final boolean USE_CHAIN_RULE = true;

    interface MultivariateChainRule {

        double[] chainGradient(double[] lhs);

        class DecomposedCorrelation implements MultivariateChainRule {

            private final double[] diagD;
            private final int dim;

            DecomposedCorrelation(double[] diagD) {
                this.diagD = diagD;
                this.dim = diagD.length;
            }

            @Override
            public double[] chainGradient(double[] vecX) {

                assert vecX.length == dim * dim;

                double[] vechuGradient = new double[getVechuDimension(dim)];

                int k = 0;
                for (int i = 0; i < dim - 1; ++i) {
                    for (int j = i + 1; j < dim; ++j) {
                        vechuGradient[k] = 2.0 * vecX[i * dim + j] * diagD[i] * diagD[j];
                        ++k;
                    }
                }

                return vechuGradient;
            }
        }

        class DecomposedDiagonals implements  MultivariateChainRule {

            private final double[] diagD;
            private final double[] vecC;
            private final int dim;

            DecomposedDiagonals(double[] diagD, double[] vecC) {
                this.diagD = diagD;
                this.vecC = vecC;
                this.dim = diagD.length;

                assert vecC.length == dim * dim;
            }

            @Override
            public double[] chainGradient(double[] vecX) {

                assert vecX.length == dim * dim;

                double[] diagGradient = new double[dim];

                for (int i = 0; i < dim; ++i) {
                    double sum = 0.0;
                    for (int k = 0; k < dim; ++k) {
                        sum += 2.0 * vecX[i * dim + k] * diagD[k] * vecC[i * dim + k];
                    }

                    diagGradient[i] = sum;
                }

                return diagGradient;
            }
        }

        class SqrtDiagonals implements MultivariateChainRule {

            private final double[] diagQ;
            private final int dim;

            SqrtDiagonals(double[] diagQ) {
                this.diagQ = diagQ;
                this.dim = diagQ.length;
            }


            @Override
            public double[] chainGradient(double[] diagX) {

                assert diagX.length == dim;

                double[] diagGradient = new double[dim];

                for (int i = 0; i < dim; ++i) {
                    diagGradient[i] = 0.5 * diagX[i] / Math.sqrt(diagQ[i]);
                }

                return diagGradient;
            }
        }

        class Inverse implements MultivariateChainRule {

            private final double[] vecP;
            private final int dim;

            Inverse(double[] vecP) {
                this.vecP = vecP;
                this.dim = (int) Math.sqrt(vecP.length);
            }

            @Override
            public double[] chainGradient(double[] lhs) {
                return new double[0];
            }
        }
    }
}
