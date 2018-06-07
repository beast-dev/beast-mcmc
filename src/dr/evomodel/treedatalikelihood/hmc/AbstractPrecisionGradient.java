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

        // TODO Handle chain-rule for parameterization

        // If necessary, apply chain rule to get the gradient w.r.t. cholesky of correlation matrix
        gradientCorrelation = parameter.updateGradientCorrelation(gradientCorrelation);

        return gradientCorrelation;
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

        // TODO Handle chain-rule for parameterization

        return gradientDiagonal;
    }

    private static final boolean CHECK_GRADIENT = false;
}
