/*
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.math;
import dr.inference.model.Parameter;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.MaxIterationsExceededException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.integration.TrapezoidIntegrator;

/**
 * @author Pratyusa Datta
 * @author Marc A. Suchard
 */

// still need to do Marc's TO DO'S

public class IntegratedTransformedSplines {


    private final Parameter coefficient;
    private final Parameter intercept;
    private final double[] knots;
    private final double lowerBoundary;
    private final double upperBoundary;
    private final double[] expandedKnots;
    private final int degree;
    private boolean expandedKnotsKnown;



    public IntegratedTransformedSplines(Parameter coefficient,
                                    Parameter intercept,
                                    double[] knots,
                                    double lowerBoundary,
                                    double upperBoundary,
                                    int degree) {

        if (coefficient.getDimension() != knots.length + degree + 1 && intercept == null) {
            throw new IllegalArgumentException("Coefficient dimension must be equal to number of knots + degree + 1 if intercept is null");
        }

        if (coefficient.getDimension() != knots.length + degree && intercept != null) {
            throw new IllegalArgumentException("Coefficient dimension must be equal to number of knots + degree if intercept is non-null");
        }

        this.coefficient = coefficient;
        this.intercept = intercept;
        this.knots = knots;
        this.lowerBoundary = lowerBoundary;
        this.upperBoundary = upperBoundary;
        this.expandedKnots = new double[knots.length + 2 * (degree + 1)];
        this.degree = degree;
        this.expandedKnotsKnown = false;

    }

    public Parameter getCoefficients() { return coefficient; }
    public Parameter getIntercept() { return intercept; }
    public int getCoefficientDim() {
        return coefficient.getDimension();
    }
    public boolean isInterceptNull() {
        return(intercept == null);
    }
    public double[] getCoefficientValues() {
        double[] values = new double[getCoefficientDim()];
        for (int i = 0; i < values.length; i++) {
            values[i] = coefficient.getParameterValue(i);
        }
        return values;
    }

    public void getExpandedKnots() {

        if (!expandedKnotsKnown) {

            for (int i = 0; i < degree + 1; i++) {
                expandedKnots[i] = lowerBoundary;
                expandedKnots[degree + knots.length + i + 1] = upperBoundary;
            }
            for (int i = 0; i < knots.length; i++) {
                expandedKnots[degree + i + 1] = knots[i];
            }
            expandedKnotsKnown = true;
        }
    }




    // TODO check
    public static double evaluatePolynomial(double x, double... coefficients) {
        int i = coefficients.length - 1;
        double value = coefficients[i];
        --i;

        for ( ; i >= 0; --i) {
            value = value * x + coefficients[i];
        }

        return value;
    }

    // TODO to integrate, could either
    // 1. transform coefficients[] --> integratedCoefficients[] and call evaluatePolynamial(), or
    // 2. fuse transformation and evaluation

    // TODO check
    public static double evaluatePolynomialIntegralEndPt(double x, double... coefficients) {
        int i = coefficients.length - 1;
        double value = coefficients[i] / (i + 1);
        --i;

        for ( ; i >= 0; --i) {
            value = value * x + coefficients[i] / (i + 1);
        }

        return value * x;
    }

    public static double evaluatePolynomialIntegral(double start, double end, double... coefficients) {
        return evaluatePolynomialIntegralEndPt(end, coefficients) -
                evaluatePolynomialIntegralEndPt(start, coefficients);
    }


    class Polynomial {

        final double[] coefficients;
        private final int degree;

        Polynomial(int degree) {
            this(new double[degree + 1], degree);
        }

        Polynomial(double[] coefficients, int degree) {

            assert coefficients.length >= degree + 1;

            this.coefficients = coefficients;
            this.degree = degree;
        }
    }

    private static final boolean TEST = false;

    public static double[] polynomialProduct1(double[] lhs, double[] rhs) {

        double[] product = new double[lhs.length + rhs.length - 1]; // TODO pass buffer

        for (int i = 0; i < lhs.length; ++i) {
            for (int j = 0; j < rhs.length; ++j) {
                product[i + j] += lhs[i] * rhs[j];
            }
        }

        return product;
    }

    public static double[] polynomialProduct2(double[] lhs, double[] rhs) {

        double[] product = new double[lhs.length + rhs.length - 1]; // TODO pass buffer

        for (int i = 0; i < product.length; ++i) {
            double sum = 0.0;
            for (int j = 0; j <= i; ++j) {
                sum += lhs[j] * rhs[i - j];
            }
            product[i] = sum;
        }

        return product;
    }

    public double getSplineBasis(int i, int d, double x, double[] knots) {
        if (d == 0) {
            if (x >= knots[i] && x < knots[i + 1]) {
                return 1.0;
            } else {
                return 0.0;
            }
        }

        double denom1 = knots[i + d] - knots[i];
        double denom2 = knots[i + d + 1] - knots[i + 1];

        double term1 = 0.0;
        double term2 = 0.0;

        if (denom1 != 0) {
            term1 = ((x - knots[i]) / denom1) * getSplineBasis(i, d - 1, x, knots);
        }
        if (denom2 != 0) {
            term2 = ((knots[i + d + 1] - x) / denom2) * getSplineBasis(i + 1, d - 1, x, knots);
        }

        return term1 + term2;
    }


    public double evaluateExpSpline(double x) {
        double sum = 0.0;
        getExpandedKnots();
        if (isInterceptNull()) {
            for (int i = 0; i < coefficient.getDimension(); i++) {
                sum += coefficient.getParameterValue(i) * getSplineBasis(i, degree, x, expandedKnots);
            }
        } else {
            for (int i = 0; i < coefficient.getDimension(); i++) {
                sum += coefficient.getParameterValue(i) * getSplineBasis(i + 1, degree, x, expandedKnots);
            }
            sum += intercept.getParameterValue(0);
        }


        return Math.exp(sum);
    }

    public double getExponentialSplinesIntegral(final double a, final double b)
            throws org.apache.commons.math.FunctionEvaluationException,
            org.apache.commons.math.MaxIterationsExceededException {

        UnivariateRealFunction f = new UnivariateRealFunction() {
            @Override
            public double value(double x) {
                return evaluateExpSpline(x);
            }
        };

        TrapezoidIntegrator integrator = new TrapezoidIntegrator();
        return integrator.integrate(f, a, b);
    }

    public double getIntegral (double start, double end) throws FunctionEvaluationException, MaxIterationsExceededException {

            return getExponentialSplinesIntegral(start, end);

    }

}





