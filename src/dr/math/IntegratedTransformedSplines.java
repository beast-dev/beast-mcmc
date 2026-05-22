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

import java.util.Arrays;
import java.util.List;

/*
 * @author Pratyusa Datta
 * @author Marc A. Suchard
 */

public class IntegratedTransformedSplines {


    public static final String TRANSFORM_EXPONENTIAL = "exponential";
    public static final String TRANSFORM_SQUARED     = "squared";
    public static final String TRANSFORM_IDENTITY    = "identity";


    private final Parameter intercept;
    private final Parameter coefficient;
    private final double[] knots;
    private final double lowerBoundary;
    private final double upperBoundary;
    private final double[] expandedKnots;
    private final int degree;
    private final String transform;

    private final List<BSpline.PPoly> basis;
    private BSpline.PPoly splineFunction;
    private BSpline.PPoly splineFunctionSquared;
    private double[] cachedValues;



    public IntegratedTransformedSplines(
            Parameter coefficient,
            Parameter intercept,
            double[] knots,
            double lowerBoundary,
            double upperBoundary,
            int degree,
            String transform) {


        if (coefficient.getDimension() != knots.length + degree) {
            throw new IllegalArgumentException("Coefficient dimension must equal number of interior knots + degree");
        }

        if (!TRANSFORM_EXPONENTIAL.equals(transform) && !TRANSFORM_SQUARED.equals(transform) &&
                !TRANSFORM_IDENTITY.equals(transform)) {
            throw new IllegalArgumentException("transform must be \"" + TRANSFORM_EXPONENTIAL +
                    "\", \"" + TRANSFORM_SQUARED + "\", or \"" + TRANSFORM_IDENTITY + "\"");
        }

        this.intercept      = intercept;
        this.coefficient    = coefficient;
        this.knots          = knots;
        this.lowerBoundary  = lowerBoundary;
        this.upperBoundary  = upperBoundary;
        this.expandedKnots  = new double[knots.length + 2 * (degree + 1)];
        this.degree         = degree;
        this.transform      = transform;


        buildExpandedKnots();

        this.basis = BSpline.bSplineBasis(expandedKnots, degree);

        this.splineFunction = null;
        this.splineFunctionSquared = null;
        this.cachedValues = null;
    }


    public Parameter getIntercept()   { return intercept; }
    public Parameter getCoefficients()   { return coefficient; }
    public int getCoefficientDim() { return coefficient.getDimension(); }



    private void buildExpandedKnots() {
        for (int i = 0; i < degree + 1; i++) {
            expandedKnots[i] = lowerBoundary;
            expandedKnots[degree + knots.length + i + 1] = upperBoundary;
        }
        for (int i = 0; i < knots.length; i++) {
            expandedKnots[degree + i + 1] = knots[i];
        }
    }


    private double[] getCombinedParameter() {
        double[] combinedParameter = new double[1 + coefficient.getDimension()];
        combinedParameter[0] = intercept.getParameterValue(0);
        for (int i = 0; i < coefficient.getDimension(); i++) {
            combinedParameter[i + 1] = coefficient.getParameterValue(i);
        }
        return combinedParameter;
    }


    private void refreshSplineFunction() {
        double[] current = getCombinedParameter();

        if (splineFunction != null && Arrays.equals(current, cachedValues)) {
            return;
        }

        double beta0 = current[0];
        BSpline.PPoly f = BSpline.PPoly.constant(expandedKnots, beta0);
        for (int i = 1; i < basis.size(); i++) {
            f = BSpline.PPoly.add(f, basis.get(i).scale(current[i]));
        }
        splineFunction = f;

        if (TRANSFORM_SQUARED.equals(transform)) {
            splineFunctionSquared = f.square();
        }

        cachedValues = current;
    }


    public double evaluateSpline(double x) {
        refreshSplineFunction();
        return splineFunction.evaluate(x);
    }


    public double evaluateExpSpline(double x) {
        refreshSplineFunction();
        return Math.exp(splineFunction.evaluate(x));
    }

    public double evaluateSquaredSpline(double x) {
        refreshSplineFunction();
        return splineFunctionSquared.evaluate(x);
    }


    public double getExponentialSplinesIntegral(final double a, final double b)
            throws FunctionEvaluationException, MaxIterationsExceededException {

        refreshSplineFunction();
        final BSpline.PPoly f = splineFunction;

        UnivariateRealFunction integrand = new UnivariateRealFunction() {
            @Override
            public double value(double x) {
                return Math.exp(f.evaluate(x));
            }
        };

        TrapezoidIntegrator integrator = new TrapezoidIntegrator();
        return integrator.integrate(integrand, a, b);
    }


    public double getSquaredSplinesIntegral(double a, double b) {
        refreshSplineFunction();
        return splineFunctionSquared.integral(a, b);
    }


    public double getIdentitySplinesIntegral(double a, double b) {
        refreshSplineFunction();
        return splineFunction.integral(a, b);
    }


    public double getIntegral(double start, double end)
            throws FunctionEvaluationException, MaxIterationsExceededException {

        refreshSplineFunction();

        if (TRANSFORM_EXPONENTIAL.equals(transform)) {
            return getExponentialIntegral(start, end);
        } else if (TRANSFORM_SQUARED.equals(transform)) {
            return getSquaredIntegral(start, end);
        } else {
            return getIdentityIntegral(start, end);
        }
    }


    private double getExponentialIntegral(double start, double end)
            throws FunctionEvaluationException, MaxIterationsExceededException {

        if (end <= upperBoundary) {
            return getExponentialSplinesIntegral(start, end);
        } else if (start < upperBoundary) {
            double boundaryRate = evaluateExpSpline(upperBoundary);
            return getExponentialSplinesIntegral(start, upperBoundary)
                    + boundaryRate * (end - upperBoundary);
        } else {
            return evaluateExpSpline(upperBoundary) * (end - start);
        }
    }

    private double getSquaredIntegral(double start, double end) {

        if (end <= upperBoundary) {
            return getSquaredSplinesIntegral(start, end);
        } else if (start < upperBoundary) {
            double boundaryRate = evaluateSquaredSpline(upperBoundary);
            return getSquaredSplinesIntegral(start, upperBoundary)
                    + boundaryRate * (end - upperBoundary);
        } else {
            return evaluateSquaredSpline(upperBoundary) * (end - start);
        }
    }

    private double getIdentityIntegral(double start, double end) {

        if (end <= upperBoundary) {
            return getIdentitySplinesIntegral(start, end);
        } else if (start < upperBoundary) {
            double boundaryRate = evaluateSpline(upperBoundary);
            return getIdentitySplinesIntegral(start, upperBoundary)
                    + boundaryRate * (end - upperBoundary);
        } else {
            return evaluateSpline(upperBoundary) * (end - start);
        }
    }



    public static double evaluatePolynomial(double x, double... coefficients) {
        int i = coefficients.length - 1;
        double value = coefficients[i];
        --i;
        for (; i >= 0; --i) {
            value = value * x + coefficients[i];
        }
        return value;
    }

    public static double evaluatePolynomialIntegralEndPt(double x, double... coefficients) {
        int i = coefficients.length - 1;
        double value = coefficients[i] / (i + 1);
        --i;
        for (; i >= 0; --i) {
            value = value * x + coefficients[i] / (i + 1);
        }
        return value * x;
    }

    public static double evaluatePolynomialIntegral(double start, double end, double... coefficients) {
        return evaluatePolynomialIntegralEndPt(end, coefficients)
                - evaluatePolynomialIntegralEndPt(start, coefficients);
    }

    public static double[] polynomialProduct1(double[] lhs, double[] rhs) {
        double[] product = new double[lhs.length + rhs.length - 1];
        for (int i = 0; i < lhs.length; ++i)
            for (int j = 0; j < rhs.length; ++j)
                product[i + j] += lhs[i] * rhs[j];
        return product;
    }

    public static double[] polynomialProduct2(double[] lhs, double[] rhs) {
        double[] product = new double[lhs.length + rhs.length - 1];
        for (int i = 0; i < product.length; ++i) {
            double sum = 0.0;
            for (int j = 0; j <= i; ++j)
                sum += lhs[j] * rhs[i - j];
            product[i] = sum;
        }
        return product;
    }


    class Polynomial {
        final double[] coefficients;
        private final int degree;

        Polynomial(int degree) { this(new double[degree + 1], degree); }

        Polynomial(double[] coefficients, int degree) {
            assert coefficients.length >= degree + 1;
            this.coefficients = coefficients;
            this.degree       = degree;
        }
    }

}