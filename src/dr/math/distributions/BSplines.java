/*
 * GaussianProcessBasisApproximation.java
 *
 * Copyright (c) 2002-2025 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.math.distributions;
import dr.inference.model.*;

/**
 * @author Marc Suchard
 * @author Pratyusa Data
 */

public class BSplines extends RandomFieldDistribution{

    public static final String TYPE = "BSplines";

    protected final double[] knots;
    protected final int degree;
    protected final double[] times;
    protected final double lowerBoundary;
    protected final double upperBoundary;
    private final Parameter coefficientParameter;
    private final Parameter precisionParameter;

    private MatrixParameter basisMatrixParameter;
    private final double[][] basisMatrix;
    private final double[] mean;
    private final double[] coefficient;
    private final double[] expandedKnots;
    private boolean basisMatrixKnown;
    private boolean timesKnown;
    private boolean expandedKnotsKnown;
    private boolean meanKnown;
    private boolean precisionKnown;


    public BSplines(String name,
                    double[] knots,
                    double[] times,
                    int degree,
                    double lowerBoundary,
                    double upperBoundary,
                    Parameter coefficient,
                    Parameter precision) {

        super(name);

        this.knots = knots;
        this.times = times;
        this.degree = degree;
        this.lowerBoundary = lowerBoundary;
        this.upperBoundary = upperBoundary;
        this.coefficientParameter = coefficient;
        this.precisionParameter = precision;

        addVariable(coefficientParameter);
        addVariable(precisionParameter);

        this.basisMatrixParameter = null;
        this.coefficient = new double[knots.length + degree - 1];
        this.expandedKnots = new double[2 * degree + knots.length];
        this.mean = new double[times.length];
        this.basisMatrix = new double[times.length][knots.length + degree - 1];
        meanKnown = false;
        basisMatrixKnown = false;
        expandedKnotsKnown = false;
        precisionKnown = false;
    }


    @Override
    public double[] getMean() {

        double sum;
        getBasisMatrix();

        if (!meanKnown) {

            for (int i = 0; i < times.length; i++) {
                sum = 0;
                for (int j = 0; j < knots.length + degree - 1; j++) {
                    sum += basisMatrix[i][j] * coefficientParameter.getParameterValue(j);
                    //System.err.println("BasisMat" + basisMatrix[i][j] + " coeff " + coefficientParameter.getParameterValue(j));
                }
                //System.err.println("Mean" + i + " at time " + mean[i]);
                mean[i] = sum;
            }

            meanKnown = true;
        }

        return mean;

    }

    public int getDegree() {
        return degree;
    }

    public double getLowerBoundary() {
        return lowerBoundary;
    }

    public double getUpperBoundary() {
        return upperBoundary;
    }

    public double[] getCoefficient() {
        int dim = coefficientParameter.getDimension();
        double[] coefficient = new double[dim];
        for (int i = 0; i < dim; i++) {
            coefficient[i] = coefficientParameter.getParameterValue(i);
        }
        return coefficient;
    }

    public double[] getExpandedKnots() {
        if (!expandedKnotsKnown) {
            for (int i = 0; i < degree; i++) {
                expandedKnots[i] = lowerBoundary;
                expandedKnots[degree + knots.length + i] = upperBoundary;
            }
            for (int i = 0; i < knots.length; i++) {
                expandedKnots[degree + i] = knots[i];
            }
            expandedKnotsKnown = true;
        }
        return expandedKnots;
    }


    public static double getSplineBasis(int i, int d, double x, double[] knots) {
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


    public void getBasisMatrix() {

        getExpandedKnots();
        if (!basisMatrixKnown) {
            int n = times.length;
            int I = knots.length + degree - 1;

            for (int r = 0; r < n; r++) {
                for (int i = 0; i < I; i++) {
                    basisMatrix[r][i] = getSplineBasis(i, degree, times[r], expandedKnots);
                }
            }
            basisMatrixKnown = true;
        }
    }

    public MatrixParameterInterface getBasisMatrixParameter() {

        getBasisMatrix();

        final int nRow = times.length;
        final int nCol = knots.length + degree - 1;

        if (basisMatrixParameter != null &&
                basisMatrixParameter.getRowDimension() == nRow &&
                basisMatrixParameter.getColumnDimension() == nCol) {
            return basisMatrixParameter;
        }

        basisMatrixParameter = new MatrixParameter("basisMatrix");
        basisMatrixParameter.setDimensions(nRow, nCol);

        for (int i = 0; i < nRow; i++) {
            for (int j = 0; j < nCol; j++) {
                basisMatrixParameter.setParameterValueQuietly(i, j, basisMatrix[i][j]);
            }
        }


        return basisMatrixParameter;
    }

    @Override
    public GradientProvider getGradientWrt(Parameter parameter) {

        if (parameter == coefficientParameter) {
            return new GradientProvider() {
                @Override
                public int getDimension() {
                    return coefficientParameter.getDimension();
                }

                @Override
                public double[] getGradientLogDensity(Object x) {
                    getBasisMatrix();
                    double[] gradient =  gradLogPdf(times.length, degree + knots.length - 1, (double[]) x,
                            getMean(), precisionParameter.getParameterValue(0),
                            basisMatrix);
                    return gradient;
                }
            };
        }
        // TO DO: gradient with respect to lengthScale, marginalVariance and precision of noise
        else {
            throw new RuntimeException("Unknown parameter");
        }
    }



    public String getType() {
        return TYPE;
    }


    @Override
    public double[][] getScaleMatrix() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public Variable<Double> getLocationVariable() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
//    public double logPdf(double[] x) {
//        getBasisMatrix(getTimes());
//        return 0.5 * ( - dim * Math.log(2 * Math.PI) + dim * Math.log(precisionParameter.getParameterValue(0))
//                    - getSSE(x, getMean(), basisMatrix) * precisionParameter.getParameterValue(0));
//    }

    public double logPdf(double[] x) {
        return MultivariateNormalDistribution.logPdf(x, getMean(), precisionParameter.getParameterValue(0), 1);
    }

    public double getSSE(double[] x, double[] mean, double[][] basis) {

        final double[] delta = new double[times.length];

        for (int i = 0; i < times.length; i++) {
            delta[i] = x[i] - mean[i];
        }

        double SSE = 0.0;
        double sum;
        for (int i = 0; i < times.length; i++) {
            sum = 0.0;
            for (int j = 0; j < knots.length + degree - 1; j++) {
                sum += basis[i][j] * coefficientParameter.getParameterValue(j);
            }
            delta[i] = sum;
        }

        for (int i = 0; i < times.length; i++) {
            SSE += delta[i] * delta[i];
        }

        return SSE;
    }

    public static double[] gradLogPdf(int dim, int knots, double[] x, double[] mean, double precision, double[][] basisMatrix) {

        final double[] gradient = new double[knots];
        final double[] delta = new double[dim];
        double sum;
        for (int i = 0; i < dim; ++i) {
            delta[i] = x[i] - mean[i];
        }


        for (int i = 0; i < knots; ++i) {
            sum = 0;
            for (int j = 0; j < dim; j++) {
                sum += basisMatrix[j][i] * delta[j];
            }
            gradient[i] = sum * precision;
        }


        return gradient;
    }


    @Override
    public int getDimension() { return times.length; }

    @Override
    public double[] getGradientLogDensity(Object x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] getDiagonalHessianLogDensity(Object x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[][] getHessianLogDensity(Object x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] nextRandom() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        throw new IllegalArgumentException("Unknown model");
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == coefficientParameter) {
            meanKnown = false;
        } else if (variable == precisionParameter) {
            precisionKnown = false;}
        else {
            throw new IllegalArgumentException("Unknown variable");
        }
    }

    @Override
    protected void storeState() { }

    @Override
    protected void restoreState() { // TODO cache mean
        meanKnown = false;
        precisionKnown = false;
    }

    @Override
    protected void acceptState() { }
}