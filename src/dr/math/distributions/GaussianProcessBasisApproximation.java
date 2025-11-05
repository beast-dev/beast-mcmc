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

public class GaussianProcessBasisApproximation extends RandomFieldDistribution {

    public static final String TYPE = "GaussianProcessBasisApproximation";

    protected final int knots;
    protected final int dim;
    protected final double degree;
    protected final double[] times;
    protected final double boundary;
    private final Parameter meanParameter;
    private final Parameter marginalVarianceParameter;
    private final Parameter lengthScaleParameter;
    private final Parameter coefficientParameter;
    private final Parameter precisionParameter;

    private final double precisionValue;

    private final double[][] basisMatrix;
    private final double[] mean;
    private final double[] coefficient;
    private double[] centeredTimes;
    private double[] storedCoefficient;
    private boolean coefficientKnown;
    private boolean storedCoefficientKnown;
    private boolean basisMatrixKnown;
    private boolean timesKnown;
    private boolean centeredTimesKnown;
    private boolean meanKnown;
    private boolean precisionKnown;


    public GaussianProcessBasisApproximation(String name,
                                             int dim,
                                             int knots,
                                             double degree,
                                             double[] times,
                                             double boundary,
                                             Parameter mean,
                                             Parameter marginalVariance,
                                             Parameter lengthScale,
                                             Parameter coefficient,
                                             Parameter precision) {

        super(name);

        this.dim = dim;
        this.knots = knots;
        this.degree = degree;
        this.times = times;
        this.boundary = boundary;
        this.meanParameter = mean;
        this.marginalVarianceParameter = marginalVariance;
        this.lengthScaleParameter = lengthScale;
        this.precisionParameter = precision;
        this.coefficientParameter = coefficient;

        addVariable(meanParameter);
        addVariable(marginalVarianceParameter);
        addVariable(lengthScaleParameter);
        addVariable(precisionParameter);
        addVariable(coefficientParameter);



        this.coefficient = new double[knots];
        this.mean = new double[dim];
        this.centeredTimes = new double[dim];
        this.basisMatrix = new double[dim][knots];
        basisMatrixKnown = false;
        coefficientKnown = false;
        timesKnown = false;
        centeredTimesKnown = false;
        meanKnown = false;
        precisionKnown = false;
        precisionValue = 0;
    }

    @Override
    public double[] getMean() {

        double sum;
        getBasisMatrix(times);

        if (!meanKnown) {

            for (int i = 0; i < dim; i++) {
                sum = 0;
                for (int j = 0; j < knots; j++) {
                    sum += basisMatrix[i][j] * coefficientParameter.getParameterValue(j);
                    //System.err.println("BasisMat" + basisMatrix[i][j] + " coeff " + coefficientParameter.getParameterValue(j));
                }
                mean[i] = sum + meanParameter.getParameterValue(0);
                //System.err.println("Mean" + i + " at time " + mean[i]);
            }

            meanKnown = true;
        }

        return mean;

    }



    private double[] getCenteredTimes (double[] times) {

        if (!centeredTimesKnown) {
            double max = times[0];
            double min = times[0];
            double sum = times[0];
            for (int i = 1; i < times.length; i++) {
                sum += times[i];
                if (times[i] > max) {
                    max = times[i];
                }
                if (times[i] < min) {
                    min = times[i];
                }
            }

            double center = sum/times.length;
            double scale = 0.5 * (max - min);

            for (int i = 0; i <times.length; i++) {
                centeredTimes[i] = (times[i] - center)/scale;
                //System.err.println("Centered Time" + i + " at time " + centeredTimes[i]);
            }
            centeredTimesKnown = true;
        }

        return centeredTimes;
    }


    public static double getSpectralDensity(double x, double marginalVariance, double lengthScale, double degree) {

        double term = 0;

        if (degree == 1.5) {
            term = (marginalVariance * 12 * Math.sqrt(3) * Math.pow(lengthScale, -3))/(Math.pow((3/(lengthScale *
                    lengthScale)) + (x * x * 4 * Math.PI * Math.PI), 2));
        }

        if (degree == 2.5) {
            term = (marginalVariance * (400/3) * Math.sqrt(5) * Math.pow(lengthScale, -5))/(Math.pow((5/(lengthScale *
                    lengthScale)) + (x * x * 4 * Math.PI * Math.PI), 3));
        }

        if (degree == 0.5) {
            term = (2/lengthScale)/(1/(lengthScale * lengthScale) + 4 * Math.PI * Math.PI * x * x);
        }

        else if (degree == -1.0) {
            term = marginalVariance * Math.sqrt(2 *Math.PI) * lengthScale * Math.exp(-0.5 * x * x * lengthScale * lengthScale);
        }
        //System.err.println("Spectral density = " + term);
        return term;
    }

    public static double getSpectralDensityEigenValue(int j, double boundary) {
        return ((j + 1) * (j + 1) * Math.PI * Math.PI)/(4 * boundary * boundary);
    }

    private double getSpectralDensityEigenFunction(int j, double x, double boundary) {

        return Math.sin(Math.sqrt(getSpectralDensityEigenValue(j, boundary)) * (x + boundary))/Math.sqrt(boundary);
    }

    private void getBasisMatrix(double[]times) {


        double marginalVariance = marginalVarianceParameter.getParameterValue(0);
        double lengthScale = lengthScaleParameter.getParameterValue(0);
        double[] centeredTimes = getCenteredTimes(times);


        if (!basisMatrixKnown) {
            for (int i = 0; i < times.length; i++) {
                for (int j = 0; j < knots; j++) {
                    basisMatrix[i][j] = Math.sqrt(getSpectralDensity(Math.sqrt(getSpectralDensityEigenValue(j, boundary)),
                            marginalVariance, lengthScale, degree)) *
                            getSpectralDensityEigenFunction(j, centeredTimes[i], boundary);
                    //System.err.println("Sqrt S_theta" + Math.sqrt(getSpectralDensity(Math.sqrt(getSpectralDensityEigenValue(j)))) +
                    //        " phi " + getSpectralDensityEigenFunction(j, times[i]));
                }
            }
            basisMatrixKnown = true;
        }

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
                    getBasisMatrix(times);
                    double[] gradient =  gradLogPdf(dim, knots, (double[]) x, getMean(), precisionParameter.getParameterValue(0),
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

        final double[] delta = new double[dim];

        for (int i = 0; i < dim; i++) {
            delta[i] = x[i] - mean[i];
        }

        double SSE = 0.0;
        double sum;
        for (int i = 0; i < dim; i++) {
            sum = 0.0;
            for (int j = 0; j < knots; j++) {
                sum += basis[i][j] * coefficientParameter.getParameterValue(j);
            }
            delta[i] = sum;
        }

        for (int i = 0; i < dim; i++) {
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

//    public static double gradLogPdfWrtPrecision(double[] x, double[] mean, SymmetricTriDiagonalMatrix Q,
//                                                double precision, boolean isImproper) {
//        final int effectiveDim = isImproper ? x.length - 1 : x.length;
//        return 0.5 * ((effectiveDim - getSSE(x, mean, precision, Q))/precision);
//    }

//    public static double[] gradLogPdf(double[] x, double[] mean, double precision, SymmetricTriDiagonalMatrix Q) {
//
//        final int dim = x.length;
//
//        final double[] gradient = new double[dim];
//        final double[] delta = new double[dim];
//
//        for (int i = 0; i < dim; ++i) {
//            delta[i] = mean[i] - x[i];
//        }
//
//        gradient[0] = precision * (Q.diagonal[0] * delta[0] + Q.offDiagonal[0] * delta[1]);
//        for (int i = 1; i < dim - 1; ++i) {
//            gradient[i] = precision * (Q.offDiagonal[i - 1] * delta[i - 1] + Q.diagonal[i] * delta[i] + Q.offDiagonal[i] * delta[i + 1]);
//        }
//        gradient[dim - 1] = precision * (Q.offDiagonal[dim - 2] * delta[dim - 2] + Q.diagonal[dim - 1] * delta[dim - 1]);
//
//        return gradient;
//    }

//    public static double[][] hessianLogPdf(double[] x, double precision, SymmetricTriDiagonalMatrix Q) { // TODO test
//        final int dim = x .length;
//        final double[][] hessian = new double[dim][dim];
//
//        hessian[0][0] = -precision * Q.diagonal[0];
//        hessian[0][1] = -precision * Q.offDiagonal[0];
//
//        for (int i = 1; i < dim - 1; ++i) {
//            hessian[i][i - 1] = -precision * Q.offDiagonal[i - 1];
//            hessian[i][i]     = -precision * Q.diagonal[i];
//            hessian[i][i + 1] = -precision * Q.offDiagonal[i];
//        }
//
//        hessian[dim - 1][dim - 2] = -precision * Q.offDiagonal[dim - 2];
//        hessian[dim - 1][dim - 1] = -precision * Q.diagonal[dim - 1];
//
//        return hessian;
//    }

//    public static double[] diagonalHessianLogPdf(double[] x, double precision, SymmetricTriDiagonalMatrix Q) {
//        final int dim = x.length;
//        final double[] hessian = new double[dim];
//
//        System.arraycopy(Q.diagonal, 0, hessian, 0, dim);
//        // TODO do we not need to negate each element of hessian?
//
//        for (int i = 0; i < dim; i++) {
//            hessian[i] = hessian[i] * precision;
//        }
//
//        return hessian;
//    }

    // TODO Below is the relevant code from GMRFMultilocusSkyrideLikelihood for building a `SymmTridiagMatrix`
    // TODO `getFieldScalar` rescaling should be handled by `WeightsProvider`

//    protected double getFieldScalar() {
//        final double rootHeight;
//        if (rescaleByRootHeight) {
//            rootHeight = tree.getNodeHeight(tree.getRoot());
//        } else {
//            rootHeight = 1.0;
//        }
//        return rootHeight;
//    }
//
//    protected void setupGMRFWeights() {
//
//        setupSufficientStatistics();
//
//        //Set up the weight Matrix
//        double[] offdiag = new double[fieldLength - 1];
//        double[] diag = new double[fieldLength];
//
//        //First set up the offdiagonal entries;
//
//        if (!timeAwareSmoothing) {
//            for (int i = 0; i < fieldLength - 1; i++) {
//                offdiag[i] = -1.0;
//            }
//        } else {
//            for (int i = 0; i < fieldLength - 1; i++) {
//                offdiag[i] = -2.0 / (coalescentIntervals[i] + coalescentIntervals[i + 1]) * getFieldScalar();
//            }
//        }
//
//        //Then set up the diagonal entries;
//        for (int i = 1; i < fieldLength - 1; i++)
//            diag[i] = -(offdiag[i] + offdiag[i - 1]);
//
//        //Take care of the endpoints
//        diag[0] = -offdiag[0];
//        diag[fieldLength - 1] = -offdiag[fieldLength - 2];
//
//        weightMatrix = new SymmTridiagMatrix(diag, offdiag);
//    }
//
//    public SymmTridiagMatrix getScaledWeightMatrix(double precision) {
//        SymmTridiagMatrix a = weightMatrix.copy();
//        for (int i = 0; i < a.numRows() - 1; i++) {
//            a.set(i, i, a.get(i, i) * precision);
//            a.set(i + 1, i, a.get(i + 1, i) * precision);
//        }
//        a.set(fieldLength - 1, fieldLength - 1, a.get(fieldLength - 1, fieldLength - 1) * precision);
//        return a;
//    }
//
//    public SymmTridiagMatrix getScaledWeightMatrix(double precision, double lambda) {
//        if (lambda == 1)
//            return getScaledWeightMatrix(precision);
//
//        SymmTridiagMatrix a = weightMatrix.copy();
//        for (int i = 0; i < a.numRows() - 1; i++) {
//            a.set(i, i, precision * (1 - lambda + lambda * a.get(i, i)));
//            a.set(i + 1, i, a.get(i + 1, i) * precision * lambda);
//        }
//
//        a.set(fieldLength - 1, fieldLength - 1, precision * (1 - lambda + lambda * a.get(fieldLength - 1, fieldLength - 1)));
//        return a;
//    }
//
//    private DenseVector getMeanAdjustedGamma() {
//        DenseVector currentGamma = new DenseVector(popSizeParameter.getParameterValues());
//        updateGammaWithCovariates(currentGamma);
//        return currentGamma;
//    }
//
//    double getLogFieldLikelihood() {
//
//        DenseVector diagonal1 = new DenseVector(fieldLength);
//        DenseVector currentGamma = getMeanAdjustedGamma();
//
//        double currentLike = handleMissingValues();
//
//        SymmTridiagMatrix currentQ = getScaledWeightMatrix(precisionParameter.getParameterValue(0), lambdaParameter.getParameterValue(0));
//        currentQ.mult(currentGamma, diagonal1);
//
//        currentLike += 0.5 * (fieldLength - 1) * Math.log(precisionParameter.getParameterValue(0)) - 0.5 * currentGamma.dot(diagonal1);
//        if (lambdaParameter.getParameterValue(0) == 1) {
//            currentLike -= (fieldLength - 1) / 2.0 * LOG_TWO_TIMES_PI;
//        } else {
//            currentLike -= fieldLength / 2.0 * LOG_TWO_TIMES_PI;
//        }
//
//        return currentLike;
//    }

    @Override
    public int getDimension() { return dim; }

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
        if (variable == lengthScaleParameter || variable == marginalVarianceParameter) {
            basisMatrixKnown = false;
            meanKnown = false;
        } else if (variable == coefficientParameter || variable == meanParameter) {
            meanKnown = false;
        } else if (variable == precisionParameter) {
            precisionKnown = false;
        } else {
            throw new IllegalArgumentException("Unknown variable");
        }
    }

    @Override
    protected void storeState() { }

    @Override
    protected void restoreState() { // TODO cache mean
        meanKnown = false;
        basisMatrixKnown = false;
        precisionKnown = false;
    }

    @Override
    protected void acceptState() { }
}