/*
 * GlobalSigmoidSmoothFunction.java
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

package dr.evomodel.coalescent.smooth;

class GlobalSigmoidSmoothFunction {

    public double getSmoothValue(double x, double stepLocation, double preStepValue, double postStepValue, double smoothRate) {
        final double exponential = Math.exp(-smoothRate * (x - stepLocation));
        return (postStepValue - preStepValue) / (1 + exponential) + preStepValue;
    }

    public double getDerivative(double x, double stepLocation, double preStepValue, double postStepValue, double smoothRate) {
        final double exponential = Math.exp(-smoothRate * (x - stepLocation));
        final double result = Double.isInfinite(exponential) ? 0.0 : (smoothRate * (postStepValue - preStepValue) / (1.0 + 1.0 / exponential) / (1 + exponential));
        return result;
    }

    public double getLogDerivative(double x, double stepLocation, double preStepValue, double postStepValue, double smoothRate) {
        final double exponential = Math.exp(smoothRate * (x - stepLocation));
        final double result = Double.isInfinite(exponential) ? 0.0 : smoothRate * (postStepValue - preStepValue) / (postStepValue + preStepValue * exponential) / (1.0 + exponential);
        return result;
    }

    public double getSingleIntegration(double startTime, double endTime,
                                       double stepLocation, double smoothRate) {
        return (getLogOnePlusExponential(endTime - stepLocation, smoothRate) - getLogOnePlusExponential(startTime - stepLocation, smoothRate)) / smoothRate;
    }

    public double getSingleIntegrationDerivative(double startTime, double endTime, double stepLocation, double smoothRate) {
        final double exponent = Math.exp(smoothRate * (endTime - stepLocation));
        if (Double.isInfinite(exponent)) {
            return -1;
        } else {
            return 1.0 / (1.0 + Math.exp(-smoothRate * (startTime - stepLocation))) - exponent / (1.0 + exponent);
        }
    }

    public double getSingleIntegrationDerivativeWrtEndTime(double endTime, double stepLocation, double smoothRate) {
        final double exponent = Math.exp(smoothRate * (endTime - stepLocation));
        if (Double.isInfinite(exponent)) {
            return 1.0;
        } else {
            return exponent / (1.0 + exponent);
        }
    }

    public double getQuadraticIntegration(double startTime, double endTime,
                                          double stepLocation, double smoothRate) {
        return (getInverseOnePlusExponential(endTime - stepLocation, smoothRate) - getInverseOnePlusExponential(startTime - stepLocation, smoothRate) +
                getLogOnePlusExponential(endTime - stepLocation, smoothRate) - getLogOnePlusExponential(startTime - stepLocation, smoothRate)) / smoothRate;
    }

    public double getPairProductIntegration(double startTime, double endTime,
                                            double stepLocation1, double stepLocation2,
                                            double smoothRate) {
        final double firstTerm = endTime - startTime;
        final double secondTermMultiplier = getInverseOneMinusExponential(stepLocation2 - stepLocation1, smoothRate);
        final double thirdTermMultiplier = getInverseOneMinusExponential(stepLocation1 - stepLocation2, smoothRate);
        return firstTerm + secondTermMultiplier * doubleProductSingleRatio(startTime, endTime, stepLocation1, smoothRate)
                + thirdTermMultiplier * doubleProductSingleRatio(startTime, endTime, stepLocation2, smoothRate);
    }

    public double getPairProductIntegrationDerivativeWrtEndTime(double startTime, double endTime,
                                            double stepLocation1, double stepLocation2,
                                            double smoothRate) {
        final double gridNodeExponential = Math.exp(smoothRate * (stepLocation2 - stepLocation1));
        final double first = 1.0;
        final double secondTermMultiplier = 1.0 / (1.0 - Math.exp(smoothRate * (stepLocation2 - stepLocation1)));
        final double thirdTermMultiplier = 1.0 / (1.0 - Math.exp(smoothRate * (stepLocation1 - stepLocation2)));

        return first + secondTermMultiplier * doubleProductSingleRatioDerivativeWrtEndTime(startTime, endTime, stepLocation1, smoothRate)
                + thirdTermMultiplier * doubleProductSingleRatioDerivativeWrtEndTime(startTime, endTime, stepLocation2, smoothRate);
    }

    public double getPairProductIntegrationDerivative(double startTime, double endTime,
                                            double stepLocation1, double stepLocation2,
                                            double smoothRate) {

        final double gridNodeExponential = Math.exp(smoothRate * (stepLocation2 - stepLocation1));
        final double startNodeExponential = Math.exp(smoothRate * (startTime - stepLocation1));
        final double endNodeExponential = Math.exp(smoothRate * (endTime - stepLocation1));
        final double first = Double.isInfinite(gridNodeExponential) ? 0.0 : - gridNodeExponential / ((1 - gridNodeExponential) * (1 - gridNodeExponential))
                * (getLogOnePlusExponential(stepLocation1, endTime, smoothRate) - getLogOnePlusExponential(stepLocation1, startTime, smoothRate));
        final double second = 1.0 / (1.0 - gridNodeExponential) * (1.0 / (1.0 + endNodeExponential) - 1.0 / (1.0 + startNodeExponential));
        final double third = Double.isInfinite(gridNodeExponential) ? 0.0 : gridNodeExponential / ((1.0 - gridNodeExponential) * (1.0 - gridNodeExponential))
                * (getLogOnePlusExponential(stepLocation2, endTime, smoothRate) - getLogOnePlusExponential(stepLocation2, startTime, smoothRate));
        return first + second + third;
    }

    private double getLogOnePlusExponential(double t1, double t2, double smoothRate) {
        final double exponential = Math.exp(smoothRate * (t1 - t2));
        if (Double.isInfinite(exponential)) {
            return smoothRate * (t1 - t2);
        } else {
            return Math.log(1.0 + Math.exp(smoothRate * (t1 - t2)));
        }
    }

    private double doubleProductSingleRatio(double startTime, double endTime,
                                            double stepLocation, double smoothRate) {
        final double exponent = Math.exp(smoothRate * (stepLocation - startTime));
        final double others = Math.log(1 + Math.exp(smoothRate * (stepLocation - endTime))) / smoothRate;
        if (Double.isInfinite(exponent)) {
            return others - (stepLocation - startTime);
        } else {
            return others - Math.log(1.0 + exponent) / smoothRate;
        }
    }

    private double doubleProductSingleRatioDerivativeWrtEndTime(double startTime, double endTime,
                                            double stepLocation, double smoothRate) {
        final double exponent = Math.exp(smoothRate * (stepLocation - endTime));
        if (Double.isInfinite(exponent)) {
            return -1.0;
        } else {
            return - exponent / (1.0 + exponent);
        }
    }


    public double getInverseOneMinusExponential(double x, double smoothRate) {
        final double exponential = Math.exp(smoothRate * x);
        if (Double.isInfinite(exponential)) {
            return 0.0;
        } else {
            return 1.0 / (1.0 - exponential);
        }
    }

    public double getInverseOnePlusExponential(double x, double smoothRate) {
        final double exponential = Math.exp(smoothRate * x);
        if (Double.isInfinite(exponential)) {
            return 0.0;
        } else {
            return 1.0 / (1.0 + exponential);
        }
    }

    public double getLogOnePlusExponential(double x, double smoothRate) {
        final double exponential = Math.exp(smoothRate * x);
        if (Double.isInfinite(exponential)) {
            return smoothRate * x;
        } else {
            return Math.log(1 + exponential);
        }
    }

    public double getTripleProductIntegration(double startTime, double endTime,
                                              double stepLocation1, double stepLocation2, double stepLocation3,
                                              double smoothRate) {
        final double result = (getInverseOneMinusExponential(stepLocation2 - stepLocation1, smoothRate)
                * getInverseOneMinusExponential(stepLocation3 - stepLocation1, smoothRate) *
                (getLogOnePlusExponential(stepLocation1 - endTime, smoothRate) - getLogOnePlusExponential(stepLocation1 - startTime, smoothRate))
                + getInverseOneMinusExponential(stepLocation3 - stepLocation2, smoothRate)
                * getInverseOneMinusExponential(stepLocation1 - stepLocation2, smoothRate) *
                (getLogOnePlusExponential(stepLocation2 - endTime, smoothRate) - getLogOnePlusExponential(stepLocation2 - startTime, smoothRate))
                + getInverseOneMinusExponential(stepLocation1 - stepLocation3, smoothRate)
                * getInverseOneMinusExponential(stepLocation2 - stepLocation3, smoothRate) *
                (getLogOnePlusExponential(stepLocation3 - endTime, smoothRate) - getLogOnePlusExponential(stepLocation3 - startTime, smoothRate))) / smoothRate
                + endTime - startTime;
        return result;
    }

    public double getTripleProductWithQuadraticIntegration(double startTime, double endTime,
                                                           double stepLocation12, double stepLocation3,
                                                           double smoothRate) {
        final double first = endTime - startTime;
        final double second = getInverseOneMinusExponential(stepLocation3 - stepLocation12, smoothRate)
                * (getInverseOnePlusExponential(stepLocation12 - startTime, smoothRate) -
                getInverseOnePlusExponential(stepLocation12 - endTime, smoothRate));
        final double thirdMultiplier = ((2.0 - getInverseOneMinusExponential(stepLocation3 - stepLocation12, smoothRate)) * getInverseOneMinusExponential(stepLocation3 - stepLocation12, smoothRate));
        final double third = thirdMultiplier * (getLogOnePlusExponential(stepLocation12 - endTime, smoothRate) - getLogOnePlusExponential(stepLocation12 - startTime, smoothRate));
        final double fourthMultiplier = getInverseOneMinusExponential(stepLocation12 - stepLocation3, smoothRate) * getInverseOneMinusExponential(stepLocation12 - stepLocation3, smoothRate);
        final double fourth = fourthMultiplier * (getLogOnePlusExponential(stepLocation3 - endTime, smoothRate) - getLogOnePlusExponential(stepLocation3 - startTime, smoothRate));
        return first + (second + third + fourth) / smoothRate;
    }


}