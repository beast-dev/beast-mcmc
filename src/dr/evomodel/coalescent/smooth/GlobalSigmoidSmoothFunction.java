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
        final double result = Double.isInfinite(exponential) ? 0.0 : (smoothRate * (postStepValue - preStepValue) * exponential / (1 + exponential) / (1 + exponential));
        return result;
    }

    public double getLogDerivative(double x, double stepLocation, double preStepValue, double postStepValue, double smoothRate) {
        final double exponential = Math.exp(smoothRate * (x - stepLocation));
        final double result = Double.isInfinite(exponential) ? 0.0 : smoothRate * (postStepValue - preStepValue) / (postStepValue + preStepValue * exponential) / (1.0 + exponential);
        return result;
    }

    public double getSingleIntegration(double startTime, double endTime,
                                       double stepLocation, double smoothRate) {
        final double exponent = Math.exp(smoothRate * (endTime - stepLocation));
        if (Double.isInfinite(exponent)) {
            return endTime - stepLocation;
        } else {
            final double ratio = (1 + exponent) / (1 + Math.exp(smoothRate * (startTime - stepLocation)));
            return Math.log(ratio) / smoothRate;
        }
    }

    public double getSingleIntegrationDerivative(double startTime, double endTime, double stepLocation, double smoothRate) {
        final double exponent = Math.exp(smoothRate * (endTime - stepLocation));
        if (Double.isInfinite(exponent)) {
            return -1;
        } else {
            return 1.0 / (1.0 + Math.exp(-smoothRate * (startTime - stepLocation))) - 1.0 / (1.0 + Math.exp(-smoothRate * (endTime - stepLocation)));
        }
    }

    public double getPairProductIntegration(double startTime, double endTime,
                                            double stepLocation1, double stepLocation2,
                                            double smoothRate) {
        final double firstTerm = endTime - startTime;
        final double secondTermMultiplier = 1.0 / (1.0 - Math.exp(smoothRate * (stepLocation2 - stepLocation1)));
        final double thirdTermMultiplier = 1.0 / (1.0 - Math.exp(smoothRate * (stepLocation1 - stepLocation2)));
        return firstTerm + secondTermMultiplier * doubleProductSingleRatio(startTime, endTime, stepLocation1, smoothRate)
                + thirdTermMultiplier * doubleProductSingleRatio(startTime, endTime, stepLocation2, smoothRate);
    }

    public double getPairProductIntegrationDerivative(double nodeTime, double startTime, double endTime,
                                            double stepLocation1, double stepLocation2,
                                            double smoothRate) {
        final double inverse = 1.0 / (1.0 - Math.exp(smoothRate * (stepLocation2 - stepLocation1)));
        final double gridNodeExponential = Math.exp(smoothRate * (stepLocation2 - nodeTime));
        final double startNodeExponential = Math.exp(smoothRate * (nodeTime - startTime));
        final double endNodeExponential = Math.exp(smoothRate * (nodeTime - endTime));
        final double first = - gridNodeExponential / ((1 - gridNodeExponential) * (1 - gridNodeExponential))
                * (getLogOnePlusExponential(nodeTime, endTime, smoothRate) - getLogOnePlusExponential(nodeTime, startTime, smoothRate));
        final double second = 1.0 / (1.0 - gridNodeExponential) * (1.0 / (1.0 + endNodeExponential) - 1.0 / (1.0 + startNodeExponential));
        final double third = gridNodeExponential / ((1.0 - gridNodeExponential) * (1.0 - gridNodeExponential))
                * (getLogOnePlusExponential(stepLocation2, endTime, smoothRate) - getLogOnePlusExponential(stepLocation2, startTime, smoothRate));
        return first + second + third;
    }

    private double getLogOnePlusExponential(double t1, double t2, double smoothRate) {
        return Math.log(1.0 + Math.exp(smoothRate * (t1 - t2)));
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

    public double getTripleProductIntegration(double startTime, double endTime,
                                              double stepLocation1, double preStepValue1, double postStepValue1,
                                              double stepLocation2, double preStepValue2, double postStepValue2,
                                              double stepLocation3, double preStepValue3, double postStepValue3,
                                              double smoothRate) {

        int numPreStepZeros = (preStepValue1 == 0 ? 1 : 0) + (preStepValue2 == 0 ? 1 : 0) + (preStepValue3 == 0 ? 1 : 0);
        int numSameStepLocations = (stepLocation1 == stepLocation2 ? 1 : 0) +
                (stepLocation2 == stepLocation3 ? 1 : 0) + (stepLocation3 == stepLocation1 ? 1 : 0);

        if (numPreStepZeros == 3) {
            if (numSameStepLocations == 0) {
                return tripleProductIntegrationCase1(endTime, stepLocation1, postStepValue1, stepLocation2, postStepValue2, stepLocation3, postStepValue3, smoothRate)
                        - tripleProductIntegrationCase1(startTime, stepLocation1, postStepValue1, stepLocation2, postStepValue2, stepLocation3, postStepValue3, smoothRate);
            } else if (numSameStepLocations == 1) { //TODO: not sure if want to generalize for any stepLocation pair being the same
                return tripleProductIntegrationCase2(endTime, stepLocation1, postStepValue1, stepLocation3, postStepValue3, smoothRate)
                        - tripleProductIntegrationCase2(startTime, stepLocation1, postStepValue1, stepLocation3, postStepValue3, smoothRate);
            }
        }
        if (numPreStepZeros == 2) {
            if (numSameStepLocations == 0) {
                return tripleProductIntegrationCase3(endTime, stepLocation1, postStepValue1, stepLocation2, postStepValue2, stepLocation3, preStepValue3, postStepValue3, smoothRate)
                        - tripleProductIntegrationCase3(startTime, stepLocation1, postStepValue1, stepLocation2, postStepValue2, stepLocation3, preStepValue3, postStepValue3, smoothRate);
            } else if (numSameStepLocations == 1) {
                return tripleProductIntegrationCase4(endTime, stepLocation1, postStepValue1, postStepValue2, stepLocation3, preStepValue3, postStepValue3, smoothRate)
                        - tripleProductIntegrationCase4(startTime, stepLocation1, postStepValue1, postStepValue2, stepLocation3, preStepValue3, postStepValue3, smoothRate);
            }
        }
        if (numPreStepZeros == 1) {
            if (numSameStepLocations == 0) {
                return tripleProductIntegrationCase5(endTime, stepLocation1, postStepValue1, stepLocation2, preStepValue2, postStepValue2,
                        stepLocation3, preStepValue3, postStepValue3, smoothRate)
                        - tripleProductIntegrationCase5(startTime, stepLocation1, postStepValue1, stepLocation2, preStepValue2, postStepValue2,
                        stepLocation3, preStepValue3, postStepValue3, smoothRate);
            } else if (numSameStepLocations == 1) {
                return tripleProductIntegrationCase6(endTime, stepLocation1, postStepValue1, preStepValue2, postStepValue2, stepLocation3, preStepValue3, postStepValue3, smoothRate)
                        - tripleProductIntegrationCase6(startTime, stepLocation1, postStepValue1, preStepValue2, postStepValue2, stepLocation3, preStepValue3, postStepValue3, smoothRate);
            }
        }
        throw new RuntimeException("Integration case not implemented!");
    }

    private double tripleProductIntegrationCase6(double time,
                                                 double stepLocation1, double postStepValue1,
                                                 double preStepValue2, double postStepValue2,
                                                 double stepLocation3, double preStepValue3, double postStepValue3,
                                                 double smoothRate) {
        final double exponentialStepLocation1 = Math.exp(smoothRate * stepLocation1);
        final double exponentialStepLocation3 = Math.exp(smoothRate * stepLocation3);
        final double smoothRateTime = smoothRate * time;
        final double exponentialTime = Math.exp(smoothRate * time);

        final double numerator = postStepValue1 * ((exponentialStepLocation1 - exponentialStepLocation3) * (preStepValue3 * postStepValue2 * exponentialTime * exponentialStepLocation3
        + preStepValue2 * (postStepValue3 * exponentialTime * exponentialStepLocation1 - preStepValue3 * exponentialTime * exponentialStepLocation3)
        + postStepValue2 * postStepValue3 * (exponentialStepLocation1 * exponentialStepLocation1 * smoothRateTime
                - exponentialTime * exponentialStepLocation3 * smoothRateTime
                - exponentialStepLocation1 * exponentialStepLocation3 * smoothRateTime
                + exponentialTime * exponentialStepLocation1 * (smoothRateTime - 1)))
        + (exponentialTime + exponentialStepLocation1) * (preStepValue3 * postStepValue2 * exponentialStepLocation3 * exponentialStepLocation3 - preStepValue2 * (preStepValue3 - postStepValue3) * exponentialStepLocation1 * exponentialStepLocation3
        + postStepValue2 * postStepValue3 * (exponentialStepLocation1 * exponentialStepLocation1 - 2 * exponentialStepLocation1 * exponentialStepLocation3)) * Math.log(exponentialStepLocation1 / exponentialTime + 1)
        - (preStepValue3 - postStepValue3) * exponentialStepLocation3 * (exponentialTime + exponentialStepLocation1)
                * (-preStepValue2 * exponentialStepLocation1 + postStepValue2 * exponentialStepLocation3) *
                Math.log(exponentialStepLocation3 / exponentialTime + 1));
        final double denominator = (exponentialTime + exponentialStepLocation1) * (exponentialStepLocation1 - exponentialStepLocation3) * (exponentialStepLocation1 - exponentialStepLocation3) * smoothRate;
        return numerator / denominator;
    }

    private double tripleProductIntegrationCase5(double time,
                                                 double stepLocation1, double postStepValue1,
                                                 double stepLocation2, double preStepValue2, double postStepValue2,
                                                 double stepLocation3, double preStepValue3, double postStepValue3,
                                                 double smoothRate) {
        final double exponentialStepLocation1 = Math.exp(smoothRate * stepLocation1);
        final double exponentialStepLocation2 = Math.exp(smoothRate * stepLocation2);
        final double exponentialStepLocation3 = Math.exp(smoothRate * stepLocation3);
        final double exponentialTime = Math.exp(smoothRate * time);

        final double numerator = postStepValue1 * (
                (exponentialStepLocation2 - exponentialStepLocation3) * (preStepValue2 * exponentialStepLocation2 - postStepValue2 * exponentialStepLocation1) * (preStepValue3 * exponentialStepLocation3 - postStepValue3 * exponentialStepLocation1) * Math.log(exponentialTime + exponentialStepLocation1)
                        + (preStepValue2 - postStepValue2) * exponentialStepLocation2 * (exponentialStepLocation1 - exponentialStepLocation3) * (postStepValue3 * exponentialStepLocation2 - preStepValue3 * exponentialStepLocation3) * Math.log(exponentialTime + exponentialStepLocation2)
                + (preStepValue3 - postStepValue3) * exponentialStepLocation3 * (exponentialStepLocation1 - exponentialStepLocation2) * (preStepValue2 * exponentialStepLocation2 - postStepValue2 * exponentialStepLocation3) * Math.log(exponentialTime + exponentialStepLocation3)
                );
        final double denominator = smoothRate * (exponentialStepLocation1 - exponentialStepLocation2) * (exponentialStepLocation1 - exponentialStepLocation3) * (exponentialStepLocation2 - exponentialStepLocation3);
        return numerator / denominator;
    }

    private double tripleProductIntegrationCase4(double time,
                                                 double stepLocation1, double postStepValue1,
                                                 double postStepValue2,
                                                 double stepLocation3, double preStepValue3, double postStepValue3,
                                                 double smoothRate) {
        final double exponentialStepLocation1 = Math.exp(smoothRate * stepLocation1);
        final double exponentialStepLocation3 = Math.exp(smoothRate * stepLocation3);
        final double exponentialTime = Math.exp(smoothRate * time);
        final double smoothRateTime = smoothRate * time;

        final double numerator = postStepValue1 * postStepValue2 * (
                (1-exponentialStepLocation3/exponentialStepLocation1)*(preStepValue3*exponentialTime/exponentialStepLocation1
                        + postStepValue3*smoothRateTime*(exponentialStepLocation1/exponentialStepLocation3 - exponentialTime/exponentialStepLocation1 - 1 + exponentialTime/exponentialStepLocation3)
                        - postStepValue3*exponentialTime/exponentialStepLocation3)
                + (1 + exponentialTime/exponentialStepLocation1)*(preStepValue3*exponentialStepLocation3/exponentialStepLocation1 + postStepValue3*(exponentialStepLocation1/exponentialStepLocation3-2))*Math.log(1+exponentialStepLocation1/exponentialTime)
                - (preStepValue3 - postStepValue3)*exponentialStepLocation3/exponentialStepLocation1*(1+exponentialTime/exponentialStepLocation1)*Math.log(1+exponentialStepLocation3/exponentialTime));
        final double denominator = smoothRate * (1 + exponentialTime/exponentialStepLocation1)*(1-exponentialStepLocation3/exponentialStepLocation1)*(exponentialStepLocation1/exponentialStepLocation3-1);
        return numerator / denominator;
    }


    private double tripleProductIntegrationCase3(double time,
                                                 double stepLocation1, double postStepValue1,
                                                 double stepLocation2, double postStepValue2,
                                                 double stepLocation3, double preStepValue3, double postStepValue3,
                                                 double smoothRate) {
        final double exponentialStepLocation1 = Math.exp(smoothRate * stepLocation1);
        final double exponentialStepLocation2 = Math.exp(smoothRate * stepLocation2);
        final double exponentialStepLocation3 = Math.exp(smoothRate * stepLocation3);
        final double exponentialTime = Math.exp(smoothRate * time);

        final double numerator = postStepValue1 * postStepValue2 *
                (exponentialStepLocation1 * (exponentialStepLocation2 - exponentialStepLocation3) * (postStepValue3 * exponentialStepLocation1 - preStepValue3 * exponentialStepLocation3) * Math.log(exponentialStepLocation1 / exponentialTime + 1)
                - exponentialStepLocation2 * (exponentialStepLocation1 - exponentialStepLocation3) * (postStepValue3 * exponentialStepLocation2 - preStepValue3 * exponentialStepLocation3) * Math.log(exponentialStepLocation2 / exponentialTime + 1)
                - exponentialStepLocation3 * exponentialStepLocation3 * (preStepValue3 - postStepValue3) * (exponentialStepLocation1 - exponentialStepLocation2) * Math.log(exponentialStepLocation3 / exponentialTime + 1));
        final double denominator = smoothRate * (exponentialStepLocation1 - exponentialStepLocation2) * (exponentialStepLocation1 - exponentialStepLocation3) * (exponentialStepLocation2 - exponentialStepLocation3);
        return numerator / denominator + postStepValue1 * postStepValue2 * postStepValue3 * time;
    }

    private double tripleProductIntegrationCase1(double time,
                                                 double stepLocation1, double postStepValue1,
                                                 double stepLocation2, double postStepValue2,
                                                 double stepLocation3, double postStepValue3,
                                                 double smoothRate) {
        final double exponentialStepLocation1 = Math.exp(smoothRate * stepLocation1);
        final double exponentialStepLocation2 = Math.exp(smoothRate * stepLocation2);
        final double exponentialStepLocation3 = Math.exp(smoothRate * stepLocation3);
        final double exponentialTime = Math.exp(smoothRate * time);

        final double postValueProduct = postStepValue1 * postStepValue2 * postStepValue3;

        final double numerator = postValueProduct * (exponentialStepLocation1 * exponentialStepLocation1 * (exponentialStepLocation2 - exponentialStepLocation3) * Math.log(exponentialStepLocation1 / exponentialTime + 1)
        + exponentialStepLocation2 * exponentialStepLocation2 * (exponentialStepLocation3 - exponentialStepLocation1) * Math.log(exponentialStepLocation2 / exponentialTime + 1)
        + exponentialStepLocation3 * exponentialStepLocation3 * (exponentialStepLocation1 - exponentialStepLocation2) * Math.log(exponentialStepLocation3 / exponentialTime + 1));
        final double denominator = smoothRate * (exponentialStepLocation1 - exponentialStepLocation2) * (exponentialStepLocation1 - exponentialStepLocation3) * (exponentialStepLocation2 - exponentialStepLocation3);
        return numerator / denominator + postValueProduct * time;
    }

    private double tripleProductIntegrationCase2(double time,
                                                 double stepLocation1, double postStepValue1,
                                                 double stepLocation3, double postStepValue3,
                                                 double smoothRate) {
        final double exponentialStepLocation1 = Math.exp(smoothRate * stepLocation1);
        final double exponentialStepLocation3 = Math.exp(smoothRate * stepLocation3);
        final double exponentialTime = Math.exp(smoothRate * time);
        final double postValueProduct = postStepValue1 * postStepValue1 * postStepValue3;
        return postValueProduct * (-exponentialStepLocation1 / (smoothRate * (exponentialStepLocation1 / exponentialTime + 1) * (exponentialStepLocation1 - exponentialStepLocation3))
        + (exponentialStepLocation1 * (exponentialStepLocation1 - 2*exponentialStepLocation3) * Math.log(exponentialStepLocation1 / exponentialTime + 1)
                + exponentialStepLocation3 * exponentialStepLocation3 * Math.log(exponentialStepLocation3 / exponentialTime + 1) )/ (smoothRate * (exponentialStepLocation1 - exponentialStepLocation3) * (exponentialStepLocation1 - exponentialStepLocation3))
        + time);
    }

}