/*
 * TimeVaryingBranchRateModel.java
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
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.MaxIterationsExceededException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.integration.TrapezoidIntegrator;
import dr.math.distributions.BSplines;

/**
 * @author Pratyusa Datta
 * @author Marc A. Suchard
 */

// still need to do Marc's TO DO'S

public class IntegratedSquaredSplines {


    public final double[] coefficient;
    private final double[] knots;
    private final double[] expandedKnots;
    private final int degree;
    private boolean expandedKnotsKnown;
    private double[][] temp;
    private final SplinesIntegral splinesIntegral;

    public IntegratedSquaredSplines(double[] coefficient,
                                    double[] knots,
                                    int degree) {

        if (degree != 2 && degree != 3) {
            throw new IllegalArgumentException("Only degree 2 and 3 supported");
        }
        this.coefficient = coefficient;
        this.knots = knots;
        this.expandedKnots = new double[knots.length + 2 * degree];
        this.degree = degree;
        this.expandedKnotsKnown = false;

        this.temp = new double[coefficient.length][coefficient.length];
        this.splinesIntegral = (degree == 2) ? new Quadratic() : new Cubic();
    }


    public void getExpandedKnots() {
        if (!expandedKnotsKnown) {
            for (int i = 0; i < degree; i++) {
                expandedKnots[i] = knots[0];
                expandedKnots[degree + knots.length + i] = knots[knots.length - 1];
            }
            for (int i = 0; i < knots.length; i++) {
                expandedKnots[degree + i] = knots[i];
            }
        }
    }

    private interface SplinesIntegral {
        double getDiagonal(int i, double start, double end);
        double getOffDiagonalOrder1(int i, double start, double end);
        double getOffDiagonalOrder2(int i, double start, double end);
        double getOffDiagonalOrder3(int i, double start, double end); // may return 0 for degree 2
    }

    private class Quadratic implements SplinesIntegral {
        public double getDiagonal(int i, double start, double end) {
            double value = 0;
            double low;
            double up;
            low = Math.max(start, expandedKnots[i]);
            up = Math.min(end, expandedKnots[i + 1]);
            if (expandedKnots[i + 2] != expandedKnots[i] && expandedKnots[i + 1] != expandedKnots[i] && low < up) {
                value += ((Math.pow(up - expandedKnots[i], 5) - Math.pow(low - expandedKnots[i], 5))/
                        (5 * Math.pow((expandedKnots[i + 2] - expandedKnots[i]) * (expandedKnots[i + 1] - expandedKnots[i]), 2)));
            }
       /* System.out.println("low:" + low);
        System.out.println("up:" + up);
        System.out.println("value:" + value);*/
            low = Math.max(start, expandedKnots[i + 1]);
            up = Math.min(end, expandedKnots[i + 2]);
            if (low < up) {
                value += getSquaredQuadraticIntegral(getQuadratic(i), low, up);
            }
     /*   System.out.println("low:" + low);
        System.out.println("up:" + up);
        System.out.println("value:" + value);*/
            low = Math.max(start, expandedKnots[i + 2]);
            up = Math.min(end, expandedKnots[i + 3]);
            if (expandedKnots[i + 3] != expandedKnots[i + 1] && expandedKnots[i + 3] != expandedKnots[i + 2] && low < up) {
                value += ((Math.pow(up - expandedKnots[i + 3], 5) - Math.pow(low - expandedKnots[i + 3], 5))/
                        (5 * Math.pow((expandedKnots[i + 3] - expandedKnots[i + 1]) * (expandedKnots[i + 3] - expandedKnots[i + 2]), 2)));
            }
        /*System.out.println("low:" + low);
        System.out.println("up:" + up);
        System.out.println("value:" + value);*/
            return value;
        }

        public double getOffDiagonalOrder1(int i, double start, double end) {
            double value = 0;
            double low;
            double up;
            double[] y = new double[3];
            double c;
            if (expandedKnots[i + 3] != expandedKnots[i + 1] && expandedKnots[i + 2] != expandedKnots[i + 1]) {
                c = 1/((expandedKnots[i + 3] - expandedKnots[i + 1]) * (expandedKnots[i + 2] - expandedKnots[i + 1]));
            } else {
                c = 0;
            }
            y[0] = (expandedKnots[i + 1] * expandedKnots[i + 1]) * c;
            y[1] = -2 * expandedKnots[i + 1] * c;
            y[2] = c;
            low = Math.max(start, expandedKnots[i + 1]);
            up = Math.min(end, expandedKnots[i + 2]);
            if (low < up) {
                value += getQuadraticProductIntegral(getQuadratic(i), y, low, up);
            }


            if (expandedKnots[i + 3] != expandedKnots[i + 1] && expandedKnots[i + 3] != expandedKnots[i + 2]) {
                c = 1/((expandedKnots[i + 3] - expandedKnots[i + 1]) * (expandedKnots[i + 3] - expandedKnots[i + 2]));
            } else {
                c = 0;
            }
            y[0] = (expandedKnots[i + 3] * expandedKnots[i + 3]) * c;
            y[1] = -2 * expandedKnots[i + 3] * c;
            y[2] = c;
            low = Math.max(start, expandedKnots[i + 2]);
            up = Math.min(end, expandedKnots[i + 3]);
            if (low < up) {
                value += getQuadraticProductIntegral(getQuadratic(i + 1), y, low, up);
            }
            return value;
        }

        public double getOffDiagonalOrder2(int i, double start, double end) {
            double value = 0;
            double low;
            double up;
            double[] y = new double[3];
            double[] x = new double[3];
            double c;
            if (expandedKnots[i + 4] != expandedKnots[i + 2] && expandedKnots[i + 3] != expandedKnots[i + 2]) {
                c = 1/((expandedKnots[i + 4] - expandedKnots[i + 2]) * (expandedKnots[i + 3] - expandedKnots[i + 2]));
            } else {
                c = 0;
            }
            x[0] = (expandedKnots[i + 2] * expandedKnots[i + 2]) * c;
            x[1] = (-2 * expandedKnots[i + 2]) * c;
            x[2] = c;
            if (expandedKnots[i + 3] != expandedKnots[i + 1] && expandedKnots[i + 3] != expandedKnots[i + 2]) {
                c = 1/((expandedKnots[i + 3] - expandedKnots[i + 1]) * (expandedKnots[i + 3] - expandedKnots[i + 2]));
            } else {
                c = 0;
            }
            y[0] = expandedKnots[i + 3] * expandedKnots[i + 3] * c;
            y[1] = (-2 * expandedKnots[i + 3]) * c;
            y[2] = c;
            low = Math.max(start, expandedKnots[i + 2]);
            up = Math.min(end, expandedKnots[i + 3]);
            if (low < up) {
                value += getQuadraticProductIntegral(x, y, low, up);
            }

            return value;
        }

        @Override
        public double getOffDiagonalOrder3(int i, double start, double end) {
            return 0;
        }

        private double[] getQuadratic(int i) {
            double[] values = new double[3];
            double x;
            double y;
            if (expandedKnots[i + 2] != expandedKnots[i] && expandedKnots[i + 2] != expandedKnots[i + 1]) {
                x = 1/((expandedKnots[i + 2] - expandedKnots[i]) * (expandedKnots[i + 2] - expandedKnots[i + 1]));
            }
            else {
                x = 0;
            }
            if (expandedKnots[i + 3] != expandedKnots[i + 1] && expandedKnots[i + 2] != expandedKnots[i + 1]) {
                y = 1/((expandedKnots[i + 3] - expandedKnots[i + 1]) * (expandedKnots[i + 2] - expandedKnots[i + 1]));
            }
            else {
                y = 0;
            }
            values[2] = -x - y;
            values[1] = x * (expandedKnots[i + 2] + expandedKnots[i]) + y * (expandedKnots[i + 3] + expandedKnots[i + 1]);
            values[0] = -x * expandedKnots[i] * expandedKnots[i + 2] - expandedKnots[i + 3] * expandedKnots[i + 1] * y;
            return values;
        }

        private double getSquaredQuadraticIntegral(double[] coeff, double low, double up) {
            java.util.function.DoubleFunction<Double> F = x ->
                    (coeff[2] * coeff[2] / 5.0) * Math.pow(x, 5)
                            + (coeff[2] * coeff[1] / 2.0) * Math.pow(x, 4)
                            + ((coeff[1] * coeff[1] + 2 * coeff[2] * coeff[0]) / 3.0) * Math.pow(x, 3)
                            + (coeff[1] * coeff[0]) * Math.pow(x, 2)
                            + (coeff[0] * coeff[0]) * x;

            if (low < up) {
                return F.apply(up) - F.apply(low);
            } else {
                return 0;
            }
        }

        private double getQuadraticProductIntegral(double[] x, double[] y, double low, double up) {
            double a0 = x[0]*y[0];
            double a1 = x[0]*y[1] + x[1]*y[0];
            double a2 = x[0]*y[2] + x[1]*y[1] + x[2]*y[0];
            double a3 = x[1]*y[2] + x[2]*y[1];
            double a4 = x[2]*y[2];

            java.util.function.DoubleFunction<Double> F = (t) ->
                    a0 * t +
                            a1 * Math.pow(t, 2) / 2.0 +
                            a2 * Math.pow(t, 3) / 3.0 +
                            a3 * Math.pow(t, 4) / 4.0 +
                            a4 * Math.pow(t, 5) / 5.0;

            if (TEST) {

                java.util.function.DoubleFunction<Double> FF = (t) -> evaluatePolynomialIntegralEndPt(t,
                        a0, a1, a2, a3, a4);

                double f1 = F.apply(up); double f0 = F.apply(low);
                double ff1 = FF.apply(up); double ff0 = FF.apply(low);
                System.err.println(f1 + " " + f0);
                System.err.println(ff1 + " " + ff0);
                System.err.println("test!");
                System.exit(-1);
            }

            if (low < up) {
                return F.apply(up) - F.apply(low);
            } else {
                return 0;
            }

        }
    }


    private class Cubic implements SplinesIntegral{
        private double getCubicProductIntegral(double[] x, double[] y, double low, double up) {

            double c0 = x[0] * y[0];

            double c1 = x[0] * y[1] + x[1] * y[0];

            double c2 = x[0] * y[2] + x[1] * y[1] + x[2] * y[0];

            double c3 = x[0] * y[3] + x[1] * y[2] + x[2] * y[1] + x[3] * y[0];

            double c4 = x[1] * y[3] + x[2] * y[2] + x[3] * y[1];

            double c5 = x[2] * y[3] + x[3] * y[2];

            double c6 = x[3] * y[3];

            java.util.function.DoubleFunction<Double> F = (t) ->
                    c0 * t +
                            c1 * Math.pow(t, 2) / 2.0 +
                            c2 * Math.pow(t, 3) / 3.0 +
                            c3 * Math.pow(t, 4) / 4.0 +
                            c4 * Math.pow(t, 5) / 5.0 +
                            c5 * Math.pow(t, 6) / 6.0 +
                            c6 * Math.pow(t, 7) / 7.0;

            if (low < up) {
                return F.apply(up) - F.apply(low);
            } else {
                return 0.0;
            }
        }

        private double[] getCubic1(int i){
            double[] a = new double[3];
            double[] cubicCoeff = new double[4];
            if (expandedKnots[i + 3] != expandedKnots[i] && expandedKnots[i + 2] != expandedKnots[i] &&
                    expandedKnots[i + 2] != expandedKnots[i + 1]) {
                a[0] = 1/((expandedKnots[i + 3] - expandedKnots[i]) * (expandedKnots[i + 2] - expandedKnots[i]) *
                        (expandedKnots[i + 2] - expandedKnots[i + 1]));
            } else {
                a[0] = 0;
            }
            if (expandedKnots[i + 3] != expandedKnots[i + 1] && expandedKnots[i + 3] != expandedKnots[i] &&
                    expandedKnots[i + 2] != expandedKnots[i + 1]) {
                a[1] = 1/((expandedKnots[i + 3] - expandedKnots[i + 1]) * (expandedKnots[i + 3] - expandedKnots[i]) *
                        (expandedKnots[i + 2] - expandedKnots[i + 1]));
            } else {
                a[1] = 0;
            }
            if (expandedKnots[i + 4] != expandedKnots[i + 1] && expandedKnots[i + 3] != expandedKnots[i + 1] &&
                    expandedKnots[i + 2] != expandedKnots[i + 1]) {
                a[2] = 1/((expandedKnots[i + 4] - expandedKnots[i + 1]) * (expandedKnots[i + 3] - expandedKnots[i + 1]) *
                        (expandedKnots[i + 2] - expandedKnots[i + 1]));
            } else {
                a[2] = 0;
            }
            cubicCoeff[3] = -(a[0] + a[1] + a[2]);
            cubicCoeff[2] = 2 * a[0] * expandedKnots[i] + a[1] * expandedKnots[i] + a[1] * expandedKnots[i + 1] +
                    2 * a[2] * expandedKnots[i + 1] + a[0] * expandedKnots[i + 2] + a[1] * expandedKnots[i + 3] +
                    a[2] * expandedKnots[i + 4];
            cubicCoeff[1] = -(a[0] * expandedKnots[i] * expandedKnots[i] + a[2] * expandedKnots[i + 1] * expandedKnots[i + 1]
                    + a[1] * expandedKnots[i] * expandedKnots[i + 1] + 2 * a[0] * expandedKnots[i] * expandedKnots[i + 2]
                    + a[1] * expandedKnots[i] * expandedKnots[i + 3] + a[1] * expandedKnots[i + 1] * expandedKnots[i + 3]
                    + 2 * a[2] * expandedKnots[i + 1] * expandedKnots[i + 4]);
            cubicCoeff[0] = a[0] * expandedKnots[i] * expandedKnots[i] * expandedKnots[i + 2] +
                    a[1] * expandedKnots[i] * expandedKnots[i + 1] * expandedKnots[i + 3] +
                    a[2] * expandedKnots[i + 1] * expandedKnots[i + 1] * expandedKnots[i + 4];
            return  cubicCoeff;
        }


        private double[] getCubic2(int i) {


            double[] cubicCoeff = new double[4];
            double c;

            if (expandedKnots[i + 3] != expandedKnots[i] &&
                    expandedKnots[i + 2] != expandedKnots[i] &&
                    expandedKnots[i + 1] != expandedKnots[i]) {

                c = 1/((expandedKnots[i + 3] - expandedKnots[i]) *
                        (expandedKnots[i + 2] - expandedKnots[i]) *
                        (expandedKnots[i + 1] - expandedKnots[i]));

                cubicCoeff[0] = -Math.pow(expandedKnots[i], 3) * c;
                cubicCoeff[1] = 3 * Math.pow(expandedKnots[i], 2) * c;
                cubicCoeff[2] = -3 * expandedKnots[i] * c;
                cubicCoeff[3] = c;
            }

            return cubicCoeff;
        }

        private double[] getCubic3(int i){
            double[] a = new double[3];
            double[] cubicCoeff = new double[4];
            if (expandedKnots[i + 3] != expandedKnots[i + 1] && expandedKnots[i + 3] != expandedKnots[i] &&
                    expandedKnots[i + 3] != expandedKnots[i + 2]) {
                a[0] = 1/((expandedKnots[i + 3] - expandedKnots[i + 1]) * (expandedKnots[i + 3] - expandedKnots[i]) *
                        (expandedKnots[i + 3] - expandedKnots[i + 2]));
            } else {
                a[0] = 0;
            }
            if (expandedKnots[i + 4] != expandedKnots[i + 1] && expandedKnots[i + 3] != expandedKnots[i + 1] &&
                    expandedKnots[i + 3] != expandedKnots[i + 2]) {
                a[1] = 1/((expandedKnots[i + 4] - expandedKnots[i + 1]) * (expandedKnots[i + 3] - expandedKnots[i + 1]) *
                        (expandedKnots[i + 3] - expandedKnots[i + 2]));
            } else {
                a[1] = 0;
            }
            if (expandedKnots[i + 4] != expandedKnots[i + 1] && expandedKnots[i + 4] != expandedKnots[i + 2] &&
                    expandedKnots[i + 3] != expandedKnots[i + 2]) {
                a[2] = 1/((expandedKnots[i + 4] - expandedKnots[i + 1]) * (expandedKnots[i + 4] - expandedKnots[i + 2]) *
                        (expandedKnots[i + 3] - expandedKnots[i + 2]));
            } else {
                a[2] = 0;
            }
            cubicCoeff[3] = a[0] + a[1] + a[2];
            cubicCoeff[2] = -(a[0] * expandedKnots[i] + a[1] * expandedKnots[i + 1] + a[2] * expandedKnots[i + 2] +
                    2 * a[0] * expandedKnots[i + 3] + a[1] * expandedKnots[i + 3] + a[1] * expandedKnots[i + 4] +
                    2 * a[2] * expandedKnots[i + 4]);
            cubicCoeff[1] = a[0] * expandedKnots[i + 3] * expandedKnots[i + 3] + a[2] * expandedKnots[i + 4] * expandedKnots[i + 4]
                    + 2 * a[0] * expandedKnots[i] * expandedKnots[i + 3] + a[1] * expandedKnots[i + 1] * expandedKnots[i + 3]
                    + a[1] * expandedKnots[i + 1] * expandedKnots[i + 4] + 2 * a[2] * expandedKnots[i + 2] * expandedKnots[i + 4]
                    + a[1] * expandedKnots[i + 3] * expandedKnots[i + 4];
            cubicCoeff[0] = -(a[0] * expandedKnots[i] * expandedKnots[i + 3] * expandedKnots[i + 3] +
                    a[2] * expandedKnots[i + 2] * expandedKnots[i + 4] * expandedKnots[i + 4] +
                    a[1] * expandedKnots[i + 1] * expandedKnots[i + 3] * expandedKnots[i + 4]);
            return  cubicCoeff;
        }

        private double[] getCubic4(int i) {


            double[] cubicCoeff = new double[4];
            double c;

            if (expandedKnots[i + 4] != expandedKnots[i + 1] &&
                    expandedKnots[i + 4] != expandedKnots[i + 2] &&
                    expandedKnots[i + 4] != expandedKnots[i + 3]) {

                c = 1/((expandedKnots[i + 4] - expandedKnots[i + 1]) *
                        (expandedKnots[i + 4] - expandedKnots[i + 2]) *
                        (expandedKnots[i + 4] - expandedKnots[i + 3]));

                cubicCoeff[0] = Math.pow(expandedKnots[i + 4], 3) * c;
                cubicCoeff[1] = - 3 * Math.pow(expandedKnots[i + 4], 2) * c;
                cubicCoeff[2] = 3 * expandedKnots[i + 4] * c;
                cubicCoeff[3] = -c;
            }

            return cubicCoeff;
        }



        public double getDiagonal(int i, double start, double end) {

            double value = 0;
            double low;
            double up;


            low = Math.max(start, expandedKnots[i]);
            up = Math.min(end, expandedKnots[i + 1]);
            if (expandedKnots[i + 3] != expandedKnots[i] && expandedKnots[i + 2] != expandedKnots[i] &&
                    expandedKnots[i + 1] != expandedKnots[i] && low < up) {
                value += ((Math.pow(up - expandedKnots[i], 7) - Math.pow(low - expandedKnots[i], 7))/
                        (7 * Math.pow((expandedKnots[i + 3] - expandedKnots[i]) * (expandedKnots[i + 2] - expandedKnots[i])
                                * (expandedKnots[i + 1] - expandedKnots[i]), 2)));
            }


            low = Math.max(start, expandedKnots[i + 1]);
            up = Math.min(end, expandedKnots[i + 2]);

            if (low < up) {
                value += getCubicProductIntegral(getCubic1(i), getCubic1(i), low, up);
            }


            low = Math.max(start, expandedKnots[i + 2]);
            up = Math.min(end, expandedKnots[i + 3]);
            if (low < up) {
                value += getCubicProductIntegral(getCubic3(i), getCubic3(i), low, up);
            }

            low = Math.max(start, expandedKnots[i + 3]);
            up = Math.min(end, expandedKnots[i + 4]);
            if (expandedKnots[i + 4] != expandedKnots[i + 1] && expandedKnots[i + 4] != expandedKnots[i + 2] &&
                    expandedKnots[i + 4] != expandedKnots[i + 3] && low < up) {
                value += ((Math.pow(up - expandedKnots[i + 4], 7) - Math.pow(low - expandedKnots[i + 4], 7))/
                        (7 * Math.pow((expandedKnots[i + 4] - expandedKnots[i + 1])
                                * (expandedKnots[i + 4] - expandedKnots[i + 2])
                                * (expandedKnots[i + 4] - expandedKnots[i + 3]), 2)));
            }

            return value;
        }

        public double getOffDiagonalOrder1(int i, double start, double end) {
            double value = 0;
            double low;
            double up;
            low = Math.max(start, expandedKnots[i + 1]);
            up = Math.min(end, expandedKnots[i + 2]);
            if (low < up) {
                value += getCubicProductIntegral(getCubic2(i + 1), getCubic1(i), low, up);
            }


            low = Math.max(start, expandedKnots[i + 2]);
            up = Math.min(end, expandedKnots[i + 3]);
            if (low < up) {
                value += getCubicProductIntegral(getCubic1(i + 1), getCubic3(i), low, up);
            }

            low = Math.max(start, expandedKnots[i + 3]);
            up = Math.min(end, expandedKnots[i + 4]);
            if (low < up) {
                value += getCubicProductIntegral(getCubic3(i + 1), getCubic4(i), low, up);
            }

            return value;
        }

        public double getOffDiagonalOrder2(int i, double start, double end) {
            double value = 0;
            double low;
            double up;
            low = Math.max(start, expandedKnots[i + 2]);
            up = Math.min(end, expandedKnots[i + 3]);
            if (low < up) {
                value += getCubicProductIntegral(getCubic2(i + 2), getCubic3(i), low, up);
            }


            low = Math.max(start, expandedKnots[i + 3]);
            up = Math.min(end, expandedKnots[i + 4]);
            if (low < up) {
                value += getCubicProductIntegral(getCubic1(i + 2), getCubic4(i), low, up);
            }

            return value;
        }

        public double getOffDiagonalOrder3(int i, double start, double end) {
            double value = 0;
            double low;
            double up;
            low = Math.max(start, expandedKnots[i + 3]);
            up = Math.min(end, expandedKnots[i + 4]);
            if (low < up) {
                value += getCubicProductIntegral(getCubic2(i + 3), getCubic4(i), low, up);
            }

            return value;
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

    public double[][] getIntegratedSquaredBasisMatrix(double start, double end) {
        int dim = coefficient.length;
//        double[][] mat = new double[dim][dim];
        double[][] mat = this.temp;
        for (int i = 0; i < dim ; i++) {
            mat[i][i] = splinesIntegral.getDiagonal(i, start, end);
        }
        for (int i = 0; i < dim - 1; i++) {
            double order1 = splinesIntegral.getOffDiagonalOrder1(i, start, end);
            mat[i][i + 1] = order1;
            mat[i + 1][i] = order1;
        }
        for (int i = 0; i < dim - 2; i++) {
            double order2 = splinesIntegral.getOffDiagonalOrder2(i, start, end);
            mat[i][i + 2] = order2;
            mat[i + 2][i] = order2;
        }
        for (int i = 0; i < dim - 3; i++) {
            double order3 = splinesIntegral.getOffDiagonalOrder3(i, start, end);
            mat[i][i + 3] = order3;
            mat[i + 3][i] = order3;
        }
        return mat;
    }

    public double getIntegral (double start, double end){
        getExpandedKnots();
        int dim = coefficient.length;
        double[][] mat = getIntegratedSquaredBasisMatrix(start, end);
        double sum = 0;
        if (end > start) {
            for (int i = 0; i < dim; i++) {
                sum += mat[i][i] * coefficient[i] * coefficient[i];
            }
            for (int i = 0; i < dim - 1; i++) {
                sum += 2 * mat[i][i + 1] * coefficient[i] * coefficient[i + 1];
            }
            for (int i = 0; i < dim - 2; i++) {
                sum += 2 * mat[i][i + 2] * coefficient[i] * coefficient[i + 2];
            }
            for (int i = 0; i < dim - 3; i++) {
                sum += 2 * mat[i][i + 3] * coefficient[i] * coefficient[i + 3];
            }
        }

        return sum;
    }



    public double[] getGradient (double start, double end){
        getExpandedKnots();
        int dim = coefficient.length;
        double[][] mat = getIntegratedSquaredBasisMatrix(start, end);
        double[] gradient = new double[dim];
        if (end > start) {
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    gradient[i] += 2 * mat[i][j] * coefficient[j];
                }
            }
        }
        return gradient;
    }


    public double evaluateExpSpline(double x) {
        double sum = 0.0;
        getExpandedKnots();

        for (int i = 0; i < coefficient.length; i++) {
            sum += coefficient[i] * BSplines.getSplineBasis(i, degree, x, expandedKnots);
        }

        return Math.exp(sum);
    }

    public double integrateExpSpline(final double a, final double b)
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


}





