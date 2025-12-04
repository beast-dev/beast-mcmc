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

/**
 * @author Pratyusa Datta
 * @author Marc A. Suchard
 */

public class IntegratedSquaredCubicSplines {


    private final double[] coefficient;
    private final double[] knots;
    private final double[] expandedKnots;
    private final int degree;
    private boolean expandedKnotsKnown;
    public IntegratedSquaredCubicSplines(double[] coefficient,
                                        double[] knots,
                                        int degree) {


        this.coefficient = coefficient;
        this.knots = knots;
        this.expandedKnots = new double[knots.length + 2 * degree];
        this.degree = degree;
        this.expandedKnotsKnown = false;
    }


    public void getExpandedKnots() {
            if(!expandedKnotsKnown) {
                for (int i = 0; i < degree; i++) {
                    expandedKnots[i] = knots[0];
                    expandedKnots[degree + knots.length + i] = knots[knots.length - 1];
                }
                for (int i = 0; i < knots.length; i++) {
                    expandedKnots[degree + i] = knots[i];
                }
            }
    }


    public double getCubicProductIntegral(double[] x, double[] y, double low, double up) {

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

    public double[] getCubic1(int i){
        double[] a = new double[3];
        double[] cubicCoeff = new double[4];
        if (expandedKnots[i + 3] != expandedKnots[i] && expandedKnots[i + 2] != expandedKnots[i] &&
                expandedKnots[i + 2] != expandedKnots[i + 1]) {
            a[0] += 1/((expandedKnots[i + 3] - expandedKnots[i]) * (expandedKnots[i + 2] - expandedKnots[i]) *
                    (expandedKnots[i + 2] - expandedKnots[i + 1]));
        }
        if (expandedKnots[i + 3] != expandedKnots[i + 1] && expandedKnots[i + 3] != expandedKnots[i] &&
                expandedKnots[i + 2] != expandedKnots[i + 1]) {
            a[1] += 1/((expandedKnots[i + 3] - expandedKnots[i + 1]) * (expandedKnots[i + 3] - expandedKnots[i]) *
                    (expandedKnots[i + 2] - expandedKnots[i + 1]));
        }
        if (expandedKnots[i + 4] != expandedKnots[i + 1] && expandedKnots[i + 3] != expandedKnots[i + 1] &&
                expandedKnots[i + 2] != expandedKnots[i + 1]) {
            a[2] += 1/((expandedKnots[i + 4] - expandedKnots[i + 1]) * (expandedKnots[i + 2] - expandedKnots[i + 1]) *
                    (expandedKnots[i + 2] - expandedKnots[i + 1]));
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


    public double[] getCubic2(int i) {


        double[] cubicCoeff = new double[4];
        double c = 0;

        if (expandedKnots[i + 3] != expandedKnots[i] &&
                expandedKnots[i + 2] != expandedKnots[i] &&
                expandedKnots[i + 1] != expandedKnots[i]) {

            c += 1/((expandedKnots[i + 3] - expandedKnots[i]) *
                    (expandedKnots[i + 2] - expandedKnots[i]) *
                    (expandedKnots[i + 1] - expandedKnots[i]));

            cubicCoeff[0] = -Math.pow(expandedKnots[i], 3) * c;
            cubicCoeff[1] = 3 * Math.pow(expandedKnots[i], 2) * c;
            cubicCoeff[2] = -3 * expandedKnots[i] * c;
            cubicCoeff[3] = c;
        }

        return cubicCoeff;
    }

    public double[] getCubic3(int i){
        double[] a = new double[3];
        double[] cubicCoeff = new double[4];
        if (expandedKnots[i + 3] != expandedKnots[i + 1] && expandedKnots[i + 3] != expandedKnots[i] &&
                expandedKnots[i + 3] != expandedKnots[i + 2]) {
            a[0] += 1/((expandedKnots[i + 3] - expandedKnots[i + 1]) * (expandedKnots[i + 3] - expandedKnots[i]) *
                    (expandedKnots[i + 3] - expandedKnots[i + 2]));
        }
        if (expandedKnots[i + 4] != expandedKnots[i + 1] && expandedKnots[i + 3] != expandedKnots[i + 1] &&
                expandedKnots[i + 3] != expandedKnots[i + 2]) {
            a[1] += 1/((expandedKnots[i + 4] - expandedKnots[i + 1]) * (expandedKnots[i + 3] - expandedKnots[i + 1]) *
                    (expandedKnots[i + 3] - expandedKnots[i + 2]));
        }
        if (expandedKnots[i + 4] != expandedKnots[i + 1] && expandedKnots[i + 4] != expandedKnots[i + 2] &&
                expandedKnots[i + 3] != expandedKnots[i + 2]) {
            a[2] += 1/((expandedKnots[i + 4] - expandedKnots[i + 1]) * (expandedKnots[i + 4] - expandedKnots[i + 2]) *
                    (expandedKnots[i + 3] - expandedKnots[i + 2]));
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

    public double[] getCubic4(int i) {


        double[] cubicCoeff = new double[4];
        double c = 0;

        if (expandedKnots[i + 4] != expandedKnots[i + 1] &&
                expandedKnots[i + 4] != expandedKnots[i + 2] &&
                expandedKnots[i + 4] != expandedKnots[i + 3]) {

            c += 1/((expandedKnots[i + 4] - expandedKnots[i + 1]) *
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

    public double[][] getIntegratedSquaredBasisMatrix(double start, double end) {
        getExpandedKnots();
       int dim = coefficient.length;
       double[][] mat = new double[dim][dim];
       for (int i = 0; i < dim ; i++) {
           mat[i][i] = getDiagonal(i, start, end);
       }
       for (int i = 0; i < dim - 1; i++) {
           mat[i][i + 1] = getOffDiagonalOrder1(i, start, end);
           mat[i + 1][i] = getOffDiagonalOrder1(i, start, end);
       }
       for (int i = 0; i < dim - 2; i++) {
           mat[i][i + 2] = getOffDiagonalOrder2(i, start, end);
           mat[i + 2][i] = getOffDiagonalOrder2(i, start, end);
       }
        for (int i = 0; i < dim - 3; i++) {
            mat[i][i + 3] = getOffDiagonalOrder3(i, start, end);
            mat[i + 3][i] = getOffDiagonalOrder3(i, start, end);
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
                    gradient[i] += 2 * mat[i][j] * coefficient[i];
                }
            }
        }
        return gradient;
    }
}





