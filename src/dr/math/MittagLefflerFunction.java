/*
 * MittagLefflerFunction.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import cern.jet.stat.Gamma;

/**
 * @author Marc Suchard
 *         <p/>
 *         Evaluates numerical approximations to the Mittag-Leffler function E_{\alpha}(-x^{\alpha}).
 *         This function forms a basis for many solutions to fractional differential equations (i.e., Levy flights)
 *         <p/>
 *         Diethelm K, Ford NJ, Freed AD, Luchko Y (2005) Algorithms for the fractional calculus: a selection
 *         of numerical methods. Comput. Methods Appl. Mech. Engrg. 194, 742-773.
 *         <p/>
 *         In my hands, the approximation from DKKL is poor for \alpha < 0.1 in the truncated series span of x
 *
 */
public class MittagLefflerFunction {

    private static final double EPS = 10E-15;
    private static final boolean DEBUG = false;

    public MittagLefflerFunction(double alpha) {
        this.alpha = alpha;
        this.beta = 1;
    }

    public MittagLefflerFunction(double alpha, double beta) {
        this.alpha = alpha;
        this.beta = beta;
    }

    public static void main(String arg[]) {

        System.out.println("Comparing power-series expansion to Pade approximation of");
        System.out.println("E_{\\alpha}(-x^{\\alpha})");

        int lenAlpha = 20;
        int lenX = 1000;

        double startAlpha = 0.05;
        double endAlpha = 0.999;

        double startX = 0.001;
        double endX = 20;

        double incAlpha = (endAlpha - startAlpha) / (lenAlpha);
        double incX = (endX - startX) / (lenX);

        double[][] padeResults = new double[lenAlpha][lenX];
        double[][] powsResults = new double[lenAlpha][lenX];

        long start, end;
        double alpha;
        double x;


        start = System.nanoTime();
        alpha = startAlpha;
        for(int i=0; i<lenAlpha; i++) {
            x = startX;
            for(int j=0; j<lenX; j++) {                
                padeResults[i][j] = evaluatePade(x,alpha);
                x += incX;
            }
            alpha += incAlpha;
        }
        end = System.nanoTime();
        System.out.println("Pade run: "+(end-start));

        start = System.nanoTime();
        alpha = startAlpha;       
        for(int i=0; i<lenAlpha; i++) {
            x = startX;
            for(int j=0; j<lenX; j++) {
                powsResults[i][j] = evaluatePowerSeries(x,alpha);

                x += incX;
            }
            alpha += incAlpha;
        }
        end = System.nanoTime();
        System.out.println("PowS run: "+(end-start));

        // Report MSE by alpha
        double MSE = 0;
        for(int i=0; i<lenAlpha; i++) {
            for(int j=0; j<lenX; j++) {
                final double delta = padeResults[i][j] - powsResults[i][j];
                MSE += delta * delta;
                if (Math.abs(delta) > 1E-1) {
                    System.err.println("Big difference at ("+i+","+j+") = "+delta);
                }

            }
        }
        MSE /= (lenAlpha*lenX);

        System.out.println("MSE = "+MSE);

    }

    public double evaluate(double x) {
        return evaluatePowerSeries(x, alpha, beta);
    }

    public static double evaluatePowerSeries(double x, double alpha, double beta) {
        return evaluatePowerSeries(x, alpha, beta, EPS);
    }

    public static double evaluatePowerSeries(double x, double alpha) {
        return evaluatePowerSeries(x,alpha,1.0,EPS);
    }

    public static double evaluatePowerSeries(double x, double alpha, double beta, double eps) {
        final double z = -Math.pow(x,alpha);
        double E = 0;
        double summand = 1;
        int k = 0;
        double XpowK = 1;
        while (Math.abs(summand) >= eps) {
            summand = XpowK / gamma(beta + alpha * k);
            XpowK *= z;
            k++;
            E += summand;
        }
        if (DEBUG) {
            System.err.println("Loops = " + k);
        }
        return E;
    }

    public static double evaluatePade(double x, double alpha) {
        return evaluatePade(x, 0, alpha);
    }

    public static double evaluatePade(double x, double y, double alpha) {

        if (alpha < 0.01 || alpha > 1) {
            System.err.println("alpha = "+alpha);
            throw new RuntimeException("Mittag-Leffler function only implemented for 0.01 <= alpha <= 1");
        }

        if (y != 0) {
            throw new RuntimeException("Mittag-Leffler function only implemented for real arguments");
        }

        double E;

        if (alpha == 1) { // Returns standard exponential function
            return Math.exp(-Math.pow(x, alpha));
        }

        if (x < 0) {

            throw new RuntimeException("Mittag-Leffler function only defined for positive quantities");

        } else if (x <= 0.1) {
            E = 0;
            final double z = -Math.pow(x, alpha);
            for (int k = 0; k <= 4; k++) {
                E += Math.pow(z, k) / gamma(1 + alpha * k);
            }

        } else if (x < 15) {

            // Do a linear interpolation between grid pts.
            // Pade approximation grid pts { 0.01, 0.02, ..., 0.99 }
            final double rescaledAlpha = alpha * 100.0;
            final int end = (int) rescaledAlpha;
            final int start = end - 1;
            double delta = rescaledAlpha - end;

            E = padeApproximation(x, start);
            if (delta > 0) {
                E *= (1 - delta);
                if (end < 99) { // Only 98 values.
                    E += padeApproximation(x, end) * delta;
                } else {
                    E += Math.exp(x) * delta;
                }
            }

        } else {

            E = 0;
            final double z = -Math.pow(x, -alpha);
            for (int k = 1; k <= 4; k++) {
                E -= Math.pow(z, k) / gamma(1 - alpha * k);
            }
        }

        return E;
    }

    private static double gamma(double z) {

        // To extend the gamma function to negative (non-integer) numbers, apply the relationship
        // \Gamma(z) = \frac{ \Gamma(z+n) }{ z(z+1)\cdots(z+n-1),
        // by choosing n such that z+n is positive
        if (z > 0) {
            return Gamma.gamma(z);
        }
        int n = ((int) -z) + 1;
        double factor = 1;
        for(int i=0; i<n; i++) {
            factor *= (z+i);
        }
        return Gamma.gamma(z+n) / factor;
    }


    private static double padeApproximation(double x, int i) {
        return (a0[i] + a1[i] * x + a2[i] * x * x) / (1.0 + b1[i] * x + b2[i] * x * x + b3[i] * x * x * x);
    }

    private static double[] a0 = {
            -4.24129e+02, -2.91078e+03, -1.59594e+03, -4.01668e+03, -4.80619e+03, -3.58700e+03, -4.31271e+03, -1.04813e+04, -1.13182e+04,
            -1.00534e+04, -7.71890e+03, -2.09949e+03, -2.01966e+04, -1.44185e+04, -3.58699e+03, -3.67384e+04, -8.34348e+03, -6.65920e+03,
            -8.80313e+03, -5.16937e+03, -1.81654e+04, -2.87739e+04, -2.83969e+04, -2.18950e+04, -1.84865e+04, -1.54930e+04, -2.52656e+04,
            -3.48286e+04, -2.23178e+04, -2.62152e+04, -1.10316e+04, -3.44887e+04, -2.53440e+04, -3.20037e+04, -7.88648e+04, -2.89781e+04,
            -1.42051e+04, -8.16551e+04, -5.75871e+04, -1.73808e+04, -1.03477e+04, -1.95086e+04, -3.33165e+04, -4.73286e+04, -1.47598e+04,
            -1.52914e+04, -1.67037e+04, -1.42611e+04, -6.54706e+04, -1.88800e+02, -5.16448e+04, -3.92154e+04, -1.35449e+04, -9.53841e+03,
            -7.93386e+02, -5.99232e+03, -3.37393e+04, -4.63960e+04, -1.51791e+03, -5.02729e+03, -4.80103e+04, -6.71512e+03, -3.12780e+04,
            -5.67954e+04, -1.49230e+03, -4.33396e+04, -7.52321e+03, -5.03149e+04, -1.13475e+04, -2.77604e+02, -5.51694e+04, -3.10305e+04,
            -1.60766e+03, +3.41935e+03, +1.73763e+03, -1.98004e+03, -2.55574e+03, -1.07310e+03, -4.94450e+03, -2.18066e+03, -5.34378e+03,
            -6.75467e+03, -7.02848e+02, -6.09980e+03, -1.03053e+04, -7.74771e+03, -8.38904e+03, -4.47750e+03, -3.90787e+03, -7.24688e+03,
            -1.52262e+04, -1.50108e+04, -1.06495e+04, -1.53906e+04, -6.00332e+03, -4.45797e+04, -3.55171e+04, -1.45397e+04, -1.05996e+04
    };

    private static double[] a1 = {
            -6.95856e+05, -2.31705e+06, -8.83786e+05, -1.51329e+06, -1.45550e+06, -9.53746e+05, -9.58824e+05, -2.20562e+06, -2.04162e+06,
            -1.65546e+06, -1.16952e+06, -3.31485e+05, -2.64628e+06, -1.77380e+06, -4.14163e+05, -3.87982e+06, -8.77610e+05, -6.61019e+05,
            -8.39225e+05, -4.63982e+05, -1.59929e+06, -2.41617e+06, -2.33523e+06, -1.74274e+06, -1.42847e+06, -1.16572e+06, -1.87367e+06,
            -2.50157e+06, -1.58819e+06, -1.80630e+06, -7.41001e+05, -2.31540e+06, -1.65943e+06, -2.06983e+06, -5.06620e+06, -1.95379e+06,
            -9.71770e+05, -5.23578e+06, -3.48533e+06, -1.17277e+06, -6.86777e+05, -1.29617e+06, -2.21252e+06, -3.11742e+06, -9.82326e+05,
            -1.04712e+06, -1.14533e+06, -9.91368e+05, -4.55618e+06, -1.35370e+04, -3.73170e+06, -2.89687e+06, -1.06392e+06, -7.48349e+05,
            -6.39333e+04, -4.43865e+05, -2.86190e+06, -4.05168e+06, -1.25852e+05, -4.64569e+05, -4.61103e+06, -7.08363e+05, -3.28986e+06,
            -6.32239e+06, -1.82752e+05, -5.49988e+06, -1.09682e+06, -7.82828e+06, -2.04800e+06, -5.35397e+04, -1.58192e+07, -1.02892e+07,
            -4.49254e+06, -2.41906e+06, -1.03357e+06, -1.70209e+05, -3.99973e+05, -1.02865e+05, -4.63767e+05, -2.19047e+05, -5.39996e+05,
            -6.88012e+05, -6.90578e+04, -5.83751e+05, -9.15547e+05, -6.67011e+05, -6.93431e+05, -3.59964e+05, -2.96304e+05, -5.09096e+05,
            -1.00158e+06, -9.42028e+05, -6.26399e+05, -8.47297e+05, -3.11566e+05, -2.14492e+06, -1.59674e+06, -6.10793e+05, -4.16113e+05
    };

    private static double[] a2 = {
            -3.91944e+05, -1.21393e+06, -4.85177e+05, -6.78102e+05, -6.43881e+05, -4.52754e+05, -4.24475e+05, -1.09649e+06, -9.31135e+05,
            -7.53470e+05, -5.30686e+05, -1.78194e+05, -1.18361e+06, -7.87251e+05, -1.80624e+05, -1.57708e+06, -3.81391e+05, -2.79508e+05,
            -3.52670e+05, -1.86660e+05, -6.55973e+05, -9.62412e+05, -9.35683e+05, -6.87682e+05, -5.55414e+05, -4.48241e+05, -7.20787e+05,
            -9.38597e+05, -5.97404e+05, -6.59350e+05, -2.63517e+05, -8.31618e+05, -5.79558e+05, -7.12936e+05, -1.72863e+06, -7.07274e+05,
            -3.54821e+05, -1.75872e+06, -1.07617e+06, -4.08516e+05, -2.31234e+05, -4.29282e+05, -7.18492e+05, -9.81295e+05, -3.04968e+05,
            -3.26140e+05, -3.46912e+05, -2.93586e+05, -1.30849e+06, -3.83675e+03, -1.02716e+06, -7.80114e+05, -2.89899e+05, -1.93214e+05,
            -1.59875e+04, -9.52614e+04, -6.63678e+05, -8.99289e+05, -2.41507e+04, -9.19187e+04, -8.60143e+05, -1.31044e+05, -5.33551e+05,
            -9.45394e+05, -2.59499e+04, -6.64464e+05, -1.26030e+05, -7.16973e+05, -1.58167e+05, -2.26639e+03, -6.73086e+05, +6.63348e+04,
            -1.16371e+05, -6.69717e+04, -5.13330e+05, -5.46410e+05, -4.08010e+05, -1.04591e+05, -3.57502e+05, -1.28115e+05, -2.58539e+05,
            -2.75764e+05, -2.40789e+04, -1.76849e+05, -2.47199e+05, -1.56824e+05, -1.41730e+05, -6.29983e+04, -4.46099e+04, -6.54940e+04,
            -1.06340e+05, -7.92810e+04, -3.98859e+04, -3.71308e+04, +7.56057e+03, -1.19855e+04, +2.06708e+04, +1.90402e+04, +2.04340e+04
    };

    private static double[] b1 = {
            -1.38637e+06, -4.60100e+06, -1.74711e+06, -2.98958e+06, -2.86638e+06, -1.86815e+06, -1.87447e+06, -4.27837e+06, -3.95568e+06,
            -3.19371e+06, -2.24658e+06, -6.28014e+05, -5.03944e+06, -3.36281e+06, -7.82040e+05, -7.32069e+06, -1.63991e+06, -1.23090e+06,
            -1.55465e+06, -8.57478e+05, -2.93380e+06, -4.41685e+06, -4.23970e+06, -3.14803e+06, -2.56670e+06, -2.08289e+06, -3.32362e+06,
            -4.41764e+06, -2.78265e+06, -3.15184e+06, -1.28658e+06, -3.98374e+06, -2.84133e+06, -3.51975e+06, -8.54865e+06, -3.24061e+06,
            -1.59385e+06, -8.60767e+06, -5.74012e+06, -1.87857e+06, -1.09392e+06, -2.04698e+06, -3.46552e+06, -4.84992e+06, -1.51324e+06,
            -1.59246e+06, -1.72774e+06, -1.48061e+06, -6.75094e+06, -1.98223e+04, -5.41357e+06, -4.15703e+06, -1.50142e+06, -1.04887e+06,
            -8.86439e+04, -6.19108e+05, -3.88535e+06, -5.44254e+06, -1.69238e+05, -6.11566e+05, -6.00297e+06, -9.06252e+05, -4.18611e+06,
            -7.94766e+06, -2.26154e+05, -6.74519e+06, -1.32094e+06, -9.32883e+06, -2.40096e+06, -6.21910e+04, -1.78833e+07, -1.15056e+07,
            -4.82429e+06, -2.52295e+06, -1.06853e+06, -2.12224e+05, -4.61840e+05, -1.25079e+05, -5.63528e+05, -2.62855e+05, -6.45497e+05,
            -8.19312e+05, -8.24283e+04, -6.97471e+05, -1.10278e+06, -8.05409e+05, -8.40712e+05, -4.37447e+05, -3.62820e+05, -6.30388e+05,
            -1.25254e+06, -1.18669e+06, -7.97901e+05, -1.09242e+06, -4.06277e+05, -2.84119e+06, -2.14673e+06, -8.34188e+05, -5.77797e+05
    };

    private static double[] b2 = {
            -7.95778e+05, -2.50385e+06, -1.01562e+06, -1.44743e+06, -1.39751e+06, -9.96694e+05, -9.52144e+05, -2.48747e+06, -2.15524e+06,
            -1.77305e+06, -1.26954e+06, -4.28531e+05, -2.92906e+06, -1.98184e+06, -4.63020e+05, -4.14017e+06, -1.01120e+06, -7.55801e+05,
            -9.70538e+05, -5.25412e+05, -1.87316e+06, -2.80701e+06, -2.77398e+06, -2.07911e+06, -1.71282e+06, -1.40906e+06, -2.30555e+06,
            -3.06951e+06, -1.98723e+06, -2.24602e+06, -9.19175e+05, -2.94522e+06, -2.10336e+06, -2.64342e+06, -6.54223e+06, -2.68137e+06,
            -1.36569e+06, -7.05295e+06, -4.50994e+06, -1.67833e+06, -9.77034e+05, -1.85582e+06, -3.18194e+06, -4.47200e+06, -1.42189e+06,
            -1.54553e+06, -1.69118e+06, -1.47062e+06, -6.75481e+06, -2.02768e+04, -5.59798e+06, -4.37101e+06, -1.64334e+06, -1.14460e+06,
            -9.79799e+04, -6.46389e+05, -4.38767e+06, -6.20515e+06, -1.85486e+05, -7.05142e+05, -6.97998e+06, -1.08843e+06, -4.94518e+06,
            -9.47107e+06, -2.75373e+05, -8.15645e+06, -1.64296e+06, -1.15421e+07, -3.01897e+06, -7.68192e+04, -2.32888e+07, -1.46317e+07,
            -6.82061e+06, -3.76790e+06, -2.12490e+06, -7.27938e+05, -8.89107e+05, -2.07143e+05, -8.02633e+05, -3.39154e+05, -7.66988e+05,
            -9.08437e+05, -8.47142e+04, -6.67838e+05, -9.63517e+05, -6.51177e+05, -6.24128e+05, -2.99094e+05, -2.22076e+05, -3.35377e+05,
            -5.74996e+05, -4.67003e+05, -2.55226e+05, -2.68144e+05, -7.10353e+04, -2.70455e+05, -4.23690e+04, +4.64083e+04, +7.55603e+04
    };

    private static double[] b3 = {
            -2.07400e+02, -1.27923e+03, -8.32687e+02, -1.38260e+03, -1.71118e+03, -1.60873e+03, -1.75012e+03, -5.94850e+03, -5.60193e+03,
            -5.27749e+03, -4.29168e+03, -1.83967e+03, -1.23451e+04, -9.23408e+03, -2.35432e+03, -2.18683e+04, -6.19007e+03, -4.95376e+03,
            -6.89880e+03, -3.91939e+03, -1.53968e+04, -2.44057e+04, -2.61404e+04, -2.08550e+04, -1.82554e+04, -1.59998e+04, -2.80571e+04,
            -3.92527e+04, -2.72201e+04, -3.21389e+04, -1.37449e+04, -4.73795e+04, -3.52888e+04, -4.67058e+04, -1.22026e+05, -5.55978e+04,
            -3.02386e+04, -1.56056e+05, -9.91720e+04, -4.29558e+04, -2.58727e+04, -5.14821e+04, -9.22518e+04, -1.34275e+05, -4.47789e+04,
            -5.17302e+04, -5.87605e+04, -5.31760e+04, -2.53173e+05, -7.97517e+02, -2.28252e+05, -1.85908e+05, -7.49106e+04, -5.30462e+04,
            -4.69973e+03, -2.88690e+04, -2.23581e+05, -3.24201e+05, -9.04156e+03, -3.77940e+04, -3.78472e+05, -6.26832e+04, -2.68704e+05,
            -5.09518e+05, -1.51317e+04, -4.07908e+05, -8.43703e+04, -5.00161e+05, -1.17538e+05, -1.50446e+03, -5.32630e+05, +2.02771e+05,
            -1.11373e+05, -8.08538e+04, -7.87594e+05, -9.06762e+05, -7.36776e+05, -2.05284e+05, -7.67703e+05, -3.02753e+05, -6.75138e+05,
            -8.00006e+05, -7.79293e+04, -6.42815e+05, -1.01397e+06, -7.33497e+05, -7.63282e+05, -3.96021e+05, -3.31368e+05, -5.84487e+05,
            -1.17309e+06, -1.12366e+06, -7.66267e+05, -1.06481e+06, -4.01715e+05, -2.85641e+06, -2.19239e+06, -8.65369e+05, -6.08788e+05
    };

    private double alpha;
    private double beta;

}
