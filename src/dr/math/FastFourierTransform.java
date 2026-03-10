/*
 * FastFourierTransform.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.math;

/**
 * Performs FFT on vectors with lengths equaling powers-of-2
 *
 * @author Marc A. Suchard
 * @author Frederik M. Andersen
 */
public class FastFourierTransform {

    /**
     * Computes the fast fourier transform
     *
     * @param data    an interleaved array of (real,complex) values
     * @param nn      data length
     * @param inverse true is performing inverse FFT
     */
    public static void fft(double[] data, int nn, boolean inverse) {
        int n, mmax, m, j, istep, i;
        double wtemp, wr, wpr, wpi, wi, theta;
        double tempr, tempi;

        final double radians;
        if (inverse) {
            radians = 2.0 * Math.PI;
        } else {
            radians = -2.0 * Math.PI;
        }

        // reverse-binary reindexing
        n = nn << 1;
        j = 1;
        for (i = 1; i < n; i += 2) {
            if (j > i) {
                swap(data, j - 1, i - 1);
                swap(data, j, i);
            }
            m = nn;
            while (m >= 2 && j > m) {
                j -= m;
                m >>= 1;
            }
            j += m;
        }

        // here begins the Danielson-Lanczos section
        mmax = 2;
        while (n > mmax) {
            istep = mmax << 1;
            theta = radians / mmax;
            wtemp = Math.sin(0.5 * theta);
            wpr = -2.0 * wtemp * wtemp;
            wpi = Math.sin(theta);
            wr = 1.0;
            wi = 0.0;
            for (m = 1; m < mmax; m += 2) {
                for (i = m; i <= n; i += istep) {
                    j = i + mmax;
                    tempr = wr * data[j - 1] - wi * data[j];
                    tempi = wr * data[j] + wi * data[j - 1];

                    data[j - 1] = data[i - 1] - tempr;
                    data[j] = data[i] - tempi;
                    data[i - 1] += tempr;
                    data[i] += tempi;
                }
                wtemp = wr;
                wr += wr * wpr - wi * wpi;
                wi += wi * wpr + wtemp * wpi;
            }
            mmax = istep;
        }
    }


    /**
     * Computes the fast fourier transform
     *
     * @param ca    an array of complex values
     * @param inverse true if performing inverse FFT
     */
    public static void fft(ComplexArray ca, boolean inverse) {

        final double[] real = ca.real;
        final double[] complex = ca.complex;

        int n, mmax, m, j, istep, i;
        double wtemp, wr, wpr, wpi, wi, theta;
        double tempr, tempi;

        final double radians;
        if (inverse) {
            radians = 2.0 * Math.PI;
        } else {
            radians = -2.0 * Math.PI;
        }

        // reverse-binary reindexing
        n = ca.length << 1;
        j = 1;
        for (i = 1; i < n; i += 2) {
            if (j > i) {
                final int halfI = i >> 1;
                final int halfJ = j >> 1;
                swap(real, halfJ, halfI);
                swap(complex, halfJ, halfI);
            }
            m = ca.length;
            while (m >= 2 && j > m) {
                j -= m;
                m >>= 1;
            }
            j += m;
        }

        // here begins the Danielson-Lanczos section
        mmax = 2;
        while (n > mmax) {
            istep = mmax << 1;
            theta = (radians / mmax);
            wtemp = Math.sin(0.5 * theta);
            wpr = -2.0 * wtemp * wtemp;
            wpi = Math.sin(theta);
            wr = 1.0;
            wi = 0.0;
            for (m = 1; m < mmax; m += 2) {
                for (i = m; i <= n; i += istep) {
                    j = i + mmax;
                    final int halfI = i >> 1;
                    final int halfJ = j >> 1;
                    tempr = wr * real[halfJ] - wi * complex[halfJ];
                    tempi = wr * complex[halfJ] + wi * real[halfJ];

                    real[halfJ] = real[halfI] - tempr;
                    complex[halfJ] = complex[halfI] - tempi;
                    real[halfI] += tempr;
                    complex[halfI] += tempi;
                }
                wtemp = wr;
                wr += wr * wpr - wi * wpi;
                wi += wi * wpr + wtemp * wpi;
            }
            mmax = istep;
        }
    }

    /**
     * Helper function to compute fft on a real signal of length nn
     * using a complex fft of half the size, then unpacks (or vice versa for inverse).
     *
     * @param data    real input of length nn (forward) or interleaved spectrum (inverse)
     * @param nn      data length (must be a power of 2)
     * @param inverse true if performing inverse FFT
     */
    public static void rfft(double[] data, int nn, boolean inverse) {
        int halfN = nn >> 1; // equivalent to nn / 2

        // phase for the twiddle factors based on direction
        double theta = inverse ? (2.0 * Math.PI / nn) : (-2.0 * Math.PI / nn);

        if (!inverse) {
            // FORWARD: Treat real data as halfN complex pairs and do complex FFT first
            fft(data, halfN, false);
        } else {
            // INVERSE: Undo DC/Nyquist packing first
            double dc = data[0];
            double nyq = data[1];
            data[0] = 0.5 * (dc + nyq);
            data[1] = 0.5 * (dc - nyq);
        }

        // Shared Twiddle Factor Setup
        double wtemp = Math.sin(0.5 * theta);
        double wpr = -2.0 * wtemp * wtemp;
        double wpi = Math.sin(theta);
        double wr = 1.0 + wpr;
        double wi = wpi;

        // Shared Unpack/Pack Loop
        for (int k = 1; k < (halfN >> 1) + 1; k++) {
            int k2 = k << 1;
            int mk2 = (halfN - k) << 1;

            double zr = data[k2], zi = data[k2 + 1];
            double cr = data[mk2], ci = data[mk2 + 1];

            double sr = 0.5 * (zr + cr), si = 0.5 * (zi - ci);
            double dr = 0.5 * (zr - cr), di = 0.5 * (zi + ci);

            double tr = wr * di + wi * dr;
            double ti = -(wr * dr - wi * di);

            if (!inverse) {
                data[k2] = sr + tr;
                data[k2 + 1] = si + ti;
                data[mk2] = sr - tr;
                data[mk2 + 1] = -(si - ti);
            } else {
                data[k2] = sr - tr;
                data[k2 + 1] = si - ti;
                data[mk2] = sr + tr;
                data[mk2 + 1] = -(si + ti);
            }

            // Advance twiddle factors
            wtemp = wr;
            wr += wr * wpr - wi * wpi;
            wi += wi * wpr + wtemp * wpi;
        }

        if (!inverse) {
            // FORWARD: Handle DC and Nyquist last
            double dc = data[0];
            double nyq = data[1];
            data[0] = dc + nyq;
            data[1] = dc - nyq;
        } else {
            // INVERSE: Do inverse complex FFT of size halfN last, then scale
            fft(data, halfN, true);

            double scale = 1.0 / halfN;
            for (int i = 0; i < nn; i++) {
                data[i] *= scale;
            }
        }
    }

    private static void swap(double[] x, int i, int j) {
        double tmp = x[i];
        x[i] = x[j];
        x[j] = tmp;
    }
}
