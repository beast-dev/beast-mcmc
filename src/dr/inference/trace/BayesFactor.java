/*
 * BayesFactor.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.inference.trace;

import static java.lang.Math.*;
import java.io.Reader;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;
import java.util.ArrayList;

import dr.math.MathUtils;
import dr.stats.DiscreteStatistics;

/**
 * Calculate the approximate Bayes Factor via importance sampling
 * (Newton & Raftery, 1994)
 * <p/>
 * This class is based on C++ code from Marc Suchard and Benjamin Redelings
 *
 * @author Alexei Drummond
 * @author Marc Suchard
 * @author Benjamin Redelings
 */
public class BayesFactor {

    static final double MAX_FLOAT = 3.40282347e+38F;
    static final double LOG_LIMIT = -MAX_FLOAT / 100;
    static final double LOG_0 = -MAX_FLOAT;

    static final double NATS = 40;

    public static double Pr_harmonic(double[] v) {
        double sum = 0;

        for (double aV : v) {
            sum += aV;
        }

        double denominator = LOG_0;

        for (double value : v) {
            denominator = logsum(denominator, sum - value);
        }

        return sum - denominator + log((double) v.length);
    }

    static double logsum(double x, double y) {
        double temp = y - x;
        if (temp > NATS || x < LOG_LIMIT)
            return y;
        else if (temp < -NATS || y < LOG_LIMIT)
            return x;
        else
            return (x + Math.log1p(Math.exp(temp)));
    }

    public static double Pr_smoothed(double[] v, double delta, double Pdata) {

        double log_delta = log(delta);
        double log_inv_delta = log(1 - delta);

        double n = v.length;
        double top = log(n) + log_delta - log_inv_delta + Pdata;
        double bottom = log(n) + log_delta - log_inv_delta;

        for (int i = 0; i < v.length; i++) {
            double weight = -logsum(log_delta, log_inv_delta + v[i] - Pdata);
            top = logsum(top, weight + v[i]);

            bottom = logsum(bottom, weight);
        }

        return top - bottom;
    }

    public static double Pr_smoothed(double[] v) {

        // Sample from the prior with probability delta...
        double delta = 0.01;

        // Use the harmonic estimator as the starting guess
        double Pdata = Pr_harmonic(v);

        // initialize this to a value which allows it to enter the loop
        double deltaP = 1.0;

        int iterations = 0;
        double dx;

        while (abs(deltaP) > 1.0e-3) {

            double g1 = Pr_smoothed(v, delta, Pdata) - Pdata;

            double Pdata2 = Pdata + g1;

            dx = g1 * 10;
            double g2 = Pr_smoothed(v, delta, Pdata + dx) - (Pdata + dx);

            double dgdx = (g2 - g1) / dx;

            double Pdata3 = Pdata - g1 / dgdx;
            if (Pdata3 < 2.0 * Pdata || Pdata3 > 0 || Pdata3 > 0.5 * Pdata) {
                Pdata3 = Pdata + 10 * g1;
            }

            double g3 = Pr_smoothed(v, delta, Pdata3) - Pdata3;

            // Try to do Newton's method
            if (abs(g3) <= abs(g2) && ((g3 > 0) || (abs(dgdx) > 0.01))) {
                deltaP = Pdata3 - Pdata;
                Pdata = Pdata3;
            }
            // Otherwise see if we can go 10 times as far as one step
            else if (abs(g2) <= abs(g1)) {
                Pdata2 += g2;
                deltaP = Pdata2 - Pdata;
                Pdata = Pdata2;
            }
            // Otherwise go one step
            else {
                deltaP = g1;
                Pdata += g1;
            }

            iterations++;

            // if we aren't converging, warn and don't give an answer
            if (iterations > 400) {
                System.err.println("Probabilities not converging!!!");
                return LOG_0;
            }
        }
        return Pdata;
    }

    static double[] bootstrap(double[] sample, int nresamples, int blocksize) {
        double[] results = new double[nresamples];

        for (int i = 0; i < nresamples; i++) {
            double[] temp = bootstrapSample(sample, blocksize);
            results[i] = Pr_smoothed(temp);
        }
        return results;
    }

    static double[] bootstrapSample(double[] sample, int blocksize) {
        assert(blocksize > 0);
        assert(blocksize <= sample.length);

        double[] resample = new double[sample.length];

        int i = 0;
        while (i < resample.length) {
            int j = MathUtils.nextInt(sample.length + 1 - blocksize);
            for (int k = 0; k < blocksize && i < resample.length; k++, i++, j++) {
                resample[i] = sample[j];
            }
        }

        return resample;
    }


    private static double[] readValues(Reader r) throws IOException {

        BufferedReader reader = new BufferedReader(r);

        List<Double> values = new ArrayList<Double>();
        String line = reader.readLine();
        while (line != null) {
            values.add(Double.parseDouble(line));
            line = reader.readLine();
        }

        double[] v = new double[values.size()];
        for (int i = 0; i < values.size(); i++) {
            v[i] = values.get(i);
        }

        return v;
    }

    public static void main(String[] argsString) throws IOException {

        double[] values = readValues(new InputStreamReader(System.in));

        double PM = Pr_smoothed(values);

        System.out.print("P(M|data) = " + PM + "  ");
        System.out.flush();

        //---------- Get bootstrap sample --------------/

        int blocksize = values.length / 100 + 2;
        if (blocksize > values.length) {
            blocksize = values.length;
        }

        double[] values2 = bootstrap(values, 100, blocksize);

        double stddev = DiscreteStatistics.stdev(values2);

        System.out.println("  +- " + stddev);
    }
}
