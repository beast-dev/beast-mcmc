/*
 * MMAlgorithm.java
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

package dr.inference.multidimensionalscaling.mm;

/**
 * Created by msuchard on 12/15/15.
 */
public abstract class MMAlgorithm {

    public static final double DEFAULT_TOLERANCE = 1E-1;
    public static final int DEFAULT_MAX_ITERATIONS = 1000;

    public double[] findMode(final double[] startingValue) throws NotConvergedException {
        return findMode(startingValue, DEFAULT_TOLERANCE, DEFAULT_MAX_ITERATIONS);
    }

    public double[] findMode(final double[] startingValue, final double tolerance,
                             final int maxIterations) throws NotConvergedException {

        if (DEBUG) {
            System.err.println("Starting findMode with " + tolerance + " " + maxIterations);
        }

        double[] buffer1 = new double[startingValue.length];
        double[] buffer2 = new double[startingValue.length];

        double[] current = buffer1;
        double[] next = buffer2;

        System.arraycopy(startingValue, 0, next, 0, startingValue.length);
        int iteration = 0;

        do {
            // Move next -> current
            double[] tmp = current;
            current = next;
            next = tmp;

            if (DEBUG) {
                System.err.println("Current: " + printArray(current));
            }

            mmUpdate(current, next);
            ++iteration;

            if (DEBUG) {
                System.err.println("Finished iteration " + iteration);
            }

        } while (convergenceCriterion(next, current) > tolerance && iteration < maxIterations);

        if (iteration >= maxIterations) {
            throw new NotConvergedException();
        }

        if (DEBUG) {
            System.err.println("Final  : " + printArray(next));
        }

        return next;
    }

    static private String format =  "%5.3e";

    protected String printArray(double[] x) {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format(format, x[0]));
        for (int i = 1; i < x.length; ++i) {
            sb.append(", ").append(String.format(format, x[i]));
        }
        return sb.toString();
    }

    protected abstract void mmUpdate(final double[] current, double[] next);

    private double convergenceCriterion(final double[] current, final double[] previous) {
        double norm = 0.0;

        for (int i = 0; i < current.length; ++i) {
            norm += (current[i] - previous[i]) * (current[i] - previous[i]);
        }

        double value = Math.sqrt(norm);

        if (DEBUG) {
            System.err.println("Convergence = " + value);
        }

        return value;
    }


    class NotConvergedException extends Exception {
        // Nothing interesting
    }

    private static final boolean DEBUG = false;
}
