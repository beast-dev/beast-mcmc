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

import dr.inference.multidimensionalscaling.MultiDimensionalScalingLikelihood;
import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.xml.*;

/**
 * Created by msuchard on 12/15/15.
 */
public abstract class MMAlgorithm {

    public static final double DEFAULT_TOLERANCE = 1E-1;
    public static final int DEFAULT_MAX_ITERATIONS = 1000;

    public double[] findMode(final double[] startingValue) throws NotConvergedException {
        return findMode(startingValue, DEFAULT_TOLERANCE, DEFAULT_MAX_ITERATIONS);
    }

    private void copyDifference(double[] dest, final double[] x, final double[] y) {
        final int len = dest.length;
        for (int i = 0; i < len; ++i) {
            dest[i] = x[i] - y[i];
        }
    }

    public double[] findMode(final double[] startingValue, final double tolerance,
                             final int maxIterations) throws NotConvergedException {

        if (DEBUG) {
            System.err.println("Starting findMode with " + tolerance + " " + maxIterations);
        }

        double[] buffer1 = new double[startingValue.length];
        double[] buffer2 = new double[startingValue.length];

        double[] previous = buffer1;
        double[] x = buffer2;

        System.arraycopy(startingValue, 0, x, 0, startingValue.length);
        iteration = 0;

        if (qnQ == 0) {
            // Standard MM
            do {
                // Move x -> previous
                double[] tmp = previous; previous = x; x = tmp;

                if (DEBUG) {
                    System.err.println("Current: " + printArray(previous));
                }

                mmUpdate(previous, x);
                ++iteration;

                if (DEBUG) {
                    System.err.println("Finished iteration " + iteration);
                }

            } while (convergenceCriterion(x, previous) > tolerance && iteration < maxIterations);

        } else {
            // Quasi-Newton acceleration
            final int J = startingValue.length;
            double[][] secantsU = new double[1][J];
            double[][] secantsV = new double[1][J];

            double[] Fx = new double[J];
            double[] C = new double[J];
            double[] x0 = new double[J];

            // Fill initial secants
            int countU = 0;
            int countV = 0;

            for (int q = 0; q < qnQ; ++q) {
                double[] tmp = previous; previous = x; x = tmp;
                mmUpdate(previous, x);
                ++iteration;

                if (countU == 0) {
                    copyDifference(secantsU[countU], x, previous);
                    ++countU;
                } else if (countV < qnQ - 1) {
                    copyDifference(secantsU[countU], x, previous);
                    System.arraycopy(secantsU[countU], 0, secantsV[countV], 0, J);
                    ++countU;
                    ++countV;
                } else {
                    copyDifference(secantsV[countV], x, previous);
                    ++countV;
                }
            }

            int newestSecant = qnQ - 1;
            int previousSecant = newestSecant - 1;

            boolean done = false;

            while (!done) {

                System.arraycopy(x, 0, x0, 0, J);

                // 2 cycles for each QN step
                double[] tmp = previous; previous = x; x = tmp;
                mmUpdate(previous, x);
                ++iteration;

                copyDifference(secantsU[newestSecant], x, previous);
                System.arraycopy(x, 0, Fx, 0, J); // TODO Remove Fx?

                tmp = previous; previous = x; x = tmp;
                mmUpdate(previous, x);
                ++iteration;

                copyDifference(secantsV[newestSecant], x, previous);

                // Do QN approximation here
//                auto M = secantsU.transpose() * (secantsU - secantsV);
//                auto Minv = M.inverse();
//                auto A = secantsU.transpose() * secantsU.col(newestSecant);
//                auto B = Minv * A;
//                auto C = secantsV * B;
//                VectorXd xqn = Fx + C;

                double M = 0;
                for (int j = 0; j < J; ++j) {
                    M += secantsU[0][j] * (secantsU[0][j] - secantsV[0][j]);
                }

                double Minv = 1.0 / M;

                double A = 0;
                for (int j = 0; j < J; ++j) {
                    A += secantsU[0][j] * secantsU[newestSecant][j];
                }

                double B = Minv * A;
                for (int j = 0; j < J; ++j) {
                    C[j] = secantsV[0][j] * B;
                }

                // New step
                for (int j = 0; j < J; ++j) {
                    x[j] = Fx[j] + C[j];
                }
                x[1] = 0.0; // Fixed point

                // Get ready for next secant-pair
                previousSecant = newestSecant;
                newestSecant = (newestSecant + 1) % qnQ;

                done = (convergenceCriterion(x, x0) < tolerance || iteration > maxIterations);
//                throw new RuntimeException("A");

            }



        }



        System.err.println("Finished in " + iteration + " iterations.");

        if (iteration >= maxIterations) {
            throw new NotConvergedException();
        }

        if (DEBUG) {
            System.err.println("Final  : " + printArray(x));
        }

       // throw new RuntimeException("out");

        return x;
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

    final int qnQ = 0;

    int iteration;

    private static final boolean DEBUG = true;
    private static final boolean PROGRESS = true;
}
