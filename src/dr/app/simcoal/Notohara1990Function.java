/*
 * Notohara1990Function.java
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

package dr.app.simcoal;

import dr.math.MultivariateFunction;
import dr.math.MultivariateMinimum;
import dr.math.ConjugateDirectionSearch;

public class Notohara1990Function implements MultivariateFunction {

    // the pairwise distances (in whatever units)
    double[] time;

    // the population of sequence 1
    double[] pop1;

    // the population of sequence 2
    double[] pop2;


    public Notohara1990Function(double[] time, double[] index1, double[] index2) {
        this.time = time;
        this.pop1 = index1;
        this.pop2 = index2;
    }

    public double[] optimize() {
        //MultivariateMinimum optimizer = new DifferentialEvolution(3,100);
        MultivariateMinimum optimizer = new ConjugateDirectionSearch();

        double[] params = new double[] { 100, 0.01, 0.01};

        optimizer.optimize(this, params, 1E-6, 1E-6);

        double score = evaluate(params);

        //System.out.println("score = " + score);

        return params;
    }

    public double evaluate(double[] argument) {
        return getSumOfSquares(argument[0], argument[0], argument[1], argument[2]);
    }

    private double getSumOfSquares(double N1, double N2, double m1, double m2) {

        double score = 0.0;

        int pop0Count, pop1Count;
        for (int i = 0; i < time.length; i++) {

            pop0Count = 0;
            pop1Count = 0;

            if (pop1[i] == 0.0) pop0Count += 1; else pop1Count += 1;
            if (pop2[i] == 0.0) pop0Count += 1; else pop1Count += 1;

            double expectedY = getExpectedT(pop0Count, pop1Count, N1, N2, m1, m2);
            double actualY = time[i];

            score += (expectedY - actualY) * (expectedY - actualY);
        }
        return score;
    }

    private double getExpectedT(int pop0Count, int pop1Count, double N1, double N2, double m1, double m2) {

        double m = m1 + m2;
        double topRight = 8.0 * N1 * N2 * m*m;
        double denom = m + 4.0*((N1*m1*m1) + (N2*m2*m2));

        double t;

        if (pop0Count == 2) {


            t = (2.0 * N1 * (3.0*m1 + m2)) + topRight;
            t /= denom;
            //System.out.println("E[T|(2,0)]=" + t);
            return t;

        } else if (pop0Count == 1 && pop1Count == 1) {

            t = (2.0 * N1 * (2.0*m1 + m2)) + (2.0 * N2 * (m1 + 2.0*m2)) + topRight + 1;
            t /= denom;
            //System.out.println("E[T|(1,1)]=" + t);
            return t;

        } else if (pop1Count == 2) {

            t = (2.0 * N2 * (m1 + 3.0*m2)) + topRight;
            t /= denom;
            //System.out.println("E[T|(0,2)]=" + t);
            return t;
        }
        throw new RuntimeException();

    }



    public int getNumArguments() {
        return 3;
    }

    public double getLowerBound(int n) {

        return 1e-8;
    }

    public double getUpperBound(int n) {

        return Double.MAX_VALUE;
    }
}
