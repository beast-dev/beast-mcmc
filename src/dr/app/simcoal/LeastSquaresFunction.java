/*
 * LeastSquaresFunction.java
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

public class LeastSquaresFunction implements MultivariateFunction {

    double[] x, y;

    public LeastSquaresFunction(double[] x, double[] y) {
        this.x = x;
        this.y = y;
    }

    public double[] optimize() {
        //MultivariateMinimum optimizer = new DifferentialEvolution(2,10);
        MultivariateMinimum optimizer = new ConjugateDirectionSearch();

        double[] params = new double[] { 0.1, 0.1};

        optimizer.optimize(this, params, 1E-6, 1E-6);

        return params;
    }

    public double evaluate(double[] argument) {
        return getSumOfSquares(argument[0], argument[1]);
    }

    private double getSumOfSquares(double m, double c) {

        double score = 0.0;

        for (int i = 0; i < x.length; i++) {
            double expectedY = m * x[i] + c;
            double actualY = y[i];

            score += (expectedY - actualY) * (expectedY - actualY);
        }
        return score;
    }

    public int getNumArguments() {
        return 2;
    }

    public double getLowerBound(int n) {

        return -Double.MAX_VALUE;
    }

    public double getUpperBound(int n) {

        return Double.MAX_VALUE;
    }
}
