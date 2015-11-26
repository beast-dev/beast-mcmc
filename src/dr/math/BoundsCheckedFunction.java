/*
 * BoundsCheckedFunction.java
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


/**
 * returns a very large number instead of the function value
 * if arguments are out of bound (useful for minimization with
 * minimizers that don't check argument boundaries)
 *
 * @author Korbinian Strimmer
 */
public class BoundsCheckedFunction implements MultivariateFunction {
    /**
     * construct bound-checked multivariate function
     * (a large number will be returned on function evaluation if argument
     * is out of bounds; default is 1000000)
     *
     * @param func unconstrained multivariate function
     */
    public BoundsCheckedFunction(MultivariateFunction func) {
        this(func, 1000000);
    }

    /**
     * construct constrained multivariate function
     *
     * @param func        unconstrained multivariate function
     * @param largeNumber value returned on function evaluation
     *                    if argument is out of bounds
     */
    public BoundsCheckedFunction(MultivariateFunction func, double largeNumber) {
        f = func;
        veryLarge = largeNumber;
    }

    /**
     * computes function value, taking into account the constraints on the
     * argument
     *
     * @param x function argument
     * @return function value (if argument is not in the predefined constrained area
     *         a very large number is returned instead of the true function value)
     */
    public double evaluate(double[] x) {
        int len = f.getNumArguments();

        for (int i = 0; i < len; i++) {
            if (x[i] < f.getLowerBound(i) ||
                    x[i] > f.getUpperBound(i)) {
                return veryLarge;
            }
        }

        return f.evaluate(x);
    }

    public int getNumArguments() {
        return f.getNumArguments();
    }

    public double getLowerBound(int n) {
        return f.getLowerBound(n);
    }

    public double getUpperBound(int n) {
        return f.getUpperBound(n);
    }

    //
    // Private stuff
    //

    private MultivariateFunction f;
    private double veryLarge;
}
