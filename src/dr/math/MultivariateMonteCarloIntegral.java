/*
 * MonteCarloIntegral.java
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

package dr.math;

/**
 * Approximates the integral of a given function using Monte Carlo integration
 *
 * @author Alexei Drummond
 *
 * @version $Id: MonteCarloIntegral.java,v 1.5 2005/05/24 20:26:01 rambaut Exp $
 */
public class MultivariateMonteCarloIntegral implements MultivariateIntegral {

    public MultivariateMonteCarloIntegral(int sampleSize) {
        this.sampleSize = sampleSize;
    }

    /**
     * @return the approximate integral of the given function
     * within the given range using simple monte carlo integration.
     * @param f the function whose integral is of interest
     * @param mins the minimum value of the function
     * @param maxes the  upper limit of the function
     */
    public double integrate(MultivariateFunction f, double[] mins, double[] maxes) {

        double integral = 0.0;

        double area = 1;

        for(int i=0; i<f.getNumArguments(); i++){
            area *= maxes[i]-mins[i];
        }

        for (int i=1; i <= sampleSize; i++) {

            double[] sample = new double[f.getNumArguments()];
            for(int j=0; j<sample.length; j++){
                sample[j] = MathUtils.nextDouble()*(maxes[j]-mins[j]);
            }

            integral += f.evaluate(sample);
        }
        integral *= area/(double)sampleSize;
        return integral;
    }

    private int sampleSize;
}
