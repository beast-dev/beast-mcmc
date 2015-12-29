/*
 * MultivariateMonteCarloIntegral.java
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

import java.util.Arrays;
import java.util.HashMap;

/**
 * Approximates the integral of a given function using Monte Carlo integration; possibly stratified
 *
 * @author Alexei Drummond
 *
 * @version $Id: MonteCarloIntegral.java,v 1.5 2005/05/24 20:26:01 rambaut Exp $
 */
public class MultivariateMonteCarloIntegral implements MultivariateIntegral {

    // bins is the number of divisions that each axis it split into (so the full number of bins is this to the
    // power of the dimension of the function
    // sampleSize is the number of samples PER BIN

    public MultivariateMonteCarloIntegral(int sampleSize, int bins) {
        this.sampleSize = sampleSize;
        this.bins = bins;
    }

    public MultivariateMonteCarloIntegral(int sampleSize) {
        this(sampleSize, 1);
    }

    /**
     * @return the approximate integral of the given function
     * within the given range using simple monte carlo integration.
     * @param f the function whose integral is of interest
     * @param mins the minimum value of the function
     * @param maxes the  upper limit of the function
     */
    public double integrate(MultivariateFunction f, double[] mins, double[] maxes) {

        int dim = f.getNumArguments();
        int totalBins = bins*dim;
        double[] steps = new double[dim];
        double totalArea=1;


        for(int i=0; i<dim; i++){
            totalArea *= (maxes[i]-mins[i]);
        }

        HashMap<Integer, double[]> binCorners = new HashMap<Integer, double[]>();
        double[] currentCorner = new double[dim];


        for(int index=0; index<totalBins; index++){
            binCorners.put(index, Arrays.copyOf(currentCorner, dim));

            int dimToCheck = 0;
            while(dimToCheck<dim){
                if(currentCorner[dimToCheck]+steps[dimToCheck]<maxes[dimToCheck]){
                    currentCorner[dimToCheck] += steps[dimToCheck];
                    break;
                } else {
                    currentCorner[dimToCheck] = mins[dimToCheck];
                }
                dimToCheck++;
            }
        }

        double integral = 0.0;

        for(int i=0; i<totalBins; i++){

            for (int j=1; j <= sampleSize; j++) {

                double[] sample = new double[dim];
                for(int k=0; k<sample.length; k++){
                    sample[k] = binCorners.get(i)[k] + MathUtils.nextDouble()*(steps[k]);
                }

                integral += f.evaluate(sample);
            }

        }
        integral *= totalArea/((double)sampleSize*totalBins);
        return integral;
    }

    protected int getSampleSize(){
        return sampleSize;
    }
    protected int getBins(){
        return bins;
    }

    private int sampleSize;
    private int bins;
}
