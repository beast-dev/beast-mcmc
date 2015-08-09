/*
 * MultivariateGammaDistribution.java
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

package dr.math.distributions;

/**
 * @author Marc Suchard
 * @author Guy Baele
 */
public class MultivariateGammaDistribution implements MultivariateDistribution {

    //TODO: Currently this implements a product of independent Gammas, need to re-code as true multivariate distribution

    public static final String TYPE = "multivariateGamma";

    public MultivariateGammaDistribution(double[] shape, double[] scale) {

        if (shape.length != scale.length)
            throw new RuntimeException("Creation error in MultivariateGammaDistribution");

        dim = shape.length;

        this.shape = shape;
        this.scale = scale;
        
        this.flags = new boolean[dim];
        for (int i = 0; i < dim; i++) {
        	flags[i] = true;
        }

    }
    
    public MultivariateGammaDistribution(double[] shape, double[] scale, boolean[] flags) {

        if (shape.length != scale.length)
            throw new RuntimeException("Creation error in MultivariateGammaDistribution");

        dim = shape.length;

        this.shape = shape;
        this.scale = scale;
        
        this.flags = flags;

    }

    public double logPdf(double[] x) {

        double logPdf = 0;

        if (x.length != dim) {
            throw new IllegalArgumentException("data array is of the wrong dimension");
        }

        for (int i = 0; i < dim; i++) {
        	if (flags[i]) { 
        		logPdf += GammaDistribution.logPdf(x[i], shape[i], scale[i]);
        	}
        }

        return logPdf;
    }

    public double[][] getScaleMatrix() {
        throw new RuntimeException("Not yet implemented");
    }

    public double[] getMean() {
        throw new RuntimeException("Not yet implemented");
    }

    public String getType() {
        return TYPE;
    }

    private double[] shape;
    private double[] scale;
    private int dim;
    
    //for each flag that is true, add the logPdf of that gamma distribution to the overall logPdf
    private boolean[] flags;
    
}
