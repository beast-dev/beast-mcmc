/*
 * MultivariateKDEDistribution.java
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
 * @author Guy Baele
 */
public class MultivariateKDEDistribution implements MultivariateDistribution {
	
	public static final String TYPE = "multivariateKDE";
    public static final boolean DEBUG = false;
	
	private Distribution[] multivariateKDE;
	private int dimension;
	//private boolean[] flags;

	public MultivariateKDEDistribution (Distribution[] multivariateKDE) {
		
		if (multivariateKDE.length <= 0) {
			throw new RuntimeException("Creation error in MultivariateKDEDistribution(Distribution[] multivariateKDE)");
		}
		
		this.multivariateKDE = multivariateKDE;
		this.dimension = multivariateKDE.length;
		/*for (int i = 0; i < dimension; i++) {
			flags[i] = true;
		}*/

	}
	
	public MultivariateKDEDistribution (Distribution[] multivariateKDE, boolean[] flags) {
		
		if (multivariateKDE.length <= 0) {
			throw new RuntimeException("Creation error in MultivariateKDEDistribution(Distribution[] multivariateKDE, boolean[] flags)");
		}
		
		this.multivariateKDE = multivariateKDE;
		this.dimension = multivariateKDE.length;
		//this.flags = flags;

	}

	public double logPdf(double[] x) {
		
		double logPdf = 0;
		
		if (x.length != dimension) {
            throw new IllegalArgumentException("data array is of the wrong dimension");
        }
		
		for (int i = 0; i < dimension; i++) {
			//if (flags[i]) {
			logPdf += multivariateKDE[i].logPdf(x[i]);
			//}
		}

        if (DEBUG){
            System.err.println("MultivariateKDEDistribution, dimension = " + dimension);
            for (int i = 0; i < dimension; i++) {
                System.err.println(i + ", " + "x[i] = " + x[i] + ", logPdf = " + multivariateKDE[i].logPdf(x[i]));
                //System.err.println("    mean = " + multivariateKDE[i].mean() + ", variance = " + multivariateKDE[i].variance());
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
	
}
