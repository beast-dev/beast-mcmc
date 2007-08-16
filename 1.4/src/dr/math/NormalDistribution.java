/*
 * NormalDistribution.java
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
 * normal distribution (pdf, cdf, quantile)
 *
 * @version $Id: NormalDistribution.java,v 1.7 2005/05/24 20:26:01 rambaut Exp $
 *
 * @author Korbinian Strimmer
 */
public class NormalDistribution implements Distribution
{
	//
	// Public stuff
	//

	/**
	 * Constructor
	 */
	public NormalDistribution(double mean, double sd) {
		this.m = mean;
		this.sd = sd;
	}
	
	public double getMean() { return m; }
	public void setMean(double value) { m = value; }
	
	public double getSD() { return sd; }
	public void setSD(double value) { sd = value; }
	
	public double pdf(double x) { return pdf(x, m, sd); }
	public double logPdf(double x) { return logPdf(x, m, sd); }
	public double cdf(double x) { return cdf(x, m, sd); }
	public double quantile(double y) { return quantile(y, m, sd); }
	public double mean() { return mean(m, sd); }
	public double variance() { return variance(m, sd); }

	public final UnivariateFunction getProbabilityDensityFunction() { return pdfFunction; }
	
	private UnivariateFunction pdfFunction = new UnivariateFunction() {
		public final double evaluate(double x) { return pdf(x); }
		public final double getLowerBound() { return Double.NEGATIVE_INFINITY; }
		public final double getUpperBound() { return Double.POSITIVE_INFINITY; }
	};

	/**
	 * probability density function
	 *
	 * @param x argument
	 * @param m mean
	 * @param sd standard deviation
	 *
	 * @return pdf at x
	 */
	public static double pdf(double x, double m, double sd)
	{
		double a = 1.0/(Math.sqrt(2.0*Math.PI)*sd);
		double b = -(x-m)*(x-m)/(2.0*sd*sd);
		
		return a*Math.exp(b);
	}

	/**
	 * the natural log of the probability density function of the distribution 
	 *
	 * @param x argument
	 * @param m mean
	 * @param sd standard deviation
	 *
	 * @return log pdf at x
	 */
	public static double logPdf(double x, double m, double sd)
	{
		double a = 1.0/(Math.sqrt(2.0*Math.PI)*sd);
		double b = -(x-m)*(x-m)/(2.0*sd*sd);
		
		return Math.log(a) + b;
	}

	/**
	 * cumulative density function
	 *
	 * @param x argument
	 * @param m mean
	 * @param sd standard deviation
	 *
	 * @return cdf at x
	 */
	public static double cdf(double x, double m, double sd)
	{
		double a = (x-m)/(Math.sqrt(2.0)*sd);
		
		return 0.5*(1.0 + ErrorFunction.erf(a));
	}
	
	/**
	 * quantiles (=inverse cumulative density function)
	 *
	 * @param z argument
	 * @param m mean
	 * @param sd standard deviation
	 *
	 * @return icdf at z
	 */
	public static double quantile(double z, double m, double sd)
	{
		return m + Math.sqrt(2.0)*sd*ErrorFunction.inverseErf(2.0*z-1.0);
	}
	
	/**
	 * mean
	 *
	 * @param m mean
	 * @param sd standard deviation
	 *
	 * @return mean
	 */
	public static double mean(double m, double sd)
	{
		return m;
	}

	/**
	 * variance
	 *
	 * @param m mean
	 * @param sd standard deviation
	 *
	 * @return variance
	 */
	public static double variance(double m, double sd)
	{
		return sd*sd;
	}
	
	// Private
	
	protected double m, sd;
	
}
