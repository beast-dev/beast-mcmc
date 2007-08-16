/*
 * DiscreteStatistics.java
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

package dr.stats;

import dr.util.HeapSort;

/**
 * simple discrete statistics (mean, variance, cumulative probability, quantiles etc.)
 *
 * @version $Id: DiscreteStatistics.java,v 1.11 2006/07/02 21:14:53 rambaut Exp $
 *
 * @author Korbinian Strimmer
 * @author Alexei Drummond
 */
public class DiscreteStatistics {

	//
	// Public stuff
	//

	/**
	 * compute mean
	 *
	 * @param x list of numbers
	 *
	 * @return mean
	 */
	public static double mean(double[] x)
	{
		double m = 0;
		int len = x.length;
		for (int i = 0; i < len; i++)
		{
			m += x[i];
		}

		return m/(double) len;
	}

    /**
     * compute median
     *
     * @param x list of numbers
     *
     * @return median
     */
    public static double median(double[] x) {

        int[] indices = new int[x.length];
        HeapSort.sort(x, indices);

        if (x == null || x.length == 0) {
	        throw new IllegalArgumentException();
        }

        int pos = x.length/2;
        if (x.length % 2 == 1) {
            return x[indices[pos]];
        } else {
            return (x[indices[pos-1]]+x[indices[pos]])/2.0;
        }
    }


	/**
	 * compute variance (ML estimator)
	 *
	 * @param x list of numbers
	 * @param mean assumed mean of x
	 *
	 * @return variance of x (ML estimator)
	 */
	public static double variance(double[] x, double mean)
	{
		double var = 0;
		int len = x.length;
		for (int i = 0; i < len; i++)
		{
			double diff = x[i]-mean;
			var += diff*diff;
		}

		int n;
		if (len < 2)
		{
			n = 1; // to avoid division by zero
		}
		else
		{
			n = len-1; // for ML estimate
		}

		return var/ (double) n;
	}


	/**
	 * compute covariance
	 *
	 * @param x list of numbers
	 * @param y list of numbers
	 *
	 * @return covariance of x and y
	 */
	public static double covariance(double[] x, double[] y) {

		return covariance(x, y, mean(x), mean(y), stdev(x), stdev(y));
	}

	/**
	 * compute covariance
	 *
	 * @param x list of numbers
	 * @param y list of numbers
	 * @param xmean assumed mean of x
	 * @param ymean assumed mean of y
	 * @param xstdev assumed stdev of x
	 * @param ystdev assumed stdev of y
	 *
	 * @return covariance of x and y
	 */
	public static double covariance(double[] x, double[] y, double xmean, double ymean, double xstdev, double ystdev) {

		if (x.length != y.length) throw new IllegalArgumentException("x and y arrays must be same length!");

		double covar = 0.0;
		for (int i =0; i < x.length; i++) {
			covar += (x[i]-xmean)*(y[i]-ymean);
		}
		covar /= x.length;
		covar /= (xstdev*ystdev);
		return covar;
	}

	/**
	 * compute fisher skewness
	 *
	 * @param x list of numbers
	 *
	 * @return skewness of x
	 */
	public static double skewness(double[] x) {

		double mean = mean(x);
		double stdev = stdev(x);
		double skew = 0.0;
		double len = x.length;

		for (int i = 0; i < x.length; i++)
		{
			double diff = x[i]-mean;
			diff /= stdev;

			skew += (diff*diff*diff);
		}

		skew *= (len / ((len - 1) * (len - 2)));

		return skew;
	}

	/**
	 * compute standard deviation
	 *
	 * @param x list of numbers
	 *
	 * @return standard deviation of x
	 */
	public static double stdev(double[] x) {
		return Math.sqrt(variance(x));
	}

	/**
	 * compute variance (ML estimator)
	 *
	 * @param x list of numbers
	 *
	 * @return variance of x (ML estimator)
	 */
	public static double variance(double[] x)
	{
		double m = mean(x);
		return variance(x, m);
	}


	/**
	 * compute variance of sample mean (ML estimator)
	 *
	 * @param x list of numbers
	 * @param mean assumed mean of x
	 *
	 * @return variance of x (ML estimator)
	 */
	public static double varianceSampleMean(double[] x, double mean)
	{
		return variance(x, mean)/(double) x.length;
	}

	/**
	 * compute variance of sample mean (ML estimator)
	 *
	 * @param x list of numbers
	 *
	 * @return variance of x (ML estimator)
	 */
	public static double varianceSampleMean(double[] x)
	{
		return variance(x)/(double) x.length;
	}


	/**
	 * compute the q-th quantile for a distribution of x
	 * (= inverse cdf)
	 *
	 * @param q quantile (0 < q <= 1)
	 * @param x discrete distribution (an unordered list of numbers)
	 * @param indices index sorting x
	 *
	 * @return q-th quantile
	 */
	public static double quantile(double q, double[] x, int[] indices)
	{
		if (q < 0.0 || q > 1.0) throw new IllegalArgumentException("Quantile out of range");

		if (q == 0.0)
		{
			// for q==0 we have to "invent" an entry smaller than the smallest x

			return x[indices[0]] - 1.0;
		}

		return x[indices[(int) Math.ceil(q*indices.length)-1]];
	}

	/**
	 * compute the q-th quantile for a distribution of x
	 * (= inverse cdf)
	 *
	 * @param q quantile (0 <= q <= 1)
	 * @param x discrete distribution (an unordered list of numbers)
	 *
	 * @return q-th quantile
	 */
	public static double quantile(double q, double[] x)
	{
		int[] indices = new int[x.length];
		HeapSort.sort(x, indices);

		return quantile(q, x, indices);
	}

    /**
     * compute the q-th quantile for a distribution of x
     * (= inverse cdf)
     *
     * @param q quantile (0 <= q <= 1)
     * @param x discrete distribution (an unordered list of numbers)
     *
     * @return q-th quantile
     */
    public static double quantile(double q, double[] x, int count)
    {
        int[] indices = new int[count];
        HeapSort.sort(x, indices);

        return quantile(q, x, indices);
    }

	/**
	 * compute the cumulative probability Pr(x <= z) for a given z
	 * and a distribution of x
	 *
	 * @param z threshold value
	 * @param x discrete distribution (an unordered list of numbers)
	 * @param indices index sorting x
	 *
	 * @return cumulative probability
	 */
	public static double cdf(double z, double[] x, int[] indices)
	{
		int i;
		for (i = 0; i < x.length; i++)
		{
			if (x[indices[i]] > z) break;
		}

		return (double) i/ (double) x.length;
	}

	/**
	 * compute the cumulative probability Pr(x <= z) for a given z
	 * and a distribution of x
	 *
	 * @param z threshold value
	 * @param x discrete distribution (an unordered list of numbers)
	 *
	 * @return cumulative probability
	 */
	public static double cdf(double z, double[] x)
	{
		int[] indices = new int[x.length];
		HeapSort.sort(x, indices);

		return cdf(z, x, indices);
	}

    public static final double max(double[] x) {
        double max = x[0];
        for (int i = 1; i < x.length; i++) {
            if (x[i] > max) max = x[i];
        }
        return max;
    }

    public static final double min(double[] x) {
        double min = x[0];
        for (int i = 1; i < x.length; i++) {
            if (x[i] < min) min = x[i];
        }
        return min;
    }
}
