/*
 * MathUtils.java
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

import java.text.NumberFormat;
import java.text.ParseException;

import dr.util.NumberFormatter;

/**
 * Handy utility functions which have some Mathematical relavance.
 *
 * @author Matthew Goode
 * @author Alexei Drummond
 * @author Gerton Lunter
 * @version $Id: MathUtils.java,v 1.13 2006/08/31 14:57:24 rambaut Exp $
 */
public class MathUtils {

	private MathUtils() {
	}

	/**
	 * A random number generator that is initialized with the clock when this
	 * class is loaded into the JVM. Use this for all random numbers.
	 * Note: This method or getting random numbers in not thread-safe. Since
	 * MersenneTwisterFast is currently (as of 9/01) not synchronized using
	 * this function may cause concurrency issues. Use the static get methods of the
	 * MersenneTwisterFast class for access to a single instance of the class, that
	 * has synchronization.
	 */
	private static final MersenneTwisterFast random = MersenneTwisterFast.DEFAULT_INSTANCE;

	// Chooses one category if a cumulative probability distribution is given
	public static int randomChoice(double[] cf) {

		double U = MathUtils.nextDouble();

		int s;
		if (U <= cf[0]) {
			s = 0;
		} else {
			for (s = 1; s < cf.length; s++) {
				if (U <= cf[s] && U > cf[s - 1]) {
					break;
				}
			}
		}

		return s;
	}


	/**
	 * @param pdf array of unnormalized probabilities
	 * @return a sample according to an unnormalized probability distribution
	 */
	public static int randomChoicePDF(double[] pdf) {

		double U = MathUtils.nextDouble() * getTotal(pdf);
		for (int i = 0; i < pdf.length; i++) {

			U -= pdf[i];
			if (U < 0.0) {
				return i;
			}

		}
		for (int i = 0; i < pdf.length; i++) {
			System.out.println(i + "\t" + pdf[i]);
		}
		throw new Error("randomChoicePDF falls through -- negative, infinite or NaN components in input " +
				"distribution, or all zeroes?");
	}

	/**
	 * @param logpdf array of unnormalised log probabilities
	 * @return a sample according to an unnormalised probability distribution
	 * <p/>
	 * Use this if probabilities are rounding to zero when converted to real space
	 */
	public static int randomChoiceLogPDF(double[] logpdf) {

		double scalingFactor = Double.NEGATIVE_INFINITY;

		for (double aLogpdf : logpdf) {
			if (aLogpdf > scalingFactor) {
				scalingFactor = aLogpdf;
			}
		}

		if (scalingFactor == Double.NEGATIVE_INFINITY) {
			throw new Error("randomChoiceLogPDF falls through -- all -INF components in input distribution");
		}

		for (int j = 0; j < logpdf.length; j++) {
			logpdf[j] = logpdf[j] - scalingFactor;
		}

		double[] pdf = new double[logpdf.length];

		for (int j = 0; j < logpdf.length; j++) {
			pdf[j] = Math.exp(logpdf[j]);
		}

		return randomChoicePDF(pdf);

	}

	/**
	 * @param array to normalize
	 * @return a new double array where all the values sum to 1.
	 * Relative ratios are preserved.
	 */
	public static double[] getNormalized(double[] array) {
		double[] newArray = new double[array.length];
		double total = getTotal(array);
		for (int i = 0; i < array.length; i++) {
			newArray[i] = array[i] / total;
		}
		return newArray;
	}


	/**
	 * @param array entries to be summed
	 * @param start start position
	 * @param end   the index of the element after the last one to be included
	 * @return the total of a the values in a range of an array
	 */
	public static double getTotal(double[] array, int start, int end) {
		double total = 0.0;
		for (int i = start; i < end; i++) {
			total += array[i];
		}
		return total;
	}

	/**
	 * @param array to sum over
	 * @return the total of the values in an array
	 */
	public static double getTotal(double[] array) {
		return getTotal(array, 0, array.length);

	}

	// ===================== (Synchronized) Static access methods to the private random instance ===========

	/**
	 * Access a default instance of this class, access is synchronized
	 */
	public static long getSeed() {
		synchronized (random) {
			return random.getSeed();
		}
	}

	/**
	 * Access a default instance of this class, access is synchronized
	 */
	public static void setSeed(long seed) {
		synchronized (random) {
			random.setSeed(seed);
		}
	}

	/**
	 * Access a default instance of this class, access is synchronized
	 */
	public static byte nextByte() {
		synchronized (random) {
			return random.nextByte();
		}
	}

	/**
	 * Access a default instance of this class, access is synchronized
	 */
	public static boolean nextBoolean() {
		synchronized (random) {
			return random.nextBoolean();
		}
	}

	/**
	 * Access a default instance of this class, access is synchronized
	 */
	public static void nextBytes(byte[] bs) {
		synchronized (random) {
			random.nextBytes(bs);
		}
	}

	/**
	 * Access a default instance of this class, access is synchronized
	 */
	public static char nextChar() {
		synchronized (random) {
			return random.nextChar();
		}
	}

	/**
	 * Access a default instance of this class, access is synchronized
	 */
	public static double nextGaussian() {
		synchronized (random) {
			return random.nextGaussian();
		}
	}

	//Mean = alpha / lambda
	//Variance = alpha / (lambda*lambda)

	public static double nextGamma(double alpha, double lambda) {
		synchronized (random) {
			return random.nextGamma(alpha, lambda);
		}
	}

	//Mean = alpha/(alpha+beta)
	//Variance = (alpha*beta)/(alpha+beta)^2*(alpha+beta+1)

	public static double nextBeta(double alpha, double beta) {
		double x = nextGamma(alpha, 1);
		double y = nextGamma(beta, 1);
		return x / (x + y);
	}


	/**
	 * Access a default instance of this class, access is synchronized
	 *
	 * @return a pseudo random double precision floating point number in [01)
	 */
	public static double nextDouble() {
		synchronized (random) {
			return random.nextDouble();
		}
	}

	/**
	 * @return log of random variable in [0,1]
	 */
	public static double randomLogDouble() {
		return Math.log(nextDouble());
	}

	/**
	 * Access a default instance of this class, access is synchronized
	 */
	public static double nextExponential(double lambda) {
		synchronized (random) {
			return -1.0 * Math.log(1 - random.nextDouble()) / lambda;
		}
	}

	/**
	 * Access a default instance of this class, access is synchronized
	 */
	public static double nextInverseGaussian(double mu, double lambda) {
		synchronized (random) {
			/* CODE TAKEN FROM WIKIPEDIA. TESTING DONE WITH RESULTS GENERATED IN R AND LOOK COMPARABLE */
			double v = random.nextGaussian();   // sample from a normal distribution with a mean of 0 and 1 standard deviation
			double y = v * v;
			double x = mu + (mu * mu * y) / (2 * lambda) - (mu / (2 * lambda)) * Math.sqrt(4 * mu * lambda * y + mu * mu * y * y);
			double test = MathUtils.nextDouble();  // sample from a uniform distribution between 0 and 1
			if (test <= (mu) / (mu + x)) {
				return x;
			} else {
				return (mu * mu) / x;
			}
		}
	}


	/**
	 * Access a default instance of this class, access is synchronized
	 */
	public static float nextFloat() {
		synchronized (random) {
			return random.nextFloat();
		}
	}

	/**
	 * Access a default instance of this class, access is synchronized
	 */
	public static long nextLong() {
		synchronized (random) {
			return random.nextLong();
		}
	}

	/**
	 * Access a default instance of this class, access is synchronized
	 */
	public static short nextShort() {
		synchronized (random) {
			return random.nextShort();
		}
	}

	/**
	 * Access a default instance of this class, access is synchronized
	 */
	public static int nextInt() {
		synchronized (random) {
			return random.nextInt();
		}
	}

	/**
	 * Access a default instance of this class, access is synchronized
	 */
	public static int nextInt(int n) {
		synchronized (random) {
			return random.nextInt(n);
		}
	}

	/**
	 * @param low
	 * @param high
	 * @return uniform between low and high
	 */
	public static double uniform(double low, double high) {
		return low + nextDouble() * (high - low);
	}

	/**
	 * Shuffles an array.
	 */
	public static void shuffle(int[] array) {
		synchronized (random) {
			random.shuffle(array);
		}
	}

	/**
	 * Shuffles an array. Shuffles numberOfShuffles times
	 */
	public static void shuffle(int[] array, int numberOfShuffles) {
		synchronized (random) {
			random.shuffle(array, numberOfShuffles);
		}
	}

	/**
	 * Returns an array of shuffled indices of length l.
	 *
	 * @param l length of the array required.
	 */
	public static int[] shuffled(int l) {
		synchronized (random) {
			return random.shuffled(l);
		}
	}


	public static int[] sampleIndicesWithReplacement(int length) {
		synchronized (random) {
			int[] result = new int[length];
			for (int i = 0; i < length; i++)
				result[i] = random.nextInt(length);
			return result;
		}
	}

	/**
	 * Permutes an array.
	 */
	public static void permute(int[] array) {
		synchronized (random) {
			random.permute(array);
		}
	}

	/**
	 * Returns a uniform random permutation of 0,...,l-1
	 *
	 * @param l length of the array required.
	 */
	public static int[] permuted(int l) {
		synchronized (random) {
			return random.permuted(l);
		}
	}


	public static double logHyperSphereVolume(int dimension, double radius) {
		return dimension * (0.5723649429247001 + Math.log(radius)) +
				-GammaFunction.lnGamma(dimension / 2.0 + 1.0);
	}

	/**
	 * Returns sqrt(a^2 + b^2) without under/overflow.
	 */
	public static double hypot(double a, double b) {
		double r;
		if (Math.abs(a) > Math.abs(b)) {
			r = b / a;
			r = Math.abs(a) * Math.sqrt(1 + r * r);
		} else if (b != 0) {
			r = a / b;
			r = Math.abs(b) * Math.sqrt(1 + r * r);
		} else {
			r = 0.0;
		}
		return r;
	}

	/**
	 * return double *.????
	 *
	 * @param value
	 * @param sf
	 * @return
	 */
	public static double round(double value, int sf) {
		NumberFormatter formatter = new NumberFormatter(sf);
		try {
			return NumberFormat.getInstance().parse(formatter.format(value)).doubleValue();
		} catch (ParseException e) {
			return value;
		}
	}

	public static int[] getRandomState() {
		synchronized (random) {
			return random.getRandomState();
		}
	}

	public static void setRandomState(int[] rngState) {
		synchronized (random) {
			random.setRandomState(rngState);
		}
	}
}