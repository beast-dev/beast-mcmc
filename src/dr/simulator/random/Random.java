/*
 * Random.java
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

package dr.simulator.random;

/**
 * Handy utility functions which have some Mathematical relavance.
 *
 * @author Matthew Goode
 * @author Alexei Drummond
 *
 * @version $Id: Random.java,v 1.2 2005/04/28 16:51:44 rambaut Exp $
 */
public class Random {

	private Random() {}

	/**
	 * A random number generator that is initialized with the clock when this
	 * class is loaded into the JVM. Use this for all random numbers.
	 * Note: This method or getting random numbers in not thread-safe. Since
	 * MersenneTwisterFast is currently (as of 9/01) not synchronized using
	 * this function may cause concurrency issues. Use the static get methods of the
	 * MersenneTwisterFast class for access to a single instance of the class, that
	 * has synchronization.
	 */
	private static final MersenneTwisterFast random = new MersenneTwisterFast();

    // Chooses one category if a cumulative probability distribution is given
	public static int randomChoice(double[] cf)
	{

		double U = random.nextDouble();

		int s;
		if (U <= cf[0])
		{
			s = 0;
		}
		else
		{
			for (s = 1; s < cf.length; s++)
			{
				if (U <= cf[s] && U > cf[s-1])
				{
					break;
				}
			}
		}

		return s;
	}

    /**
     * Shuffles an array.
     */
    public static final void shuffle(int[] array) {
        int l = array.length;
        for (int i = 0; i < l; i++) {
            int index = nextInt(l-i) + i;
            int temp = array[index];
            array[index] = array[i];
            array[i] = temp;
        }
    }
    /**
     * Shuffles an array. Shuffles numberOfShuffles times
     */
    public static final void shuffle(int[] array, int numberOfShuffles) {
        int i, j, temp, l = array.length;
        for (int shuffle = 0; shuffle < numberOfShuffles; shuffle++) {
            do {
                i = nextInt(l);
                j = nextInt(l);
            } while(i!=j);
            temp = array[j];
            array[j] = array[i];
            array[i] = temp;
        }
    }
    /**
     * Returns an array of shuffled indices of length l.
     * @param l length of the array required.
     */
    public static int[] shuffled(int l) {

        int[] array = new int[l];

        // initialize array
        for (int i = 0; i < l; i++) {
            array[i] = i;
        }
        shuffle(array);

        return array;
    }

    /**
     * Returns an integer value that is a
     * random deviate drawn from a Poisson distribution of mean xm.
     */
    public static int nextPoisson(double xm) {

        double em; // expected mean
        double t, y;

        if (xm < 12.0) { // use direct method
            if (xm != oldm) {
                oldm = xm;
                g = Math.exp(-xm);  // if xm is new compute the exponential
            }
            em = -1.0;
            t = 1.0;
            do {
                ++em;
                t *= nextDouble();
            } while (t > g);
        } else {
            if (xm != oldm) {
                oldm = xm;
                sq = Math.sqrt(2.0 * xm);
                alxm = Math.log(xm);
                g = xm * alxm - gammln(xm + 1.0);
                // The function gammln is the natural log of the gamma function
            }
            do {
                do {
                    y = Math.tan(Math.PI * nextDouble());
                    em = sq * y + xm; // em is y shifted and scaled
                } while (em < 0.0); // reject if in realm of zero probability
                em = Math.floor(em);
                t = 0.9 *(1.0 + y*y) * Math.exp(em*alxm-gammln(em + 1.0)-g);
                // The ratio of the desired distribution to the comparison function;
                // we accept or reject by comparing it to another uniform deviate.
                // The factor 0.9 is chosen so that t never exceeds 1.
            } while (nextDouble() > t);
        }
        return (int)em;
    }

    // oldm is a flag for whether xm has been changed since last call.
    private static double oldm = -1.0;
    private static double sq, alxm, g;

    // used by gammaln
    private static final double[] cof = {   76.18009172947146, -86.50532032941677,
                                            24.01409824083091,  -1.231739572450155,
                                            0.1208650973866179e-2, -0.5395239384953e-5};


    /**
     * Returns the value of log gamma(xx) for xx > 0.
     */
    private static double gammln(double xx) {
        double x, y, tmp, ser;

        int j;

        y = x = xx;
        tmp = x + 5.5;
        tmp -= (x + 0.5) * Math.log(tmp);
        ser = 1.000000000190015;
        for (j = 0; j <= 5; j++) {
            ser += cof[j] / ++y;
        }
        return -tmp + Math.log(2.5066282746310005 * ser / x);
    }

    // ===================== Static access methods to the private random instance ===========

    /** Access a default instance of this class, access is synchronized */
    public static final void setSeed(long seed) {
        random.setSeed(seed);
    }
    /** Access a default instance of this class, access is synchronized */
    public static final byte nextByte() {
        return random.nextByte();
    }
    /** Access a default instance of this class, access is synchronized */
    public static final boolean nextBoolean() {
        return random.nextBoolean();
    }
    public static final boolean nextBoolean(final float probability) {
        return random.nextBoolean(probability);
    }
    public static final boolean nextBoolean(final double probability) {
        return random.nextBoolean(probability);
    }
    /** Access a default instance of this class, access is synchronized */
    public static final void nextBytes(byte[] bs) {
        random.nextBytes(bs);
    }
    /** Access a default instance of this class, access is synchronized */
    public static final char nextChar() {
        return random.nextChar();
    }
    /** Access a default instance of this class, access is synchronized */
    public static final double nextGaussian() {
        return random.nextGaussian();
    }
    /** Access a default instance of this class, access is synchronized */
    public static final double nextDouble() {
        return random.nextDouble();
    }
    /** Access a default instance of this class, access is synchronized */
    public static final float nextFloat() {
        return random.nextFloat();
    }
    /** Access a default instance of this class, access is synchronized */
    public static final long nextLong() {
        return random.nextLong();
    }
    /** Access a default instance of this class, access is synchronized */
    public static final long nextLong(final long n) {
        return random.nextLong(n);
    }
    /** Access a default instance of this class, access is synchronized */
    public static final short nextShort() {
        return random.nextShort();
    }
    /** Access a default instance of this class, access is synchronized */
    public static final int nextInt() {
        return random.nextInt();
    }
    /** Access a default instance of this class, access is synchronized */
    public static final int nextInt(int n) {
        return random.nextInt(n);
    }


}
