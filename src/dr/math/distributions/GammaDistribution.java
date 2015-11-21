/*
 * GammaDistribution.java
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

import cern.jet.random.Gamma;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;
import dr.math.GammaFunction;
import dr.math.MathUtils;
import dr.math.UnivariateFunction;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.GammaDistributionImpl;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * gamma distribution.
 * <p/>
 * (Parameters: shape, scale; mean: scale*shape; variance: scale^2*shape)
 *
 * @author Korbinian Strimmer
 * @author Gerton Lunter
 * @version $Id: GammaDistribution.java,v 1.9 2006/03/30 11:12:47 rambaut Exp $
 */
public class GammaDistribution implements Distribution {
    //
    // Public stuff
    //

    /**
     * Constructor
     */
    public GammaDistribution(double shape, double scale) {
        this.shape = shape;
        this.scale = scale;
        this.samples = 0;

        if (TRY_COLT) {
            randomEngine = new MersenneTwister(MathUtils.nextInt());
            System.out.println("Colt Gamma(" + shape + "," + scale + ")");
            coltGamma = new Gamma(shape, 1.0/scale, randomEngine);
        }
    }

    public double getShape() {
        return shape;
    }

    public void setShape(double value) {
        shape = value;
    }

    public double getScale() {
        return scale;
    }

    public void setScale(double value) {
        scale = value;
    }

    public double pdf(double x) {
        if (TRY_COLT) {
            return coltGamma.pdf(x);
        } else {
            return pdf(x, shape, scale);
        }
    }

    public double logPdf(double x) {
        if (TRY_COLT) {
            return Math.log(coltGamma.pdf(x));
        } else {
            return logPdf(x, shape, scale);
        }
    }

    public double cdf(double x) {
        if (TRY_COLT) {
            return coltGamma.cdf(x);
        } else {
            return cdf(x, shape, scale);
        }
    }

    public double quantile(double y) {
        return quantile(y, shape, scale);
    }

    public double mean() {
        return mean(shape, scale);
    }

    public double variance() {
        return variance(shape, scale);
    }

    public double nextGamma() {
        return nextGamma(shape, scale);
    }

    public final UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private final UnivariateFunction pdfFunction = new UnivariateFunction() {
        public final double evaluate(double x) {
            return pdf(x);
        }

        public final double getLowerBound() {
            return 0.0;
        }

        public final double getUpperBound() {
            return Double.POSITIVE_INFINITY;
        }
    };

    /**
     * probability density function of the Gamma distribution
     *
     * @param x     argument
     * @param shape shape parameter
     * @param scale scale parameter
     * @return pdf value
     */
    public static double pdf(double x, double shape, double scale) {
        // return Math.pow(scale,-shape)*Math.pow(x, shape-1.0)/
        // Math.exp(x/scale + GammaFunction.lnGamma(shape));
        if (x < 0)  return 0; // to make BEAUti plot continue
//            throw new IllegalArgumentException();
        if (x == 0) {
            if (shape == 1.0)
                return 1.0 / scale;
            else
                return 0.0;
        }

        if (shape == 0.0)  // uninformative
            return 1.0 / x;

        if (shape == -0.5) { // Gelman 2008, hierarchical variance, -1 degrees of freedom
            return Math.sqrt(x);
        }

        final double xs = x / scale;

        if (shape == 1.0) {
            return Math.exp(-xs) / scale;
        }

        final double a = Math.exp((shape - 1.0) * Math.log(xs) - xs
                - GammaFunction.lnGamma(shape));

        return a / scale;
    }

    /**
     * the natural log of the probability density function of the distribution
     *
     * @param x     argument
     * @param shape shape parameter
     * @param scale scale parameter
     * @return log pdf value
     */
    public static double logPdf(double x, double shape, double scale) {
        // double a = Math.pow(scale,-shape) * Math.pow(x, shape-1.0);
        // double b = x/scale + GammaFunction.lnGamma(shape);
        // return Math.log(a) - b;

        // AR - changed this to return -ve inf instead of throwing an
        // exception... This makes things
        // much easier when using this to calculate log likelihoods.
        // if (x < 0) throw new IllegalArgumentException();
        if (x < 0)
            return Double.NEGATIVE_INFINITY;

        if (x == 0) {
            if (shape == 1.0)
                return Math.log(1.0 / scale);
            else
                return Double.NEGATIVE_INFINITY;
        }
        if (shape == 1.0) {
            return (-x / scale) - Math.log(scale);
        }
        if (shape == 0.0)  // uninformative
            return -Math.log(x);

        if (shape == -0.5) { // Gelman 2008, hierarchical variance, -1 degrees of freedom
            return 0.5 * Math.log(x);
        }
        
        /*return ((shape - 1.0) * Math.log(x/scale) - x / scale - GammaFunction
                .lnGamma(shape))
                - Math.log(scale);*/
        
        return ((shape - 1.0) * (Math.log(x) - Math.log(scale)) - x / scale - GammaFunction
                .lnGamma(shape))
                - Math.log(scale);
    }

    /**
     * cumulative density function of the Gamma distribution
     *
     * @param x     argument
     * @param shape shape parameter
     * @param scale scale parameter
     * @return cdf value
     */
    public static double cdf(double x, double shape, double scale) {
        if (x < 0.0 || shape <= 0.0) {
            return 0;
        }
        return GammaFunction.incompleteGammaP(shape, x / scale);
    }

    /**
     * quantile (inverse cumulative density function) of the Gamma distribution
     *
     * @param y     argument
     * @param shape shape parameter
     * @param scale scale parameter
     * @return icdf value
     */
    public static double quantile(double y, double shape, double scale) {
        return 0.5 * scale * pointChi2(y, 2.0 * shape);
    }

    /**
     * mean of the Gamma distribution
     *
     * @param shape shape parameter
     * @param scale scale parameter
     * @return mean
     */
    public static double mean(double shape, double scale) {
        return scale * shape;
    }

    /**
     * variance of the Gamma distribution.
     *
     * @param shape shape parameter
     * @param scale scale parameter
     * @return variance
     */
    public static double variance(double shape, double scale) {
        return scale * scale * shape;
    }

    /**
     * sample from the Gamma distribution. This could be calculated using
     * quantile, but current algorithm is faster.
     *
     * @param shape shape parameter
     * @param scale scale parameter
     * @return sample
     */
    public static double nextGamma(double shape, double scale) {
        if (TRY_COLT) {
            return coltGamma.nextDouble(shape, 1.0/scale);
        }
        return nextGamma(shape, scale, false);
    }

    public static double nextGamma(double shape, double scale, boolean slowCode) {

        double sample = 0.0;

        if (shape < 0.00001) {

            if (shape < 0) {
                System.out.println("Negative shape parameter");
                throw new IllegalArgumentException("Negative shape parameter");
            }

            /*
                * special case: shape==0.0 is an improper distribution; but
                * sampling works if very small values are ignored (v. large ones
                * don't happen) This is useful e.g. for sampling from the truncated
                * Gamma(0,x)-distribution.
                */

            double minimum = 1.0e-20;
            double maximum = 50;
            double normalizingConstant = Math.log(maximum) - Math.log(minimum);
            // Draw from 1/x (with boundaries), and shape by exp(-x)
            do {
                sample = Math.exp(Math.log(minimum) + normalizingConstant
                        * MathUtils.nextDouble());
            } while (Math.exp(-sample) < MathUtils.nextDouble());
            // This distribution is actually scale-free, so multiplying by
            // 'scale' is not necessary
            return sample;
        }

        if (slowCode && Math.floor(shape) == shape && shape > 4.0) {
            for (int i = 0; i < shape; i++)
                sample += -Math.log(MathUtils.nextDouble());
            return sample * scale;
        } else {

            // Fast special cases
            if (shape == 1.0) {
                return -Math.log(MathUtils.nextDouble()) * scale;
            }
            if (shape == 2.0) {
                return -Math.log(MathUtils.nextDouble()
                        * MathUtils.nextDouble())
                        * scale;
            }
            if (shape == 3.0) {
                return -Math.log(MathUtils.nextDouble()
                        * MathUtils.nextDouble() * MathUtils.nextDouble())
                        * scale;
            }
            if (shape == 4.0) {
                return -Math.log(MathUtils.nextDouble()
                        * MathUtils.nextDouble() * MathUtils.nextDouble()
                        * MathUtils.nextDouble())
                        * scale;
            }
        }

        // general case
        do {
            try {
                sample = quantile(MathUtils.nextDouble(), shape, scale);
            } catch (IllegalArgumentException e) {
                // random doubles do go outside the permissible range 0.000002 <
                // q < 0.999998
                sample = 0.0;
            }
        } while (sample == 0.0);
        return sample;
    }

    /**
     * Sample from the gamma distribution, modified by a factor exp(
     * -(x*bias)^-1 ), i.e. from x^(shape - 1) exp(-x/scale) exp(-1/(bias*x))
     * <p/>
     * Works by rejection sampling, using a shifted ordinary Gamma as proposal
     * distribution.
     * <p/>
     * (For an older, less efficient algorithm that breaks down for shape <= 1.0,
     * see revision 287)
     */
    public static double nextExpGamma(double shape, double scale, double bias) {
        return nextExpGamma(shape, scale, bias, false);
    }

    public static double nextExpGamma(double shape, double scale, double bias,
                                      boolean slowCode) {

        double sample;
        double accept;
        int iters = 0;

        if (slowCode) {

            // for testing purposes -- this can get stuck completely for small bias parameters
            do {
                sample = nextGamma(shape, scale);
                accept = Math.exp(-1.0 / (bias * sample));
            } while (MathUtils.nextDouble() > accept);

        } else {

            if (shape < 0) {

                //return scale / (nextExpGamma(-shape, scale, bias) * bias);
                return 1.0 / nextExpGamma(-shape, bias, scale);

            }

            if (shape == 0) {

                // sample from the restriction to x >= median, and sample the other half by using
                // the self-similarity transformation x -> scale / (bias*x)

                double median = Math.sqrt(scale / bias);
                double rejection_mode = 1.0 / bias;      // mode of rejection distribution, x^-1 e^(-1/bias x)
                if (rejection_mode < median)
                    rejection_mode = median;
                double rejection_norm = (1.0 / rejection_mode) * Math.exp(-1.0 / (bias * rejection_mode));
                do {
                    sample = nextGamma(1.0, scale) + median;
                    accept = (1.0 / sample) * Math.exp(-1.0 / (sample * bias)) / rejection_norm;
                    iters += 1;
                } while (MathUtils.nextDouble() > accept && iters < 10000);
                if (iters == 10000) {
                    System.out.println("Severe Warning: nextExpGamma (shape=0) failed to generate a sample - returning bogus value!");
                }
                if (MathUtils.nextDouble() > 0.5) {
                    sample = scale / (bias * sample);
                }
                return sample;
            }

            // Current algorithm works for all shape parameters > 0, but it becomes very inefficient
            // for v. small shape parameters.
            if (shape <= 0.0) {
                System.out.println("nextExpGamma: Illegal argument (shape parameter is must be positive)");
                throw new IllegalArgumentException("");
            }

            // the function -1/(bias*x) is bounded by x/(bias x0^2) + C, with C = -2/(bias*x0), so that these functions
            // coincide precisely when x=x0. This gives the scale parameter of the majorating Gamma distribution
            //
            // First calculate the value that maximizes the acceptance probability at the mean of the proposal distribution
            double x0 = (shape * scale + Math.sqrt(4 * scale / bias + shape * shape * scale * scale)) / 2.0;

            // calculate the scale parameter of the majorating Gamma distribution
            double majorandScale = 1.0 / ((1.0 / scale) - 1.0 / (bias * x0 * x0));

            // now do rejection sampling
            do {
                sample = nextGamma(shape, majorandScale);
                accept = Math.exp(-(sample / x0 - 1) * (sample / x0 - 1) / (bias * sample));
                iters += 1;
            } while (MathUtils.nextDouble() > accept && iters < 10000);

            if (accept > 1.0) {
                System.out.println("PROBLEM!!  This should be impossible!!  Contact the authors.");
            }
            if (majorandScale < 0.0) {
                System.out.println("PROBLEM!! This should be impossible too!!  Contact the authors.");
            }
            if (iters == 10000) {
                System.out.println("Severe Warning: nextExpGamma failed to generate a sample - returning bogus value!");
            }
        }

        return sample;

    }

    // Private

    private static double pointChi2(double prob, double v) {
        // Returns z so that Prob{x<z}=prob where x is Chi2 distributed with df
        // = v
        // RATNEST FORTRAN by
        // Best DJ & Roberts DE (1975) The percentage points of the
        // Chi2 distribution. Applied Statistics 24: 385-388. (AS91)

        final double e = 0.5e-6, aa = 0.6931471805, p = prob;
        double ch, a, q, p1, p2, t, x, b, s1, s2, s3, s4, s5, s6;
        double epsi = .01;
        if( p < 0.000002 || p > 1 - 0.000002)  {
            epsi = .000001;
        }
        // if (p < 0.000002 || p > 0.999998 || v <= 0) {
        //      throw new IllegalArgumentException("Arguments out of range p" + p + " v " + v);
        //  }
        double g = GammaFunction.lnGamma(v / 2);
        double xx = v / 2;
        double c = xx - 1;
        if (v < -1.24 * Math.log(p)) {
            ch = Math.pow((p * xx * Math.exp(g + xx * aa)), 1 / xx);
            if (ch - e < 0) {
                return ch;
            }
        } else {
            if (v > 0.32) {
                x = NormalDistribution.quantile(p, 0, 1);
                p1 = 0.222222 / v;
                ch = v * Math.pow((x * Math.sqrt(p1) + 1 - p1), 3.0);
                if (ch > 2.2 * v + 6) {
                    ch = -2 * (Math.log(1 - p) - c * Math.log(.5 * ch) + g);
                }
            } else {
                ch = 0.4;
                a = Math.log(1 - p);


                do {
                    q = ch;
                    p1 = 1 + ch * (4.67 + ch);
                    p2 = ch * (6.73 + ch * (6.66 + ch));
                    t = -0.5 + (4.67 + 2 * ch) / p1
                            - (6.73 + ch * (13.32 + 3 * ch)) / p2;
                    ch -= (1 - Math.exp(a + g + .5 * ch + c * aa) * p2 / p1)
                            / t;
                } while (Math.abs(q / ch - 1) - epsi > 0);
            }
        }
        do {
            q = ch;
            p1 = 0.5 * ch;
            if ((t = GammaFunction.incompleteGammaP(xx, p1, g)) < 0) {
                throw new IllegalArgumentException(
                        "Arguments out of range: t < 0");
            }
            p2 = p - t;
            t = p2 * Math.exp(xx * aa + g + p1 - c * Math.log(ch));
            b = t / ch;
            a = 0.5 * t - b * c;

            s1 = (210 + a * (140 + a * (105 + a * (84 + a * (70 + 60 * a))))) / 420;
            s2 = (420 + a * (735 + a * (966 + a * (1141 + 1278 * a)))) / 2520;
            s3 = (210 + a * (462 + a * (707 + 932 * a))) / 2520;
            s4 = (252 + a * (672 + 1182 * a) + c * (294 + a * (889 + 1740 * a))) / 5040;
            s5 = (84 + 264 * a + c * (175 + 606 * a)) / 2520;
            s6 = (120 + c * (346 + 127 * c)) / 5040;
            ch += t
                    * (1 + 0.5 * t * s1 - b
                    * c
                    * (s1 - b
                    * (s2 - b
                    * (s3 - b
                    * (s4 - b * (s5 - b * s6))))));
        } while (Math.abs(q / ch - 1) > e);

        return (ch);
    }

    public static void main(String[] args) {

        testQuantile(1e-10, 0.878328435043444, 0.0013696236839573005);
        testQuantile(0.5, 0.878328435043444, 0.0013696236839573005);
        testQuantile(1.0 - 1e-10, 0.878328435043444, 0.0013696236839573005);

        testQuantileCM(1e-10, 0.878328435043444, 0.0013696236839573005);
        testQuantileCM(0.5, 0.878328435043444, 0.0013696236839573005);
        testQuantileCM(1.0 - 1e-10, 0.878328435043444, 0.0013696236839573005);

        for (double i = 0.0125; i < 1.0; i += 0.025) {
        	System.out.print(i + ": ");
        	try {
        		System.out.println(new GammaDistributionImpl(0.878328435043444, 0.0013696236839573005).inverseCumulativeProbability(i));
        	} catch (MathException e) {
        	System.out.println(e.getMessage());
        	}
        }

        GammaDistribution gamma = new GammaDistribution(0.01,100.0);
        double[] samples = new double[100000];
        double sum = 0.0;
        for (int i = 0; i < samples.length; i++) {
            samples[i] = gamma.nextGamma();
            sum += samples[i];
        }
        double mean = sum/(double)samples.length;
        System.out.println("Mean = " + mean);
        double variance = 0.0;
        for (int i = 0; i < samples.length; i++) {
            variance += Math.pow((samples[i] - mean),2.0);
        }
        variance = variance/(double)samples.length;
        System.out.println("Variance = " + variance);

//        System.out
//                .println("K-S critical values: 1.22(10%), 1.36(5%), 1.63(1%)\n");
//
//        int iters = 30000;
//        testExpGamma(1.0, 0.01, 7, iters);
//        testExpGamma(1.0, 0.01, 5, iters);
//        testExpGamma(2.0, 0.01, 10000, iters);
//        testExpGamma(1.0, 0.01, 10000, iters);
//        testExpGamma(0.1, 0.01, 10000, iters);
//        testExpGamma(0.01, 0.01, 10000, iters);
//
//        testExpGamma(2.0, 0.01, 10, iters);
//        testExpGamma(1.5, 0.01, 10, iters);
//        testExpGamma(1.0, 0.01, 10, iters);
//        testExpGamma(0.9, 0.01, 10, iters);
//        testExpGamma(0.5, 0.01, 10, iters);
//        testExpGamma(0.4, 0.01, 10, iters);
//        testExpGamma(0.3, 0.01, 10, iters);
//        testExpGamma(0.2, 0.01, 10, iters);
//        testExpGamma(0.1, 0.01, 10, iters);
//
//        // test distributions with severe bias, where rejection sampling doesn't
//        // work anymore
//        testExpGamma2(2.0, 0.01, 1, iters, 0.112946);
//        testExpGamma2(2.0, 0.01, 0.1, iters, 0.328874);
//        testExpGamma2(2.0, 0.01, 0.01, iters, 1.01255);
//        testExpGamma2(1.0, 0.01, 0.0003, iters, 5.781);
//        testExpGamma2(4.0, 0.01, 0.0003, iters, 5.79604);
//        testExpGamma2(20.0, 0.01, 0.0003, iters, 5.87687);
//        testExpGamma2(10.0, 0.01, 0.01, iters, 1.05374);
//        testExpGamma2(1.0, 0.01, 0.05, iters, 0.454734);
//        // test the basic Gamma distribution
//        test(1.0, 1.0, iters);
//        test(2.0, 1.0, iters);
//        test(3.0, 1.0, iters);
//        test(4.0, 1.0, iters);
//        test(100.0, 1.0, iters);
//        testAddition(0.5, 1.0, 2, iters);
//        testAddition(0.25, 1.0, 4, iters);
//        testAddition(0.1, 1.0, 10, iters);
//        testAddition(10, 1.0, 10, iters);
//        testAddition(20, 1.0, 10, iters);
//        test(0.001, 1.0, iters);
//        test(1.0, 2.0, iters);
//        test(10.0, 1.0, iters);
//        test(16.0, 1.0, iters);
//        test(16.0, 0.1, iters);
//        test(100.0, 1.0, iters);
//        test(0.5, 1.0, iters);
//        test(0.5, 0.1, iters);
//        test(0.1, 1.0, iters);
//        test(0.9, 1.0, iters);
//        // test distributions with milder biases, and compare with results from
//        // simple rejection sampling
//        testExpGamma(2.0, 0.000001, 1000000, iters);
//        testExpGamma(2.0, 0.000001, 100000, iters);
//        testExpGamma(2.0, 0.000001, 70000, iters);
//        testExpGamma(10.0, 0.01, 7, iters);
//        testExpGamma(10.0, 0.01, 5, iters);
//        testExpGamma(1.0, 0.01, 100, iters);
//        testExpGamma(1.0, 0.01, 10, iters);
//        testExpGamma(1.0, 0.01, 7, iters / 3);
//        testExpGamma(1.0, 0.01, 5, iters / 3);
//        testExpGamma(1.0, 0.00001, 1000000, iters);
//        testExpGamma(1.0, 0.00001, 100000, iters);
//        testExpGamma(1.0, 0.00001, 10000, iters);
//        testExpGamma(1.0, 0.00001, 5000, iters / 3); /*
//														 * this one takes some
//														 * time
//														 */
//        testExpGamma(2.0, 1.0, 0.5, iters);
//        testExpGamma(2.0, 1.0, 1.0, iters);
//        testExpGamma(2.0, 1.0, 2.0, iters);
//        testExpGamma(3.0, 3.0, 2.0, iters);
//        testExpGamma(10.0, 3.0, 5.0, iters);
//        testExpGamma(1.0, 3.0, 5.0, iters);
//        testExpGamma(1.0, 10.0, 5.0, iters);
//        testExpGamma(2.0, 10.0, 5.0, iters);
//        // test the basic Gamma distribution
//        test(1.0, 1.0, iters);
//        test(2.0, 1.0, iters);
//        test(3.0, 1.0, iters);
//        test(4.0, 1.0, iters);
//        test(100.0, 1.0, iters);
//        testAddition(0.5, 1.0, 2, iters);
//        testAddition(0.25, 1.0, 4, iters);
//        testAddition(0.1, 1.0, 10, iters);
//        testAddition(10, 1.0, 10, iters);
//        testAddition(20, 1.0, 10, iters);
//        test(0.001, 1.0, iters);
//        test(1.0, 2.0, iters);
//        test(10.0, 1.0, iters);
//        test(16.0, 1.0, iters);
//        test(16.0, 0.1, iters);
//        test(100.0, 1.0, iters);
//        test(0.5, 1.0, iters);
//        test(0.5, 0.1, iters);
//        test(0.1, 1.0, iters);
//        test(0.9, 1.0, iters);

    }

    private static void testQuantile(double y, double shape, double scale) {

        long time = System.currentTimeMillis();

        double value = 0;
        for (int i = 0; i < 1000; i++) {
            value = quantile(y, shape, scale);
        }
        value = quantile(y, shape, scale);
        long elapsed = System.currentTimeMillis() - time;


        System.out.println("Quantile, "+ y +", for shape=" + shape + ", scale=" + scale
                + " : " + value + ", time=" + elapsed + "ms");

    }
    private static void testQuantileCM(double y, double shape, double scale) {

        long time = System.currentTimeMillis();

        double value = 0;
        try {
            for (int i = 0; i < 1000; i++) {
                value = (new org.apache.commons.math.distribution.GammaDistributionImpl(shape, scale)).inverseCumulativeProbability(y);
            }
            value = (new org.apache.commons.math.distribution.GammaDistributionImpl(shape, scale)).inverseCumulativeProbability(y);
        } catch (MathException e) {
            e.printStackTrace();
        }
        long elapsed = System.currentTimeMillis() - time;


        System.out.println("commons.maths inverseCDF, "+ y +", for shape=" + shape + ", scale=" + scale
                + " : " + value + ", time=" + elapsed + "ms");

    }

    /* assumes the two arrays are the same size */
    private static double KolmogorovSmirnov(List<Double> l1, List<Double> l2) {

        int idx2 = 0;
        int max = 0;
        for (int i = 0; i < l1.size(); i++) {
            while (idx2 < l2.size() && l2.get(idx2) < l1.get(i) ) {
                idx2++;
            }
            max = Math.max(max, idx2 - i);
        }
        return max / Math.sqrt(2.0 * l1.size());
    }

    private static void testExpGamma2(double shape, double scale, double bias,
                                      int iterations, double mean) {

        double s0 = 0;
        double s1 = 0;
        double s2 = 0;
        List<Double> fast = new ArrayList<Double>(0);

        for (int i = 0; i < iterations; i++) {
            double sample = nextExpGamma(shape, scale, bias, false);
            s0 += 1;
            s1 += sample;
            s2 += sample * sample;
            fast.add(sample);
        }
        Collections.sort(fast);
        double expmean = s1 / s0;
        double expvar = (s2 - s1 * s1 / s0) / s0;
        double z = (mean - expmean) / Math.sqrt(expvar / iterations);
        System.out.println("Equal-mean test: (shape=" + shape + " scale="
                + scale + " bias=" + bias + " mean=" + expmean + " expected="
                + mean + " var=" + expvar + " median="
                + fast.get(iterations / 2) + "): z=" + z);
    }

    private static void testExpGamma(double shape, double scale, double bias,
                                     int iterations) {

        List<Double> slow = new ArrayList<Double>(0);
        List<Double> fast = new ArrayList<Double>(0);

        long time = System.currentTimeMillis();
        for (int i = 0; i < iterations; i++) {
            slow.add(nextExpGamma(shape, scale, bias, true));
        }
        long slowtime = System.currentTimeMillis();
        for (int i = 0; i < iterations; i++) {
            fast.add(nextExpGamma(shape, scale, bias, false));
        }
        long fasttime = System.currentTimeMillis() - slowtime;
        slowtime -= time;

        Collections.sort(slow);
        Collections.sort(fast);

        System.out.println("KS test for shape=" + shape + ", bias=" + bias
                + " : " + KolmogorovSmirnov(slow, fast) + " and "
                + KolmogorovSmirnov(fast, slow) + " slow=" + slowtime + "ms, fast=" + fasttime + "ms");

    }

    private static void test(double shape, double scale, int iterations) {

        List<Double> slow = new ArrayList<Double>(0);
        List<Double> fast = new ArrayList<Double>(0);

        for (int i = 0; i < iterations; i++) {
            slow.add(nextGamma(shape, scale, true));
            fast.add(nextGamma(shape, scale, false));
        }

        Collections.sort(slow);
        Collections.sort(fast);

        System.out.println("KS test for shape=" + shape + " : "
                + KolmogorovSmirnov(slow, fast) + " and "
                + KolmogorovSmirnov(fast, slow));

    }

    private static void testAddition(double shape, double scale, int N,
                                     int iterations) {

        List<Double> slow = new ArrayList<Double>(0);
        List<Double> fast = new ArrayList<Double>(0);
        List<Double> test = new ArrayList<Double>(0);

        for (int i = 0; i < iterations; i++) {
            double s = 0.0;
            for (int j = 0; j < N; j++) {
                s += nextGamma(shape, scale, true);
            }
            slow.add(s);
            s = 0.0;
            for (int j = 0; j < N; j++) {
                s += nextGamma(shape, scale, false);
            }
            fast.add(s);
            test.add(nextGamma(shape * N, scale, true));
        }

        Collections.sort(slow);
        Collections.sort(fast);
        Collections.sort(test);

        System.out.println("KS test for shape=" + shape + " : slow="
                + KolmogorovSmirnov(slow, test) + " & "
                + KolmogorovSmirnov(test, slow) + "; fast="
                + KolmogorovSmirnov(fast, test) + " & "
                + KolmogorovSmirnov(test, fast));

    }

    protected double shape, scale;
    protected int samples;

    private static final boolean TRY_COLT = false;
    private static RandomEngine randomEngine;
    private static Gamma coltGamma;

}

