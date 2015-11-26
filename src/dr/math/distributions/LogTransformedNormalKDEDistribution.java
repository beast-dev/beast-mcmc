/*
 * LogTransformedNormalKDEDistribution.java
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

import dr.math.ComplexArray;
import dr.math.FastFourierTransform;
import dr.stats.DiscreteStatistics;
import dr.util.HeapSort;

import java.util.Arrays;

/**
 * @author Guy Baele
 */
public class LogTransformedNormalKDEDistribution extends KernelDensityEstimatorDistribution {

    public static final int MINIMUM_GRID_SIZE = 2048;
    public static final boolean DEBUG = false;

    //the samples should not already be log transformed (the log transformation is done in this class)
    public LogTransformedNormalKDEDistribution(Double[] sample) {
    	this(sample, null, null, null);
    }
    
    public LogTransformedNormalKDEDistribution(Double[] sample, int n) {
    	this(sample, null, null, null, 3.0, n);
    }
    
    public LogTransformedNormalKDEDistribution(Double[] sample, Double lowerBound, Double upperBound, Double bandWidth) {
        this(sample, lowerBound, upperBound, bandWidth, 3.0, MINIMUM_GRID_SIZE);
    }

    public LogTransformedNormalKDEDistribution(Double[] sample, Double lowerBound, Double upperBound, Double bandWidth,
                                 int n) {
        this(sample, lowerBound, upperBound, bandWidth, 3.0, n);
    }

    public LogTransformedNormalKDEDistribution(Double[] sample, Double lowerBound, Double upperBound, Double bandWidth, double cut, int n) {
    	
    	//first call the super constructor, but immediately overwrite the stored information
    	/* code in super constructor
    	this.sample = new double[sample.length];
        for (int i = 0; i < sample.length; i++) {
            this.sample[i] = sample[i];
        }
        this.N = sample.length;
        processBounds(lowerBound, upperBound);
        setBandWidth(bandWidth); 
    	*/

        super(sample, lowerBound, upperBound, bandWidth);
        //transform the data to the log scale and store in logSample
        if (DEBUG) {
            System.out.println("Creating the KDE in log space");
            System.out.println("lowerBound = " + lowerBound);
            System.out.println("upperBound = " + upperBound);
        }

        this.logSample = new double[sample.length];
        for (int i = 0; i < logSample.length; i++) {
        	this.logSample[i] = Math.log(sample[i]);
        }
        //keep a backup copy of the samples in normal space
        this.backupSample = new double[sample.length];
        for (int i = 0; i < sample.length; i++) {
        	this.backupSample[i] = sample[i];
        }
        //overwrite the stored samples, sample.length stays the same
        this.sample = logSample;
        processBounds(lowerBound, upperBound);
        setBandWidth(bandWidth);
        
        this.gridSize = Math.max(n, MINIMUM_GRID_SIZE);
        if (this.gridSize > MINIMUM_GRID_SIZE) {
            this.gridSize = (int) Math.pow(2, Math.ceil(Math.log(this.gridSize) / Math.log(2.0)));
        }
        this.cut = cut;

        from = DiscreteStatistics.min(this.sample) - this.cut * this.bandWidth;
        to = DiscreteStatistics.max(this.sample) + this.cut * this.bandWidth;

        if (DEBUG) {
            System.out.println("bandWidth = " + this.bandWidth);
            System.out.println("cut = " + this.cut);
            System.out.println("from = " + from);
            System.out.println("to = " + to);
        }
        
        lo = from - 4.0 * this.bandWidth;
        up = to + 4.0 * this.bandWidth;

        if (DEBUG) {
            System.out.println("lo = " + lo);
            System.out.println("up = " + up);
        }

        densityKnown = false;
        
        //run computeDensity to estimate the KDE on the log scale
        //and afterwards return to the normal scale
        computeDensity();
        
    }
    
    public double getFromPoint() {
        return from;
    }

    public double getToPoint() {
        return to;
    }

    /**
     * Returns a linear approximation evaluated at pt
     * @param x data (assumed sorted increasingly
     * @param y data
     * @param pt evaluation point
     * @param low return value if pt < x
     * @param high return value if pt > x
     * @return  evaluated coordinate
     */
    private double linearApproximate(double[] x, double[] y, double pt, double low, double high) {

        int i = 0;
        int j = x.length - 1;

        if (pt < x[i]) {
            return low;
        }
        if (pt > x[j]) {
            return high;
        }

        // Bisection search
        while (i < j - 1) {
            int ij = (i + j) / 2;
            if (pt < x[ij]) {
                j = ij;
            } else {
                i = ij;
            }
        }

        if (pt == x[j]) {
            return y[j];
        }
        if (pt == x[i]) {
            return y[i];
        }
        //System.out.println("return value: "+ (y[i] + (y[j] - y[i]) * ((pt - x[i]) / (x[j] - x[i]))));
        return y[i] + (y[j] - y[i]) * ((pt - x[i]) / (x[j] - x[i]));
    }

    private double[] rescaleAndTrim(double[] x) {
        final int length = x.length / 2;
        final double scale = 1.0 / x.length;
        double[] out = new double[length];
        for (int i = 0; i < length; ++i) {
            out[i] = x[i] * scale;
            if (out[i] < 0) {
                out[i] = 0;
            }
        }
        return out;
    }

    private double[] massdist(double[] x, double xlow, double xhigh, int ny) {

        int nx = x.length;
        double[] y = new double[ny * 2];

        final int ixmin = 0;
        final int ixmax = ny - 2;
        final double xdelta = (xhigh - xlow) / (ny - 1);

        for (int i = 0; i < ny; ++i) {
            y[i] = 0.0;
        }

        final double xmi = 1.0 / nx;
        for (int i = 0; i < nx; ++i) {
            final double xpos = (x[i] - xlow) /  xdelta;
            final int ix = (int) Math.floor(xpos);
            final double fx = xpos - ix;
//            final double xmi = xmass[i];

            if (ixmin <= ix && ix <= ixmax) {
                y[ix] += (1 - fx) * xmi;
                y[ix + 1] += fx * xmi;
            } else if (ix == -1) {
                y[0] += fx * xmi;
            } else if (ix == ixmax + 1) {
                y[ix] += (1 - fx) * xmi;
            }
        }
        return y;
    }

    /**
     * Override for different kernels
     * @param ordinates the points in complex space
     * @param bandWidth predetermined bandwidth
     */
    protected void fillKernelOrdinates(ComplexArray ordinates, double bandWidth) {
        final int length = ordinates.length;
        final double a = 1.0 / (Math.sqrt(2.0 * Math.PI) * bandWidth);
        final double precision = -0.5 / (bandWidth * bandWidth);
        for (int i = 0; i < length; i++) {
            final double x = ordinates.real[i];
            ordinates.real[i] = a * Math.exp(x * x * precision);
        }
    }

    protected void computeDensity() {
    	//transformData calls massdist and rescaleAndTrim
        makeOrdinates();
        //makeOrdinates calls fillKernelOrdinates
        transformData();
        //we're still in log space and need to return to normal space
        //preferably before setting densityKnown to true
        //stored values are in xPoints and densityPoints
        transformEstimator();
        densityKnown = true;
    }
    
    private void transformEstimator() {
    	
    	if (DEBUG) {
            System.out.println("\nCreating the KDE in normal space");
            System.out.println("lowerBound = " + lowerBound);
            System.out.println("upperBound = " + upperBound);
        }

    	this.sample = backupSample;
    	//processBounds(lowerBound, upperBound);
        setBandWidth(null);
        
        from = DiscreteStatistics.min(this.sample) - this.cut * this.bandWidth;
        to = DiscreteStatistics.max(this.sample) + this.cut * this.bandWidth;

        if (DEBUG) {
            System.out.println("min: " + DiscreteStatistics.min(this.sample));
            System.out.println("max: " + DiscreteStatistics.max(this.sample));
            System.out.println("bandWidth = " + this.bandWidth);
            System.out.println("cut = " + this.cut);
            System.out.println("from = " + from);
            System.out.println("to = " + to);
        }
        
        lo = from - 4.0 * this.bandWidth;
        up = to + 4.0 * this.bandWidth;
        
        if (DEBUG) {
            System.out.println("lo = " + lo);
            System.out.println("up = " + up);
        }

        if (lo < 0.0) {
            //small hack, but our asymmetric kernel estimators are terribly slow
            lo = DiscreteStatistics.min(this.sample);
        }
        
        //make new ordinates for the transformation back to normal space
        //need a backup of the xPoints for the log scale KDE
        this.backupXPoints = new double[xPoints.length];
        System.arraycopy(xPoints, 0, backupXPoints, 0, xPoints.length);
        makeOrdinates();
        
        //the KDE on log scale is contained in the xPoints and densityPoints arrays
    	//copy them to finalXPoints and finalDensityPoints

    	//this.finalXPoints = new double[xPoints.length - numberOfNegatives];
        //this.finalXPoints = new double[xPoints.length];

    	//System.arraycopy(xPoints, numberOfNegatives, finalXPoints, 0, xPoints.length - numberOfNegatives);
        //System.arraycopy(xPoints, numberOfNegatives, finalXPoints, 0, xPoints.length);
    	if (DEBUG) {
            for (int i = 0; i < xPoints.length; i++) {
                System.out.println(xPoints[i] + "   " + backupXPoints[i] + " : " + densityPoints[i]);
            }
            //System.out.println("\nfinalXPoints length = " + finalXPoints.length);
        }

    	//this.finalDensityPoints = new double[densityPoints.length - numberOfNegatives];
        this.finalDensityPoints = new double[densityPoints.length];
    	
    	for (int i = 0; i < xPoints.length; i++) {
    		finalDensityPoints[i] = linearApproximate(backupXPoints, densityPoints, Math.log(xPoints[i]), 0.0, 0.0)*(1.0/xPoints[i]);
    		if (DEBUG) {
                System.out.println(xPoints[i] + "\t" + finalDensityPoints[i]);
            }
    	}
    	
    	//System.exit(0);
    	
    }

    private void transformData() {
        ComplexArray Y  = new ComplexArray(massdist(this.logSample, lo, up, this.gridSize));
        FastFourierTransform.fft(Y, false);

        ComplexArray product = Y.product(kOrdinates);
        FastFourierTransform.fft(product, true);

        densityPoints = rescaleAndTrim(product.real); 
        /*System.out.println("Y");
        for (int i = 0; i < gridSize; i++) {
        	System.out.println(densityPoints[i]);
        }*/
    }

    private void makeOrdinates() {

        final int length = 2 * gridSize;
        if (kOrdinates == null) {
            kOrdinates = new ComplexArray(new double[length]);
        }

        // Fill with grid values
        final double max = 2.0 * (up - lo);
        double value = 0;
        final double inc = max / (length - 1);
        for (int i = 0; i <= gridSize; i++) {
            kOrdinates.real[i] = value;
            value += inc;
        }
        for (int i = gridSize + 1; i < length; i++) {
            kOrdinates.real[i] = -kOrdinates.real[length - i];
        }
        fillKernelOrdinates(kOrdinates, bandWidth);

        FastFourierTransform.fft(kOrdinates, false);
        kOrdinates.conjugate();

        // Make x grid
        xPoints = new double[gridSize];
        double x = lo;
        double delta = (up - lo) / (gridSize - 1);
        if (DEBUG) {
            System.out.println("X");
        }
        for (int i = 0; i < gridSize; i++) {
            xPoints[i] = x;
            x += delta;
            if (DEBUG) {
                System.out.println(xPoints[i]);
            }
        }
    }

    @Override
    protected double evaluateKernel(double x) {        
        if (!densityKnown) {
        	//computeDensity() calls makeOrdinates and transformData
        	computeDensity();
        }
        //xPoints and densityPoints are now back in normal space
        return linearApproximate(xPoints, finalDensityPoints, x, 0.0, 0.0);
    }

    @Override
    protected void processBounds(Double lowerBound, Double upperBound) {
        if ((lowerBound != null && lowerBound != Double.NEGATIVE_INFINITY) ||
                (upperBound != null && upperBound != Double.POSITIVE_INFINITY)) {
            throw new RuntimeException("LogTransformedNormalKDEDistribution must be unbounded");
        }
    }

    @Override
    protected void setBandWidth(Double bandWidth) {
        if (bandWidth == null) {
            // Default bandwidth
            this.bandWidth = bandwidthNRD(sample);
        } else
            this.bandWidth = bandWidth;
                    
        densityKnown = false;
    }

    public double bandwidthNRD(double[] x) {
        int[] indices = new int[x.length];
        HeapSort.sort(x, indices);

        final double h =
                (DiscreteStatistics.quantile(0.75, x, indices) - DiscreteStatistics.quantile(0.25, x, indices)) / 1.34;
        return 1.06 *
                Math.min(Math.sqrt(DiscreteStatistics.variance(x)), h) *
                Math.pow(x.length, -0.2);
    }

    private ComplexArray kOrdinates;
    private double[] xPoints, backupXPoints;
    private double[] densityPoints, finalDensityPoints;
    private double[] backupSample, logSample;

    private int gridSize;
    private double cut;
    private double from;
    private double to;
    private double lo;
    private double up;

    private boolean densityKnown = false;

    public static void main(String[] args) {

        /*Double[] testfivehundred = {9.169378e-06,9.169378e-06,9.169378e-06,9.169378e-06,4.116020e-06,4.116020e-06,4.116020e-06,4.161038e-07,3.644982e-07,2.195354e-07,2.195354e-07,2.195354e-07,2.195354e-07,7.462043e-08,7.462043e-08,8.368972e-08,9.283378e-08,9.283378e-08,2.866075e-08,6.532907e-09,6.532907e-09,3.936993e-08,3.936993e-08,2.257473e-08,2.257473e-08,1.525075e-08,1.525075e-08,1.525075e-08,1.525075e-08,1.525075e-08,2.370723e-08,2.370723e-08,1.961463e-08,7.642434e-08,8.838916e-08,8.838916e-08,8.838916e-08,8.838916e-08,1.149915e-07,1.149915e-07,7.226055e-08,7.226055e-08,6.276125e-08,5.330820e-08,5.330820e-08,4.626315e-08,4.626315e-08,4.709987e-08,1.011533e-07,1.011533e-07,1.011533e-07,1.011533e-07,9.176941e-08,1.876641e-08,1.876641e-08,1.582991e-08,1.582991e-08,6.121650e-08,5.823668e-08,5.823668e-08,5.823668e-08,5.823668e-08,5.823668e-08,1.138725e-08,1.327368e-08,1.327368e-08,1.327368e-08,5.123223e-08,5.123223e-08,5.123223e-08,3.812422e-08,2.592858e-08,1.584334e-07,1.584334e-07,1.584334e-07,4.091027e-07,4.091027e-07,4.091027e-07,3.180726e-07,5.671995e-07,5.671995e-07,5.671995e-07,5.671995e-07,3.718532e-07,2.140636e-07,1.141818e-07,9.974324e-08,8.885531e-08,8.885531e-08,1.336255e-08,1.336255e-08,1.336255e-08,2.066815e-09,1.359399e-09,2.433474e-09,2.433474e-09,2.621958e-09,2.621958e-09,2.621958e-09,2.621958e-09,1.760426e-09,5.016878e-09,1.297173e-09,1.297173e-09,2.729714e-09,2.729714e-09,4.096775e-09,4.066638e-09,4.938423e-09,2.279718e-09,2.279718e-09,2.279718e-09,2.279718e-09,7.927200e-09,7.927200e-09,5.882404e-09,3.039284e-09,3.039284e-09,7.933658e-09,1.577847e-08,1.577847e-08,1.577847e-08,1.577847e-08,7.421679e-09,9.966059e-09,9.966059e-09,2.978775e-08,2.978775e-08,2.978775e-08,2.037757e-08,2.037757e-08,3.838478e-08,6.434436e-09,6.434436e-09,6.434436e-09,1.330451e-09,1.330451e-09,1.330451e-09,1.330451e-09,1.330451e-09,1.821259e-09,1.821259e-09,1.821259e-09,1.821259e-09,1.136283e-09,1.136283e-09,1.136283e-09,8.538424e-09,1.115792e-08,1.115792e-08,1.115792e-08,6.708989e-09,6.708989e-09,6.708989e-09,6.708989e-09,6.708989e-09,6.171120e-09,6.171120e-09,6.171120e-09,4.324598e-09,4.324598e-09,4.086923e-09,1.018727e-08,1.076947e-08,7.300160e-09,1.882574e-08,7.825507e-09,5.768726e-09,5.768726e-09,5.768726e-09,3.647371e-09,5.801127e-09,5.801127e-09,5.801127e-09,1.061525e-08,1.232341e-08,1.232341e-08,1.232341e-08,1.232341e-08,3.058357e-09,3.058357e-09,3.058357e-09,4.222149e-09,4.222149e-09,4.222149e-09,9.101976e-09,9.101976e-09,9.101976e-09,9.101976e-09,1.272342e-08,5.057809e-08,8.124001e-08,2.988208e-07,2.988208e-07,2.208343e-06,2.208343e-06,2.563228e-06,7.319180e-06,5.100615e-06,5.100615e-06,3.908437e-05,9.561673e-06,9.561673e-06,4.051668e-07,1.165789e-06,1.165789e-06,6.355679e-07,1.267624e-06,1.267624e-06,1.269425e-06,1.269425e-06,6.906194e-07,3.027074e-06,4.954935e-06,4.954935e-06,3.069373e-06,3.516738e-05,3.516738e-05,1.337480e-04,2.621304e-04,2.621304e-04,2.621304e-04,2.621304e-04,2.621304e-04,8.714378e-04,8.714378e-04,8.714378e-04,8.714378e-04,8.714378e-04,8.714378e-04,8.714378e-04,8.714378e-04,3.789225e-03,2.911953e-03,2.911953e-03,2.911953e-03,2.911953e-03,2.911953e-03,2.911953e-03,2.889193e-03,2.741145e-03,2.741145e-03,2.741145e-03,2.741145e-03,1.039219e-03,1.039219e-03,1.039219e-03,4.519776e-04,5.756779e-04,5.756779e-04,1.608451e-04,2.319337e-04,1.004998e-04,1.004998e-04,1.004998e-04,2.479749e-04,2.479749e-04,2.479749e-04,2.290535e-04,2.290535e-04,2.290535e-04,5.562921e-05,5.562921e-05,5.562921e-05,2.263067e-05,2.924862e-05,1.625127e-05,1.625127e-05,1.625127e-05,5.130025e-06,5.130025e-06,5.130025e-06,2.601581e-06,3.884423e-05,3.884423e-05,3.884423e-05,1.899082e-05,1.128966e-05,1.128966e-05,1.128966e-05,1.457212e-05,8.006272e-06,8.006272e-06,1.007007e-05,1.007007e-05,1.007007e-05,1.007007e-05,8.494570e-06,8.494570e-06,8.494570e-06,1.394617e-05,1.394617e-05,1.394617e-05,7.943013e-06,1.523964e-05,1.523964e-05,2.775473e-05,2.775473e-05,2.775473e-05,2.775473e-05,2.775473e-05,3.772894e-05,3.637147e-05,7.268067e-06,3.481332e-06,3.481332e-06,2.401154e-06,2.401154e-06,8.197940e-06,8.197940e-06,3.997811e-06,3.997811e-06,3.997811e-06,3.997811e-06,1.207494e-05,1.207494e-05,1.963691e-05,2.055311e-05,2.178261e-06,3.398210e-06,3.398210e-06,3.175174e-06,3.175174e-06,1.975710e-05,5.860156e-05,5.860156e-05,5.860156e-05,7.028432e-05,7.028432e-05,7.028432e-05,7.028432e-05,7.028432e-05,7.028432e-05,8.770771e-05,8.770771e-05,3.941266e-04,3.941266e-04,1.002434e-04,4.501147e-04,9.140456e-05,9.140456e-05,9.140456e-05,1.216764e-04,1.216764e-04,1.216764e-04,1.216764e-04,1.216764e-04,1.358293e-04,1.358293e-04,1.358293e-04,1.358293e-04,2.433390e-04,2.128913e-04,6.047829e-05,6.047829e-05,6.047829e-05,6.047829e-05,6.047829e-05,6.047829e-05,1.559573e-04,6.756065e-04,6.756065e-04,1.252698e-03,1.252698e-03,1.252698e-03,1.722452e-03,2.974541e-03,2.788099e-03,2.788099e-03,2.788099e-03,2.788099e-03,2.788099e-03,2.788099e-03,9.704677e-04,6.524029e-04,6.524029e-04,6.524029e-04,6.524029e-04,1.986166e-04,5.148064e-04,3.935172e-05,3.935172e-05,3.935172e-05,3.596480e-05,3.596480e-05,3.840476e-04,3.265822e-04,3.265822e-04,1.699738e-03,3.382742e-04,3.382742e-04,4.291126e-04,4.291126e-04,4.291126e-04,4.291126e-04,3.012033e-04,7.874703e-05,4.009423e-04,4.009423e-04,4.009423e-04,4.009423e-04,7.015827e-04,7.015827e-04,8.998811e-04,3.750951e-04,8.001801e-05,8.001801e-05,8.001801e-05,8.001801e-05,8.769844e-05,8.769844e-05,2.336436e-05,8.662544e-05,4.360443e-05,4.360443e-05,4.360443e-05,4.360443e-05,4.360443e-05,4.360443e-05,3.684793e-05,1.296607e-05,2.434277e-06,4.454218e-06,4.454218e-06,5.346013e-06,5.346013e-06,5.346013e-06,9.166112e-07,9.166112e-07,9.166112e-07,2.914328e-07,4.793011e-08,4.793011e-08,2.518825e-08,2.518825e-08,2.559798e-08,2.559798e-08,2.559798e-08,2.559798e-08,1.786016e-08,1.786016e-08,1.786016e-08,1.786016e-08,1.786016e-08,2.581946e-08,2.581946e-08,7.651924e-09,7.651924e-09,7.651924e-09,2.654439e-09,2.654439e-09,1.160599e-08,1.089639e-08,4.441808e-09,2.993354e-08,2.993354e-08,6.597157e-09,1.756244e-08,3.663998e-08,3.663998e-08,3.663998e-08,3.663998e-08,2.087983e-08,7.219497e-09,7.219497e-09,7.219497e-09,8.681694e-09,8.681694e-09,8.681694e-09,8.681694e-09,8.681694e-09,2.434356e-10,5.152142e-10,1.627307e-10,1.627307e-10,1.627307e-10,3.166602e-10,3.166602e-10,3.166602e-10,3.166602e-10,1.984753e-10,1.984753e-10,6.607674e-10,6.607674e-10,7.945952e-10,7.945952e-10,1.049938e-09,1.049938e-09,1.049938e-09,1.049938e-09,1.049938e-09,1.049938e-09,1.049938e-09,7.585986e-10,7.585986e-10,7.585986e-10,8.310453e-10,3.369425e-10,3.369425e-10,3.369425e-10,3.369425e-10};
        LogTransformedNormalKDEDistribution five = new LogTransformedNormalKDEDistribution(testfivehundred);

        for (int i = 0; i < 100; i++) {
            System.out.println(((double) i / 1E7) + " : " + five.evaluateKernel((double)i/1E7));
        }

        System.exit(0);*/


        /*Double[] testvpu = {3.264474e-10,3.264474e-10,2.270974e-10,4.053443e-10,1.440543e-10,1.948030e-11,2.160727e-11,2.160727e-11,2.160727e-11,2.798423e-11,1.220722e-10,1.220722e-10,6.669941e-11,6.669941e-11,4.142276e-11,5.325259e-11,5.325259e-11,5.325259e-11,9.190992e-11,1.165668e-10,4.202000e-11,4.202000e-11,4.202000e-11,7.824317e-12,7.824317e-12,7.824317e-12,7.824317e-12,7.824317e-12,7.824317e-12,1.245222e-12,5.797409e-14,5.797409e-14,5.797409e-14,5.797409e-14,5.797409e-14,5.797409e-14,6.524378e-14,3.446118e-14,3.446118e-14,3.446118e-14,1.406773e-14,3.089455e-14,3.089455e-14,3.089455e-14,4.719081e-14,2.636496e-14,2.854360e-14,2.854360e-14,3.201851e-14,3.201851e-14,5.526584e-14,5.526584e-14,5.526584e-14,5.526584e-14,2.943399e-13,4.180499e-13,4.180499e-13,4.180499e-13,3.678937e-13,3.678937e-13,2.121552e-13,8.475159e-13,3.201119e-13,2.478820e-13,2.478820e-13,2.927839e-13,2.927839e-13,1.586398e-13,1.586398e-13,8.238769e-14,1.329486e-13,9.046881e-13,9.046881e-13,9.046881e-13,9.046881e-13,1.075469e-13,7.662459e-14,9.494873e-14,9.494873e-14,1.415190e-13,1.415190e-13,1.415190e-13,1.415190e-13,1.180839e-13,1.180839e-13,1.180839e-13,1.180839e-13,7.610445e-14,7.150592e-14,3.058209e-14,2.428279e-14,2.428279e-14,2.428279e-14,2.428279e-14,2.428279e-14,3.412281e-14,1.005870e-13,1.005870e-13,3.600161e-14,5.032761e-14,4.792479e-14,3.137045e-14,3.137045e-14,7.775396e-15,9.260730e-15,2.632486e-14,1.542887e-14,1.542887e-14,1.542887e-14,3.920010e-16,2.705095e-16,2.705095e-16,2.705095e-16,1.611363e-16,1.611363e-16,1.611363e-16,1.611363e-16,1.611363e-16,1.611363e-16,1.611363e-16,5.833794e-17,9.684154e-18,9.684154e-18,2.140219e-17,2.140219e-17,2.140219e-17,2.131286e-17,1.997783e-18,1.003044e-18,1.003044e-18,1.003044e-18,1.003044e-18,1.252530e-18,1.252530e-18,1.252530e-18,5.360969e-19,4.049063e-20,8.391021e-20,8.391021e-20,8.391021e-20,8.391021e-20,8.391021e-20,1.656486e-20,1.656486e-20,2.592875e-20,2.592875e-20,2.592875e-20,6.944018e-20,6.944018e-20,9.193503e-21,3.153529e-20,2.002192e-21,9.662736e-22,9.662736e-22,2.366541e-21,2.366541e-21,9.746925e-22,5.145195e-22,5.145195e-22,2.040818e-21,9.600368e-22,9.600368e-22,9.600368e-22,9.600368e-22,9.600368e-22,6.477488e-22,8.845651e-23,8.845651e-23,8.612247e-23,4.819019e-23,5.332775e-23,5.332775e-23,1.612997e-23,3.568781e-24,4.478014e-25,4.478014e-25,4.608863e-25,2.532690e-25,2.532690e-25,4.184380e-26,4.184380e-26,4.833258e-26,1.996298e-26,7.691123e-27,7.691123e-27,7.577914e-27,7.506433e-27,7.506433e-27,7.506433e-27,4.145722e-28,2.791080e-28,1.087618e-27,1.087618e-27,1.152689e-27,1.152689e-27,1.152689e-27,1.152689e-27,9.221644e-28,8.193212e-29,8.193212e-29,8.193212e-29,2.450405e-29,2.591711e-29,2.591711e-29,2.591711e-29,7.476938e-30,7.476938e-30,7.476938e-30,7.476938e-30,7.476938e-30,7.476938e-30,7.476938e-30,7.476938e-30,4.053158e-30,3.491727e-29,1.156184e-28,1.156184e-28,1.379535e-28,9.951936e-28,9.951936e-28,9.951936e-28,9.951936e-28,9.951936e-28,9.951936e-28,1.452129e-27,6.233493e-27,6.233493e-27,1.890821e-27,1.890821e-27,4.573558e-27,4.573558e-27,1.852501e-26,1.852501e-26,3.534462e-26,6.457303e-26,6.457303e-26,6.457303e-26,9.789591e-26,1.722045e-26,2.146792e-27,2.146792e-27,2.146792e-27,2.146792e-27,2.146792e-27,2.163014e-27,2.163014e-27,4.770303e-28,4.770303e-28,4.770303e-28,9.266854e-28,9.266854e-28,9.266854e-28,7.038816e-28,2.887066e-28,2.887066e-28,2.887066e-28,3.215952e-28,1.694957e-28,1.787078e-28,1.787078e-28,1.533689e-28,8.813121e-29,8.813121e-29,8.813121e-29,8.813121e-29,8.813121e-29,8.813121e-29,8.813121e-29,8.813121e-29,8.813121e-29,8.813121e-29,8.813121e-29,8.202488e-29,8.202488e-29,2.341551e-28,2.341551e-28,3.304727e-29,3.127826e-29,3.127826e-29,5.295062e-30,5.295062e-30,5.295062e-30,9.226137e-30,9.226137e-30,9.226137e-30,5.203193e-30,1.288022e-29,1.288022e-29,2.330187e-29,9.787334e-29,9.787334e-29,1.559790e-29,1.559790e-29,2.492078e-30,2.492078e-30,2.492078e-30,4.966093e-30,4.966093e-30,4.966093e-30,4.966093e-30,1.350492e-29,2.104433e-29,8.483873e-29,5.165878e-29,3.600579e-29,3.600579e-29,3.600579e-29,2.624000e-28,3.463846e-28,3.463846e-28,1.481569e-28,1.481569e-28,3.577404e-28,3.914783e-28,4.493977e-28,4.493977e-28,4.493977e-28,4.493977e-28,4.493977e-28,4.493977e-28,3.432940e-28,3.432940e-28,3.432940e-28,3.432940e-28,3.432940e-28,4.480180e-28,4.390359e-28,5.724157e-28,5.724157e-28,8.549744e-28,8.549744e-28,4.305970e-28,4.305970e-28,1.393338e-27,8.152519e-28,4.189267e-27,5.074637e-26,5.074637e-26,5.074637e-26,7.426017e-26,6.834243e-26,3.218294e-26,3.218294e-26,4.572465e-26,4.572465e-26,4.572465e-26,2.120649e-25,2.120649e-25,2.120649e-25,2.120649e-25,8.192616e-25,1.464765e-25,1.464765e-25,1.464765e-25,1.464765e-25,1.464765e-25,3.165377e-25,3.165377e-25,3.165377e-25,3.165377e-25,2.657489e-25,2.657489e-25,2.657489e-25,3.619027e-25,2.767848e-25,1.162370e-24,1.162370e-24,3.354422e-24,3.967748e-24,3.967748e-24,3.967748e-24,3.967748e-24,3.967748e-24,2.275277e-24,4.737946e-24,4.737946e-24,3.155432e-24,3.155432e-24,3.155432e-24,3.155432e-24,3.155432e-24,5.108613e-24,2.534356e-23,9.638765e-23,9.638765e-23,1.082050e-22,3.284043e-22,2.544953e-22,2.561433e-22,2.561433e-22,2.561433e-22,1.788474e-22,1.788474e-22,3.070120e-22,6.679008e-22,4.254458e-23,4.254458e-23,6.737668e-24,9.744693e-24,6.508964e-24,3.432448e-23,3.432448e-23,3.997224e-23,3.997224e-23,3.997224e-23,3.997224e-23,3.997224e-23,3.997224e-23,3.380456e-23,1.697625e-23,1.697625e-23,3.928130e-24,9.920672e-25,9.920672e-25,9.920672e-25,6.972963e-25,6.972963e-25,6.972963e-25,6.972963e-25,1.245301e-24,4.815021e-25,2.025882e-24,2.025882e-24,2.025882e-24,2.025882e-24,3.271768e-24,3.271768e-24,7.527654e-24,7.527654e-24,4.701377e-24,2.168193e-24,2.168193e-24,3.201564e-24,3.201564e-24,3.201564e-24,3.320860e-25,3.320860e-25,3.320860e-25,3.320860e-25,4.065828e-25,3.251797e-25,3.251797e-25,3.251797e-25,3.251797e-25,3.251797e-25,3.251797e-25,4.751037e-25,4.751037e-25,4.751037e-25,4.751037e-25,2.470807e-25,2.470807e-25,1.533374e-25,1.533374e-25,8.054821e-26,5.435195e-25,1.292397e-24,1.270255e-24,1.270255e-24,1.270255e-24,1.270255e-24,1.013316e-24,1.242894e-24,2.091064e-24,2.091064e-24,2.729537e-24,9.825103e-25,9.825103e-25,9.825103e-25,9.825103e-25,9.825103e-25,2.213826e-25,9.499452e-25,1.038643e-24,1.038643e-24,3.891512e-24,2.886182e-24,9.920124e-25,9.920124e-25,9.920124e-25,2.813333e-24,1.258997e-23,8.722277e-24,8.722277e-24,8.722277e-24,2.575166e-24,1.972378e-24,1.972378e-24,1.972378e-24,1.925702e-24,5.637562e-24,5.637562e-24,5.637562e-24,1.047603e-23,1.245956e-23,1.245956e-23,1.034809e-23,1.034809e-23,1.034809e-23,1.034809e-23,1.218613e-23};
        LogTransformedNormalKDEDistribution five = new LogTransformedNormalKDEDistribution(testvpu);

        for (int i = 0; i < 100; i++) {
            System.out.println(((double) i / 1E12) + " : " + five.evaluateKernel((double)i/1E12));
        }
        System.out.println();
        for (int i = 0; i < 100; i++) {
            System.out.println(((double) i / 1E14) + " : " + five.evaluateKernel((double)i/1E14));
        }
        System.out.println();
        for (int i = 0; i < 100; i++) {
            System.out.println(((double) i / 1E16) + " : " + five.evaluateKernel((double)i/1E16));
        }
        System.out.println();
        for (int i = 0; i < 100; i++) {
            System.out.println(((double) i / 1E26) + " : " + five.evaluateKernel((double)i/1E26));
        }

        System.exit(0);*/


        //Normal distribution
        Double[] samples = new Double[10000];
        NormalDistribution dist = new NormalDistribution(0.5, 0.10);
        for (int i = 0; i < samples.length; i++) {
            samples[i] = (Double)dist.nextRandom();
            if (samples[i] < 0.0 && DEBUG) {
                System.err.println("Negative value generated!");
            }
        }
        Arrays.sort(samples);
        if (DEBUG) {
            System.out.println("min: " + samples[0]);
            System.out.println("max: " + samples[samples.length - 1] + "\n");
        }

        LogTransformedNormalKDEDistribution ltn = new LogTransformedNormalKDEDistribution(samples);
        NormalKDEDistribution nKDE = new NormalKDEDistribution(samples);

        if (DEBUG) {
            for (int i = 0; i < 50; i++) {
                Double test = (Double) dist.nextRandom();
                System.out.println("random draw: " + test);
                System.out.println("normal KDE: " + nKDE.evaluateKernel(test));
                System.out.println("log transformed normal KDE: " + ltn.evaluateKernel(test) + "\n");
            }
        }


        //LogNormal distribution
        samples = new Double[2000];
        dist = new NormalDistribution(0.0, 1.0);
        for (int i = 0; i < samples.length; i++) {
            samples[i] = (Double)dist.nextRandom();
            while (samples[i] < 0.0) {
                samples[i] = (Double)dist.nextRandom();
            }
            samples[i] = Math.exp(samples[i]-Math.exp(1.0));
        }

        //Generate R code for visualisation
        System.out.print("par(mfrow=c(2,2))\n\nsamples <- c(");
        for (int i = 0; i < samples.length-1; i++) {
            System.out.print(samples[i] + ",");
        }
        System.out.println(samples[samples.length-1] + ")\n");
        System.out.println("hist(samples, 200)\nminimum=min(samples)\nabline(v=minimum,col=2,lty=2)\n");
        System.out.println("plot(density(samples))\nabline(v=minimum,col=2,lty=2)\n");

        Arrays.sort(samples);
        if (DEBUG) {
            System.out.println("min: " + samples[0]);
            System.out.println("max: " + samples[samples.length - 1] + "\n");
        }

        ltn = new LogTransformedNormalKDEDistribution(samples);
        nKDE = new NormalKDEDistribution(samples);

        System.out.print("normalKDE <- c(");
        for (int i = 0; i < 1999; i++) {
            Double test = 0.0 + ((double)i)/((double)1000);
            System.out.print(nKDE.evaluateKernel(test) + ",");
        }
        System.out.println(nKDE.evaluateKernel(((double)1999)/((double)1000)) + ")\n");
        System.out.println("index <- seq(0.0,1.999,by=0.001)");
        System.out.println("plot(index,normalKDE,type=\"l\")\nabline(v=minimum,col=2,lty=2)\n");

        System.out.print("TransKDE <- c(");
        for (int i = 0; i < 1999; i++) {
            Double test = 0.0 + ((double)i)/((double)1000);
            System.out.print(ltn.evaluateKernel(test) + ",");
        }
        System.out.println(ltn.evaluateKernel(((double)1999)/((double)1000)) + ")\n");
        System.out.println("plot(index,TransKDE,type=\"l\")\nabline(v=minimum,col=2,lty=2)");

    }

}