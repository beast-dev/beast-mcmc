/*
 * LogitTransformedNormalKDEDistribution.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

/**
 * @author Guy Baele
 */
public class LogitTransformedNormalKDEDistribution extends KernelDensityEstimatorDistribution {

    public static final int MINIMUM_GRID_SIZE = 2048;
    public static final boolean DEBUG = false;

    //the samples should not already be logit transformed (the logit transformation is done in this class)
    public LogitTransformedNormalKDEDistribution(Double[] sample) {
        this(sample, 1.0, null, null, null);
    }

    public LogitTransformedNormalKDEDistribution(Double[] sample, Double upperLimit) {
        this(sample, upperLimit, null, null, null);
    }

    public LogitTransformedNormalKDEDistribution(Double[] sample, int n) {
        this(sample, 1.0, null, null, null, 3.0, n);
    }

    public LogitTransformedNormalKDEDistribution(Double[] sample, Double upperLimit, int n) {
        this(sample, upperLimit, null, null, null, 3.0, n);
    }

    public LogitTransformedNormalKDEDistribution(Double[] sample, Double upperLimit, Double lowerBound, Double upperBound, Double bandWidth) {
        this(sample, upperLimit, lowerBound, upperBound, bandWidth, 3.0, MINIMUM_GRID_SIZE);
    }

    public LogitTransformedNormalKDEDistribution(Double[] sample, Double upperLimit, Double lowerBound, Double upperBound, Double bandWidth,
                                                 int n) {
        this(sample, upperLimit, lowerBound, upperBound, bandWidth, 3.0, n);
    }

    public LogitTransformedNormalKDEDistribution(Double[] sample, Double upperLimit, Double lowerBound, Double upperBound, Double bandWidth, double cut, int n) {

        super(sample, lowerBound, upperBound, bandWidth);
        //transform the data to the logit scale and store in logSample
        this.upperLimit = upperLimit;
        if (DEBUG) {
            System.out.println("Creating the KDE in logit space");
            System.out.println("lowerBound = " + lowerBound);
            System.out.println("upperBound = " + upperBound);
            System.out.println("upperlimit = " + upperLimit);
        }

        this.logitSample = new double[sample.length];
        for (int i = 0; i < logitSample.length; i++) {
            //this.logitSample[i] = Math.log(sample[i] / (1.0 - sample[i]));
            this.logitSample[i] = Math.log(sample[i] / (upperLimit - sample[i]));
        }
        //keep a backup copy of the samples in normal space
        this.backupSample = new double[sample.length];
        for (int i = 0; i < sample.length; i++) {
            this.backupSample[i] = sample[i];
        }
        //overwrite the stored samples, sample.length stays the same
        this.sample = logitSample;
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
        //we're still in logit space and need to return to normal space
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
            System.out.println("upperlimit = " + upperLimit);
        }

        this.sample = backupSample;
        //processBounds(lowerBound, upperBound);
        setBandWidth(null);

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

        /*if (from < 0.0) {
        	from = 0.0;
        }
        if (lo < 0.0) {
        	lo = 0.0;
        }*/

        //make new ordinates for the transformation back to normal space
        //need a backup of the xPoints for the logit scale KDE
        this.backupXPoints = new double[xPoints.length];
        System.arraycopy(xPoints, 0, backupXPoints, 0, xPoints.length);
        makeOrdinates();

        int numberOfNegatives = 0;
        if (DEBUG) {
            System.out.println("\nxPoints length = " + xPoints.length);
        }
        for (int i = 0; i < xPoints.length; i++) {
            if (xPoints[i] < 0.0) {
                numberOfNegatives++;
            }
            //System.out.println(xPoints[i]);
        }
        if (DEBUG) {
            System.out.println("number of negative xPoints = " + numberOfNegatives);
        }

        //the KDE on logit scale is contained in the xPoints and densityPoints arrays
        //copy them to finalXPoints and finalDensityPoints
        this.finalXPoints = new double[xPoints.length - numberOfNegatives];
        System.arraycopy(xPoints, numberOfNegatives, finalXPoints, 0, xPoints.length - numberOfNegatives);
        if (DEBUG) {
            for (int i = 0; i < xPoints.length; i++) {
                System.out.println(backupXPoints[i] + " : " + densityPoints[i]);
            }
            System.out.println("\nfinalXPoints length = " + finalXPoints.length);
        }

        this.finalDensityPoints = new double[densityPoints.length - numberOfNegatives];

        for (int i = 0; i < finalXPoints.length; i++) {
            finalDensityPoints[i] = linearApproximate(backupXPoints, densityPoints, Math.log(finalXPoints[i]/(upperLimit - finalXPoints[i])), 0.0, 0.0)*(1.0/(finalXPoints[i]*(upperLimit-finalXPoints[i])));
            if (DEBUG) {
                System.out.println(finalXPoints[i] + "\t" + finalDensityPoints[i]);
            }
        }

        //System.exit(0);
    	
    }

    private void transformData() {
        ComplexArray Y  = new ComplexArray(massdist(this.logitSample, lo, up, this.gridSize));
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
        //System.out.println("X");
        for (int i = 0; i < gridSize; i++) {
            xPoints[i] = x;
            x += delta;
            //System.out.println(xPoints[i]);
        }
        //System.out.println();
    }

    @Override
    protected double evaluateKernel(double x) {        
        if (!densityKnown) {
        	//computeDensity() calls makeOrdinates and transformData
        	computeDensity();
        }
        //xPoints and densityPoints are now back in normal space
        return linearApproximate(finalXPoints, finalDensityPoints, x, 0.0, 0.0);
    }

    @Override
    protected void processBounds(Double lowerBound, Double upperBound) {
        if ((lowerBound != null && lowerBound != Double.NEGATIVE_INFINITY) ||
                (upperBound != null && upperBound != Double.POSITIVE_INFINITY)) {
            throw new RuntimeException("LogitTransformedNormalKDEDistribution must be unbounded");
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
    private double[] xPoints, finalXPoints, backupXPoints;
    private double[] densityPoints, finalDensityPoints;
    private double[] backupSample, logitSample;

    private int gridSize;
    private double cut;
    private double from;
    private double to;
    private double lo;
    private double up;

    private double upperLimit;

    private boolean densityKnown = false;

    public static void main(String[] args) {

        double upperlimit = 3.0;

        //Generate a normal distribution, truncated at 0.0
        Double[] samples = new Double[2001];
        NormalDistribution dist = new NormalDistribution(0.5, 0.08);
        for (int i = 0; i < samples.length; i++) {
            samples[i] = upperlimit-Math.abs((Double)dist.nextRandom()-0.5);
        }

        //Generate R code for visualisation in [0,1]
        System.out.print("par(mfrow=c(2,2))\n\nsamples <- c(");
        for (int i = 0; i < samples.length-1; i++) {
            System.out.print(samples[i] + ",");
        }
        System.out.println(samples[samples.length-1] + ")\n");
        System.out.println("hist(samples,200,xlim=c(" + (upperlimit-0.5) + "," + (upperlimit+0.5) + "))\n");
        System.out.println("plot(density(samples),xlim=c(" + (upperlimit-0.5) + "," + (upperlimit+0.5) + "))\n");

        LogitTransformedNormalKDEDistribution ltn = new LogitTransformedNormalKDEDistribution(samples, upperlimit);
        NormalKDEDistribution nKDE = new NormalKDEDistribution(samples);

        System.out.print("normalKDE <- c(");
        for (int i = ((int)upperlimit-1)*1000; i < ((int)upperlimit+1)*1000; i++) {
            Double test = 0.0 + ((double)i)/(1000);
            System.out.print(nKDE.evaluateKernel(test) + ",");
        }
        System.out.println(nKDE.evaluateKernel((((int)upperlimit+1)*1000)/((upperlimit*1000))) + ")\n");
        System.out.println("index <- seq(" + (upperlimit-1) + "," + (upperlimit+1) + ",by=0.001)");
        System.out.println("plot(index,normalKDE,type=\"l\",xlim=c(" + (upperlimit-0.5) + "," + (upperlimit+0.5) + "))\n");

        System.out.print("TransKDE <- c(");
        for (int i = ((int)upperlimit-1)*1000; i < ((int)upperlimit+1)*1000; i++) {
            Double test = 0.0 + ((double)i)/(1000);
            System.out.print(ltn.evaluateKernel(test) + ",");
        }
        System.out.println(ltn.evaluateKernel((((int)upperlimit+1)*1000)/((upperlimit*1000))) + ")\n");
        System.out.println("plot(index,TransKDE,type=\"l\",xlim=c(" + (upperlimit-0.5) + "," + (upperlimit+0.5) + "))\n");

    }

}