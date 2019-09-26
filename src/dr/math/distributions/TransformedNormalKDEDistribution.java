/*
 * TransformedNormalKDEDistribution.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.util.Transform;

import java.util.Arrays;
import java.util.Random;


/**
 * @author Guy Baele
 * @author Paul Bastide
 */
public class TransformedNormalKDEDistribution extends NormalKDEDistribution {

    public static final int MINIMUM_GRID_SIZE = 2048;
    public static final boolean DEBUG = false;

    private final Transform transform;

    //the samples should not already be transformed (the transformation is done in this class)
    public TransformedNormalKDEDistribution(Double[] sample, Transform transform) {
        this(sample, transform, null, null, null);
    }

    private TransformedNormalKDEDistribution(Double[] sample, Transform transform, Double lowerBound, Double upperBound, Double bandWidth) {
        this(sample, transform, lowerBound, upperBound, bandWidth, 3.0, MINIMUM_GRID_SIZE);
    }

    private TransformedNormalKDEDistribution(Double[] sample, Transform transform, Double lowerBound, Double upperBound, Double bandWidth, double cut, int n) {

        super(getTransform(sample, transform), lowerBound, upperBound, bandWidth, cut, n);
        if (DEBUG) {
            System.out.println("Creating the KDE in transformed space");
            System.out.println("lowerBound = " + lowerBound);
            System.out.println("upperBound = " + upperBound);
        }
        this.transform = transform;
    }

    private static Double[] getTransform(Double[] x, Transform transform) {
        Double[] logX = new Double[x.length];
        for (int i = 0; i < logX.length; i++) {
            logX[i] = transform.transform(x[i]);
        }
        return logX;
    }

    public double getFromPoint() {
        return from;
    }

    public double getToPoint() {
        return to;
    }

    @Override
    void makeXGrid() {
        super.makeXGrid();
        for (int i = 0; i < xPoints.length; i++) {
            xPoints[i] = transform.inverse(xPoints[i]);
        }
    }

    @Override
    void transformData() {
        super.transformData();
        for (int i = 0; i < densityPoints.length; i++) {
            densityPoints[i] = densityPoints[i] / transform.gradient(xPoints[i]);
        }
    }


    @Override
    protected void processBounds(Double lowerBound, Double upperBound) {
        if ((lowerBound != null && lowerBound != Double.NEGATIVE_INFINITY) ||
                (upperBound != null && upperBound != Double.POSITIVE_INFINITY)) {
            throw new RuntimeException("TransformedNormalKDEDistribution must be unbounded");
        }
    }

    public static TransformedNormalKDEDistribution getLogTransformedNormalKDEDistribution(Double[] sample, Double lowerBound, Double upperBound, Double bandWidth, double cut, int n) {
        return new TransformedNormalKDEDistribution(sample, new Transform.LogTransform(), lowerBound, upperBound, bandWidth, cut, n);
    }

    public static TransformedNormalKDEDistribution getLogTransformedNormalKDEDistribution(Double[] sample) {
        return new TransformedNormalKDEDistribution(sample, new Transform.LogTransform());
    }

    public static TransformedNormalKDEDistribution getLogitTransformedNormalKDEDistribution(Double[] sample, Double upperLimit, Double lowerBound, Double upperBound, Double bandWidth, double cut, int n) {
        return new TransformedNormalKDEDistribution(sample, new Transform.ScaledLogitTransform(upperLimit, 0.0), lowerBound, upperBound, bandWidth, cut, n);
    }

    public static TransformedNormalKDEDistribution getLogitTransformedNormalKDEDistribution(Double[] sample, Double upperLimit) {
        return new TransformedNormalKDEDistribution(sample, new Transform.ScaledLogitTransform(upperLimit, 0.0));
    }


    public static void main(String[] args) {

//        //Normal distribution
//        Double[] samples = new Double[10000];
//        NormalDistribution dist = new NormalDistribution(0.5, 0.10);
//        for (int i = 0; i < samples.length; i++) {
//            samples[i] = (Double) dist.nextRandom();
//            if (samples[i] < 0.0 && DEBUG) {
//                System.err.println("Negative value generated!");
//            }
//        }
//        Arrays.sort(samples);
//        if (DEBUG) {
//            System.out.println("min: " + samples[0]);
//            System.out.println("max: " + samples[samples.length - 1] + "\n");
//        }
//
//        LogTransformedNormalKDEDistribution ltn = new LogTransformedNormalKDEDistribution(samples);
//        NormalKDEDistribution nKDE = new NormalKDEDistribution(samples);

//        if (DEBUG) {
//            for (int i = 0; i < 50; i++) {
//                Double test = (Double) dist.nextRandom();
//                System.out.println("random draw: " + test);
//                System.out.println("normal KDE: " + nKDE.evaluateKernel(test));
//                System.out.println("log transformed normal KDE: " + ltn.evaluateKernel(test) + "\n");
//            }
//        }


//        //LogNormal distribution
//        int size = 100;
//        samples = new Double[2 * size];
//        dist = new NormalDistribution(0.0, 1.0);
//        for (int i = 0; i < samples.length; i++) {
//            samples[i] = (Double)dist.nextRandom();
//            while (samples[i] < 0.0) {
//                samples[i] = (Double)dist.nextRandom();
//            }
//            samples[i] = Math.exp(samples[i]-Math.exp(1.0));
//        }

        //Gamma distribution
        Random randUnif = new Random();
        randUnif.setSeed(17920920); // The day before
        int size = 1000;
        Double[] samples = new Double[2 * size];
        double scale = 0.001;
        for (int i = 0; i < samples.length; i++) {
            samples[i] = randUnif.nextDouble();
            samples[i] = -Math.log(samples[i]) * scale; // Gamma(1, scale)
        }

        //Generate R code for visualisation
        System.out.print("par(mfrow=c(2,3))\n\nsamples <- c(");
        for (int i = 0; i < samples.length - 1; i++) {
            System.out.print(samples[i] + ",");
            if (i % 10 == 0) System.out.print("\n");
        }
        System.out.println(samples[samples.length - 1] + ")\n");
        System.out.println("ii <- 0:" + 2 * size + "/" + size);
        System.out.println("hist(samples," + 2 * size + ",freq=FALSE)\nminimum=min(samples)\nabline(v=minimum,col=2,lty=2)");
        System.out.println("lines(ii, dgamma(ii, shape=1,scale=" + scale + "),col=\"green\")\n");
        System.out.println("plot(density(samples))\nabline(v=minimum,col=2,lty=2)");
        System.out.println("lines(ii, dgamma(ii, shape=1,scale=" + scale + "),col=\"green\")\n");

        Arrays.sort(samples);
        if (DEBUG) {
            System.out.println("min: " + samples[0]);
            System.out.println("max: " + samples[samples.length - 1] + "\n");
        }

        LogTransformedNormalKDEDistribution ltnOld = new LogTransformedNormalKDEDistribution(samples);
        TransformedNormalKDEDistribution ltnNew = getLogTransformedNormalKDEDistribution(samples);
        NormalKDEDistribution nKDE = new NormalKDEDistribution(samples);

        System.out.print("index <- c(");
        double[] testPoints = new double[2 * size];
        double low = 0.0;
        double up = 1.5 * java.util.Collections.max(Arrays.asList(samples));
        double x = low;
        double delta = (up - low) / (2 * size - 1);
        for (int i = 0; i < 2 * size - 1; i++) {
            System.out.print(x + ",");
            if (i % 10 == 0) System.out.print("\n");
            testPoints[i] = x;
            x += delta;
        }
        testPoints[2 * size - 1] = x;
        System.out.println(x + ")\n");

        System.out.print("normalKDE <- c(");
        for (int i = 0; i < 2 * size - 1; i++) {
            Double test = testPoints[i];
            System.out.print(nKDE.evaluateKernel(test) + ",");
            if (i % 10 == 0) System.out.print("\n");
        }
        System.out.println(nKDE.evaluateKernel(((double) 2 * size - 1) / ((double) size)) + ")\n");

        System.out.print("TransKDEOld <- c(");
        for (int i = 0; i < 2 * size - 1; i++) {
            Double test = testPoints[i];
            System.out.print(ltnOld.evaluateKernel(test) + ",");
            if (i % 10 == 0) System.out.print("\n");
        }
        System.out.println(ltnOld.evaluateKernel(((double) 2 * size - 1) / ((double) size)) + ")\n");

        System.out.print("TransKDENew <- c(");
        for (int i = 0; i < 2 * size - 1; i++) {
            Double test = testPoints[i];
            System.out.print(ltnNew.evaluateKernel(test) + ",");
            if (i % 10 == 0) System.out.print("\n");
        }
        System.out.println(ltnNew.evaluateKernel(((double) 2 * size - 1) / ((double) size)) + ")\n");


        System.out.println("plot(index,normalKDE,type=\"l\",xlab=\"Bandwidth=" + nKDE.bandWidth + "\")\nabline(v=minimum,col=2,lty=2)");
        System.out.println("lines(ii, dgamma(ii, shape=1,scale=" + scale + "),col=\"green\")");

        System.out.println("plot(index,TransKDEOld,type=\"l\",xlab=\"Bandwidth=" + ltnOld.bandWidth + "\")\nabline(v=minimum,col=2,lty=2)");
        System.out.println("lines(ii, dgamma(ii, shape=1,scale=" + scale + "),col=\"green\")");

        System.out.println("plot(index,TransKDENew,type=\"l\",xlab=\"Bandwidth=" + ltnNew.bandWidth + "\")\nabline(v=minimum,col=2,lty=2)");
        System.out.println("lines(ii, dgamma(ii, shape=1,scale=" + scale + "),col=\"green\")");

        // Timing
        long start = System.currentTimeMillis();

        Random random = new Random(1234);

        Double[] samplesTiming = new Double[10000000];
        for (int i = 0; i < samplesTiming.length; i++) {
            samplesTiming[i] = random.nextDouble();
        }
        TransformedNormalKDEDistribution nKDETiming = getLogTransformedNormalKDEDistribution(samplesTiming);

        for (int i = 0; i < 100; i++) {
            nKDETiming.evaluateKernel(random.nextDouble());
        }

        long end = System.currentTimeMillis();

        System.out.println("Time: " + (end - start));

    }

}
