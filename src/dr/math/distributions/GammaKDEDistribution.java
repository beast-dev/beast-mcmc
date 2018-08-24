/*
 * GammaKDEDistribution.java
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

import dr.stats.DiscreteStatistics;

import java.util.Random;


/**
 * @author Jennifer Tom
 *         Based on S. X. Chen. Probability density function estimation using gamma kernels.
 *         Annals of the Institute of Statistical Mathematics, 52(3):471-480, 2000.
 *         Use to create KDE for positive valued functions
 *         Assumes limits are (0, inf)
 *         Must provide with a bandwidth, or defaults to Scott's Rule
 *         Univariate distribution only
 */
public class GammaKDEDistribution extends KernelDensityEstimatorDistribution {


    public GammaKDEDistribution(Double[] sample) {
        this(sample, null);
    }

    public GammaKDEDistribution(Double[] sample, Double bandWidth) {
        super(sample, 0.0, Double.POSITIVE_INFINITY, bandWidth);

    }

    @Override
    public double getFromPoint() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double getToPoint() {
        throw new RuntimeException("Not yet implemented");
    }

    protected void processBounds(Double lowerBound, Double upperBound) {
        if (lowerBound > DiscreteStatistics.min(sample)) {
            throw new RuntimeException("Sample min out of bounds.  Gamma kernel for use with positive data only: " + DiscreteStatistics.min(sample));
        } else if (upperBound < DiscreteStatistics.max(sample)) {
            throw new RuntimeException("Sample max out of bounds" + DiscreteStatistics.max(sample));
        }
        this.lowerBound = lowerBound;
        this.upperBound = upperBound;

    }

    protected void setBandWidth(Double bandWidth) {
        if (bandWidth == null) {
            double sigma = DiscreteStatistics.stdev(sample);
            //Scott's rule  (Hardle, 2004, Nonparametric and Semiparameteric Models)
            this.bandWidth = sigma * Math.pow(N, -0.2);
        } else
            this.bandWidth = bandWidth;

    }

    protected double evaluateKernel(double x) {

        double shape;
        double scale;

        if (x >= 2 * bandWidth) {
            shape = x / bandWidth;
        } else {
            shape = .25 * Math.pow(x / bandWidth, 2) + 1;
        }
        scale = bandWidth;
        double pdf = 0;
        for (int i = 0; i < N; i++) {
            pdf +=
                    Math.pow(sample[i], shape - 1) * Math.exp(-sample[i] / scale) / (Math.pow(scale, shape) * gamma(shape));
        }
        return pdf / N;
    }

    private double gamma(double value) {
        return cern.jet.stat.Gamma.gamma(value);

    }

    public static void main(String[] args) {

        long start = System.currentTimeMillis();

        Random random = new Random(1234);

        Double[] samples = new Double[10000000];
        for (int i = 0; i < samples.length; i++) {
            samples[i] = random.nextDouble();
        }
        GammaKDEDistribution nKDE = new GammaKDEDistribution(samples);

        for (int i = 0; i < 100; i++) {
            nKDE.evaluateKernel(random.nextDouble());
        }

        long end = System.currentTimeMillis();

        System.out.println("Time: " + (end-start));

    }

//    private double sampleMean() {return DiscreteStatistics.mean(sample);}


// public static void main(String[] args) {
//    String fileName = "/Users/jen/School/Programs/BEAST/kdeTest/simulUExp.txt";
//    //String fileName = "/Users/jen/School/Programs/BEAST/kdeTest/simulUNorm.txt";
//    double[] values = null;  //vector of the training set
//    try{
//    BufferedReader br = new BufferedReader(new FileReader(fileName));
//    StringBuilder AllDataSb = new StringBuilder();
//    String s;
//    String myLineSeparator = "\r";
//   while (( s = br.readLine()) != null){
//    AllDataSb.append(s);
//        AllDataSb.append(myLineSeparator);
//    }
//    //convert from stringbuilder to string
//    String AllDataS = AllDataSb.toString();
//    String[] result = AllDataS.split("\t|\r|,");
//
//    //convert string to double
//     values = new double[result.length];
//     for (int i = 0; i < result.length; i++) {values[i] = Double.parseDouble(result[i]);}
//
//    } //close try
//    catch(Exception e){
//        System.out.println(e.getMessage());
//        System.exit(-1);
//    }
//
//
//     GammaKDEDistribution kde;
//
//
//     //expecting 2.02177 -> .073; .4046729 -> 1.0257; 0.1502078 -> 1.4021
//     //normal-reference bandwidth              0.1939
//     kde = new GammaKDEDistribution(values, null);
//     System.err.println("prediction: at 2.02: "+kde.pdf(2.02177)+" at 0.405: "+kde.pdf(0.4046729)+" at 0.15: "+kde.pdf(0.1502078));
//     System.err.println("sm: "+kde.sampleMean());
// }
}
