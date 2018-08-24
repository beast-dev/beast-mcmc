/*
 * BetaKDEDistribution.java
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

import dr.math.GammaFunction;
import dr.stats.DiscreteStatistics;

//import java.io.BufferedReader;
//import java.io.FileReader;
//import java.io.FileNotFoundException;
//import java.io.IOException;
//import java.util.StringTokenizer;

/**
 * @author Marc Suchard
 */
public class BetaKDEDistribution extends KernelDensityEstimatorDistribution {

     public BetaKDEDistribution(Double[] sample, Double lowerBound, Double upperBound, Double bandWidth) {
         super(sample, lowerBound, upperBound, bandWidth);

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
         if (lowerBound == null || upperBound == null || upperBound - lowerBound <= 0) {
             throw new RuntimeException("BetaKDEDistribution must be bounded");
         }
         if (lowerBound > DiscreteStatistics.min(sample) || upperBound < DiscreteStatistics.max(sample))
             throw new RuntimeException("Sample range outside bounds: "+DiscreteStatistics.min(sample)+" -> "+DiscreteStatistics.max(sample));

         this.lowerBound = lowerBound;
         this.upperBound = upperBound;

         // Make a copy because we are translating values into [0,1)
         double[] oldSample = sample;
         sample = new double[sample.length];

         range = upperBound - lowerBound;
         for(int i=0; i<N; i++)
             sample[i] = (oldSample[i] - this.lowerBound) / range;
     }

     protected void setBandWidth(Double bandWidth) {
         if (bandWidth == null) {
         // Default bandwidth
            double sigma = DiscreteStatistics.stdev(sample);
            this.bandWidth = sigma * Math.pow(N,-0.4); // Renault & Scaillet (2004); Chen (1999)
         } else
             this.bandWidth = bandWidth;
     }

     /* Beta kernel density estimator based on:
      * Chen SX (1999) Beta kernel estimators for density functions. Comp. Statist. Data. Anal. 31, 131-145.
      *
      */
     protected double evaluateKernel(double x) {

         double xPrime = (x - lowerBound) / range;
         double alphaMinus1 = xPrime/bandWidth - 1.0;
         double betaMinus1 = (1.0 - xPrime)/bandWidth - 1.0;

         if (xPrime < 2*bandWidth)           // Removing these two cases reduces the kernel to C1
             alphaMinus1 = getRho(xPrime,bandWidth) - 1.0;
         else if (xPrime > 1 - 2*bandWidth)
             betaMinus1 = getRho(1.0-xPrime,bandWidth) - 1.0;

         double logK = GammaFunction.lnGamma(alphaMinus1+betaMinus1+2.0) - GammaFunction.lnGamma(alphaMinus1+1.0) - GammaFunction.lnGamma(betaMinus1+1.0);

         double pdf = 0;
         for(int i=0; i<N; i++)
             pdf +=  Math.pow(sample[i],alphaMinus1) * Math.pow(1.0-sample[i],betaMinus1);

         return pdf * Math.exp(logK) / (double)N / range;
     }

     private double getRho(double x, double bandWidth) {
         return 2*bandWidth*bandWidth + 2.5 - Math.sqrt(4*bandWidth*bandWidth*bandWidth*bandWidth
                                                            + 6*bandWidth*bandWidth
                                                            + 2.25
                                                            - x*x
                                                            - x/bandWidth);
     }

    public double sampleMean() { return DiscreteStatistics.mean(sample); }

    private double range;

// public static void main(String[] args) {
//
//     String fileName = "out.txt";
//     double[] values = null;
//            try {
//
//         BufferedReader reader = new BufferedReader(new FileReader(fileName));
//
//         String line1 = reader.readLine();
//         StringTokenizer st = new StringTokenizer(line1," ");
//         values = new double[st.countTokens()];
//         for(int i=0; i<values.length; i++)
//             values[i] = Double.valueOf(st.nextToken());
//
//         reader.close();
//
//     } catch (FileNotFoundException e) {
//         System.err.println("File not found: "+fileName);
//         System.exit(-1);
//     } catch (IOException e) {
//         System.err.println("IO exception reading: "+fileName);
//         System.exit(-1);
//     }
//
////        System.err.println("v: "+new Vector(values));
//     BetaKDEDistribution kde;
//     kde = new BetaKDEDistribution(values,0.0, 10.0, 0.1);
//     System.err.println("r: "+kde.pdf(3.0));
//
//     kde = new BetaKDEDistribution(values,1.0, 9.0, 0.25);
//     System.err.println("r: "+kde.pdf(7.0));
//
//     System.err.println("sm: "+kde.sampleMean());
//
// }
}
