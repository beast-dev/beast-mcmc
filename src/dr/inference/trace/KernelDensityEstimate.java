/*
 * KernelDensityEstimate.java
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

package dr.inference.trace;

import dr.math.distributions.GammaKDEDistribution;
import dr.math.distributions.KernelDensityEstimatorDistribution;
import dr.math.distributions.NormalKDEDistribution;
import dr.stats.Variate;
import dr.util.FrequencyDistribution;

/**
 * @author Marc A. Suchard
 * @author Andrew Rambaut
 * @version $Id$
 */
public class KernelDensityEstimate extends DensityEstimate {

    public KernelDensityEstimate(Double[] samples, int minimumBinCount) {
        this(KernelDensityEstimatorDistribution.Type.GAUSSIAN, samples, minimumBinCount);
    }

    public KernelDensityEstimate(KernelDensityEstimatorDistribution.Type type, Double[] samples, int minimumBinCount) {
        super(samples, minimumBinCount);
        
        switch (type) {
            case GAUSSIAN: this.kde = new NormalKDEDistribution(samples);
            case GAMMA: this.kde = new GammaKDEDistribution(samples);
            default:
                throw new RuntimeException("Unknown type");
        }
    }
    
    public KernelDensityEstimate(String type, Double[] samples, int minimumBinCount) {
        super(samples, minimumBinCount);
       
        this.kde = new NormalKDEDistribution(samples);

    }
    
    public KernelDensityEstimatorDistribution getKernelDensityEstimatorDistribution() {
    	return this.kde;
    }
        
    protected void calculateDensity(Variate data, int minimumBinCount) {
    	
        FrequencyDistribution frequency = calculateFrequencies(data, minimumBinCount);

        xCoordinates = new Variate.D();
        yCoordinates = new Variate.D();

        double x = frequency.getLowerBound() - (frequency.getBinSize() / 2.0);
        int extraEdgeCount = 0;
        while (kde.pdf(x) > minDensity && x > lowerBoundary) {
            x -= frequency.getBinSize();
            extraEdgeCount += 1;
        }
        xCoordinates.add(x);
        yCoordinates.add(0.0);
        x += frequency.getBinSize();
        int count = 0;
        while (count < (frequency.getBinCount() + extraEdgeCount)) {// ||
//                (kde.pdf(x) > minDensity && x < upperBoundary)) {
            xCoordinates.add(x);
            yCoordinates.add(kde.pdf(x));
            x += frequency.getBinSize();
            count++;
        }
//        System.err.println("kde = " + kde.pdf(x));
        while (kde.pdf(x) > minDensity ) {
//            System.err.println("add bit on end!!!");
            xCoordinates.add(x);
            yCoordinates.add(kde.pdf(x));
            x += frequency.getBinSize();
        }
        xCoordinates.add(x);
        yCoordinates.add(0.0);

//
//
//        int extraBinsOnEdges = 5;
//        double x = frequency.getLowerBound() - extraBinsOnEdges * frequency.getBinSize();
//        for (int i = 0; i < frequency.getBinCount() + 2 * extraBinsOnEdges; i++) {
//            double xMidPoint = x + (frequency.getBinSize() / 2.0);
//            xData.add(xMidPoint);
//            yData.add(kde.pdf(xMidPoint));
//            x += frequency.getBinSize();
//        }
    }

    private KernelDensityEstimatorDistribution kde;

    private double lowerBoundary = 0;
    private double upperBoundary = Double.POSITIVE_INFINITY;
    private static final double minDensity = 10E-6;

    public static void main(String[] args) {
    	Double[] samples = new Double[101];
    	for (int i = 0; i <= 100; i++) {
    		samples[i] = ((double)i)/100.0;
    		System.out.println(samples[i]);
    	}
    	KernelDensityEstimate kdeTest = new KernelDensityEstimate("Gaussian", samples, 512);
    	Variate.D samplesTwo = new Variate.D(samples);
    	//kdeTest.calculateDensitykde(samplesTwo, 512);
    	KernelDensityEstimatorDistribution kde = kdeTest.getKernelDensityEstimatorDistribution();
    	System.out.println("Bandwidth: " + kde.getBandWidth());
    	System.out.println();
    }
    
}
