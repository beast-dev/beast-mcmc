/*
 * KDENumericalDensityPlot.java
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

package dr.app.gui.chart;

import dr.inference.trace.TraceDistribution;
import dr.math.distributions.GammaKDEDistribution;
import dr.math.distributions.KernelDensityEstimatorDistribution;
import dr.math.distributions.LogTransformedNormalKDEDistribution;
import dr.math.distributions.NormalKDEDistribution;
import dr.stats.Variate;
import dr.util.FrequencyDistribution;

import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class KDENumericalDensityPlot extends NumericalDensityPlot {
    private final static boolean DEBUG = false;

    protected static final int DEFAULT_KDE_BINS = 5000;

    public KDENumericalDensityPlot(List<Double> data) {
        super(data, DEFAULT_KDE_BINS);
    }

    public KDENumericalDensityPlot(List<Double> data, int minimumBinCount) {
        super(data, minimumBinCount);
    }

    private KernelDensityEstimatorDistribution getKDE(Double[] samples) {
//        System.err.println("samples is null? " + (samples == null ? "yes" : "no"));
//        System.err.println("type is null? " + (type == null ? "yes" : "no"));
        type = KernelDensityEstimatorDistribution.Type.GAUSSIAN;
        switch (type) {
            case GAUSSIAN: return new NormalKDEDistribution(samples);
            case GAMMA: return new GammaKDEDistribution(samples);
            case LOG_TRANSFORMED_GAUSSIAN: return new LogTransformedNormalKDEDistribution(samples);
            default:
                throw new RuntimeException("Unknown type");
        }
    }

    /**
     * Set data
     */
    public void setData(Variate.D data, int minimumBinCount) {

        setRawData(data);
        Double[] samples = new Double[data.getCount()];
        for (int i = 0; i < data.getCount(); i++) {
            samples[i] = data.get(i);
        }
        kde = getKDE(samples);

        FrequencyDistribution frequency = getFrequencyDistribution(data, minimumBinCount);

        Variate.D xData = new Variate.D();
        Variate.D yData = new Variate.D();

//        double x = frequency.getLowerBound() - frequency.getBinSize();
//        double maxDensity = 0.0;
//        // TODO Compute KDE once
////        for (int i = 0; i < frequency.getBinCount(); i++) {
////            double density = frequency.getFrequency(i) / frequency.getBinSize() / data.getValuesSize();
////            if (density > maxDensity) maxDensity = density;
////        }
//
//        xData.add(x + (frequency.getBinSize() / 2.0));
//        yData.add(0.0);
//        x += frequency.getBinSize();
//
//        for (int i = 0; i < frequency.getBinCount(); i++) {
//            double xPoint = x + (frequency.getBinSize() / 2.0);
//            xData.add(xPoint);
////            double density = frequency.getFrequency(i) / frequency.getBinSize() / data.getValuesSize();
//            double density = kde.pdf(xPoint);
//            if (relativeDensity) {
//                yData.add(density / maxDensity);
//            } else {
//                yData.add(density);
//            }
//            x += frequency.getBinSize();
//        }
//
//        xData.add(x + (frequency.getBinSize() / 2.0));
//        yData.add(0.0);
        double x = frequency.getLowerBound() - (frequency.getBinSize() / 2.0);
        int extraEdgeCount = 0;
        while (kde.pdf(x) > minDensity && x > lowerBoundary) {
            x -= frequency.getBinSize();
            extraEdgeCount += 1;
        }
        xData.add(x);
        yData.add(0.0);
        x += frequency.getBinSize();
        int count = 0;
        while (count < (frequency.getBinCount() + extraEdgeCount)) {// ||
//                (kde.pdf(x) > minDensity && x < upperBoundary)) {
            xData.add(x);
            yData.add(kde.pdf(x));
            x += frequency.getBinSize();
            count++;
        }
        if (DEBUG) {
            System.err.println("kde = " + kde.pdf(x));
        }
        while (kde.pdf(x) > minDensity ) {
            if (DEBUG) {
                System.err.println("add bit on end!!!");
            }
            xData.add(x);
            yData.add(kde.pdf(x));
            x += frequency.getBinSize();
        }
        xData.add(x);
        yData.add(0.0);



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

        setData(xData, yData);
    }

    protected Variate getXCoordinates(int numPoints) {
        Double[] points = new Double[numPoints];
        for (int i = 0; i < numPoints; i++) {
            points[i] = (double) i;
        }
        return new Variate.D(points);
    }

    protected Variate getYCoordinates(Variate.D xData) {
        final int length = xData.getCount();
        Double[] points = new Double[length];
        for (int i = 0; i < length; i++) {
            points[i] = kde.pdf(xData.get(i));
        }
        return new Variate.D(points);
    }

    protected double getQuantile(double y) {
        return kde.quantile(y);
    }

    protected double getDensity(double x) {
        return kde.pdf(x);
    }

    private KernelDensityEstimatorDistribution kde;
    private NumericalDensityPlot densityPlot;

    private KernelDensityEstimatorDistribution.Type type;

    private double lowerBoundary = 0;
    private double upperBoundary = Double.POSITIVE_INFINITY;
    private static final double minDensity = 10E-6;

//    @Override
//
//    protected void paintData(Graphics2D g2, Variate xData, Variate yData) {
//
//        double x = transformX(xData.get(0));
//        double y = transformY(yData.get(0));
//
//        GeneralPath path = new GeneralPath();
//        path.moveTo((float) x, (float) y);
//
//        int n = xData.getValuesSize();
//        boolean failed = false;
//        for (int i = 1; i < n; i++) {
//            x = transformX(xData.get(i));
//            y = transformY(yData.get(i));
//            if (x == Double.NEGATIVE_INFINITY || y == Double.NEGATIVE_INFINITY ||
//                    Double.isNaN(x) || Double.isNaN(y)) {
//                failed = true;
//            } else if (failed) {
//                failed = false;
//                path.moveTo((float) x, (float) y);
//            } else {
//                path.lineTo((float) x, (float) y);
//            }
//        }
//
//        g2.setPaint(linePaint);
//        g2.setStroke(lineStroke);
//
//        g2.draw(path);
//
//    }
}
