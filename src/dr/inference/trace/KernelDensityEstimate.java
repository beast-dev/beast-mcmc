package dr.inference.trace;

import dr.math.distributions.GammaKDEDistribution;
import dr.math.distributions.KernelDensityEstimatorDistribution;
import dr.math.distributions.NormalKDEDistribution;
import dr.stats.Variate;
import dr.util.FrequencyDistribution;

import java.awt.*;
import java.awt.geom.GeneralPath;

/**
 * @author Marc A. Suchard
 * @author Andrew Rambaut
 * @version $Id$
 */
public class KernelDensityEstimate extends DensityEstimate {

    public KernelDensityEstimate(double[] samples, int minimumBinCount) {
        this(KernelDensityEstimatorDistribution.Type.GAUSSIAN, samples, minimumBinCount);
    }

    public KernelDensityEstimate(KernelDensityEstimatorDistribution.Type type, double[] samples, int minimumBinCount) {
        super(samples, minimumBinCount);
        
        switch (type) {
            case GAUSSIAN: this.kde = new NormalKDEDistribution(samples);
            case GAMMA: this.kde = new GammaKDEDistribution(samples);
            default:
                throw new RuntimeException("Unknown type");
        }
    }
        
    protected void calculateDensity(Variate data, int minimumBinCount) {

        FrequencyDistribution frequency = calculateFrequencies(data, minimumBinCount);

        xCoordinates = new Variate.Double();
        yCoordinates = new Variate.Double();

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

}
