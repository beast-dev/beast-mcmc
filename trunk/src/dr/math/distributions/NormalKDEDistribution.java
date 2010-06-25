package dr.math.distributions;

import dr.stats.DiscreteStatistics;
import dr.util.HeapSort;

/**
 * @author Marc A. Suchard
 */
public class NormalKDEDistribution extends KernelDensityEstimatorDistribution {

    public NormalKDEDistribution(double[] sample, Double lowerBound, Double upperBound, Double bandWidth) {
        super(sample, lowerBound, upperBound, bandWidth);
    }

    @Override
    protected double evaluateKernel(double x) {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected void processBounds(Double lowerBound, Double upperBound) {
        if ((lowerBound != null && lowerBound != Double.NEGATIVE_INFINITY) ||
                (upperBound != null && upperBound != Double.POSITIVE_INFINITY)) {
            throw new RuntimeException("NormalKDEDistribution must be unbounded");
        }
    }

    @Override
    protected void setBandWidth(Double bandWidth) {
        if (bandWidth == null) {
            // Default bandwidth
            this.bandWidth = bandwidthNRD(sample);
        } else
            this.bandWidth = bandWidth;
    }

//   bandwidth.nrd =
//   function (x)
//   {
//       r <- quantile(x, c(0.25, 0.75))
//       h <- (r[2] - r[1])/1.34
//       4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-1/5)
//   }

    public double bandwidthNRD(double[] x) {

        int[] indices = new int[x.length];
        HeapSort.sort(x, indices);

        final double h =
                (DiscreteStatistics.quantile(0.75, x, indices) - DiscreteStatistics.quantile(0.25, x, indices)) / 1.34;
        return 1.06 *
                Math.min(Math.sqrt(DiscreteStatistics.variance(x)), h) *
                Math.pow(x.length, -0.2);
    }
}
