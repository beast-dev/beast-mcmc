package dr.app.tracer.traces;

import dr.math.distributions.KernelDensityEstimatorDistribution;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class TraceOptions {

    public enum Type {
        REAL,
        INTEGER,
        CATEGORICAL
    }

    private int histogramBinCount = 50;
    private double lowerBound = Double.NEGATIVE_INFINITY;
    private double upperBound = Double.POSITIVE_INFINITY;

    private KernelDensityEstimatorDistribution.Type kdeType = KernelDensityEstimatorDistribution.Type.GAUSSIAN;
}
