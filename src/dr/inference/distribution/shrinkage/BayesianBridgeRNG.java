package dr.inference.distribution.shrinkage;

import dr.math.MathUtils;
import dr.math.distributions.GammaDistribution;

public class BayesianBridgeRNG {

    private static double drawUnregularizedCoefficient(double globalScale, double exponent) {
        // This is how rgnorm works in the R package gnorm
        double lambda = Math.pow((1.0 / globalScale),exponent);
        double u = MathUtils.nextDouble();
        double unsigned = Math.pow(GammaDistribution.quantile(u, 1.0 / exponent, 1.0 / lambda),1.0 / exponent);
        return MathUtils.nextBoolean() ? unsigned : -unsigned;
    }

    private static double drawRegularizedCoefficientBridgeProposal(double globalScale, double exponent, double slabWidth) {
        double x = Double.NaN;
        boolean done = false;
        double twoSlabSquared = 2.0 * slabWidth * slabWidth;
        while (!done) {
            double prop = drawUnregularizedCoefficient(globalScale, exponent);
            double logAcceptProb = -(prop * prop)/twoSlabSquared;
            if ( Math.log(MathUtils.nextDouble()) <= logAcceptProb) {
                x = prop;
                done = true;
            }
        }
        return x;
    }

    public static double drawRegularizedCoefficientNormalProposal(double globalScale, double exponent, double slabWidth) {
        double x = Double.NaN;
        boolean done = false;
        double logSlab2SqrtPi = Math.log(slabWidth * Math.sqrt(2.0 * Math.PI));
        while (!done) {
            double prop = MathUtils.nextGaussian() * slabWidth;
            double logAcceptProb = -Math.pow(Math.abs(prop / globalScale),exponent) - logSlab2SqrtPi;;
            if ( Math.log(MathUtils.nextDouble()) <= logAcceptProb) {
                x = prop;
                done = true;
            }
        }
        return x;
    }

    public static double drawRegularizedCoefficient(double globalScale, double exponent, double slabWidth) {
        double x;
        if ( globalScale > slabWidth ) {
            x = drawRegularizedCoefficientNormalProposal(globalScale, exponent, slabWidth);
        } else {
            x = drawRegularizedCoefficientBridgeProposal(globalScale, exponent, slabWidth);
        }
        return x;
    }

    // For Bridge without slab
    public static double[] nextRandom(double globalScale, double exponent, int dim) {
        double[] draws = new double[dim];
        for (int i = 0; i < dim; i++) {
            draws[i] = drawUnregularizedCoefficient(globalScale, exponent);
        }
        return draws;
    }

    // For Bridge with slab
    public static double[] nextRandom(double globalScale, double exponent, double slabWidth, int dim) {
        double[] draws = new double[dim];
        for (int i = 0; i < dim; i++) {
            draws[i] = drawRegularizedCoefficient(globalScale, exponent, slabWidth);
        }
        return draws;
    }
}
