package dr.evomodel.clock;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.distributions.InverseGaussianDistribution;
import dr.math.distributions.LogNormalDistribution;
import dr.math.distributions.NormalDistribution;

/**
 * Calculates the likelihood of a set of rate changes in a tree, assuming a (log)normal or inverse gaussian distributed
 * change in rate at each node, with a mean of the previous (log) rate and a variance proportional to branch length.
 * cf Yang and Rannala 2006
 *
 * @author Michael Defoin Platel
 * @author Wai Lok Sibon Li
 */
public class ACLikelihood extends RateEvolutionLikelihood {
    public static final String LOGNORMAL = "logNormal";
    public static final String NORMAL = "normal";
    public static final String INVERSEGAUSSIAN = "inverseGaussian";

    public ACLikelihood(TreeModel tree, Parameter ratesParameter, Parameter variance, Parameter rootRate,
                        boolean isEpisodic, String distribution) {

        //super((isLogSpace) ? "LogNormally Distributed" : "Normally Distributed", tree, ratesParameter, rootRate, isEpisodic);
        super(distribution, tree, ratesParameter, rootRate, isEpisodic);

        //this.isLogSpace = isLogSpace;
        this.variance = variance;
        this.distribution = distribution;

        addVariable(variance);
    }

    /**
     * @return the log likelihood of the rate change from the parent to the child.
     */
    double branchRateChangeLogLikelihood(double parentRate, double childRate, double time) {
        double var = variance.getParameterValue(0);

        if (!isEpisodic())
            var *= time;

        //if (isLogSpace) {
        //    double logParentRate = Math.log(parentRate);
        //    double logChildRate = Math.log(childRate);

        //    return NormalDistribution.logPdf(logChildRate, logParentRate - (var / 2.), Math.sqrt(var)) - logChildRate;

        //} else {
        //    return NormalDistribution.logPdf(childRate, parentRate, Math.sqrt(var));
        //}
        if(distribution.equals(LOGNORMAL)) {
            return LogNormalDistribution.logPdf(childRate, Math.log(parentRate) - (var / 2.), Math.sqrt(var));
        }
        else if(distribution.equals(NORMAL)) {
            return NormalDistribution.logPdf(childRate, parentRate, Math.sqrt(var));
        }
        else if(distribution.equals(INVERSEGAUSSIAN)) { /* Inverse Gaussian */
            double shape = (parentRate * parentRate * parentRate) / var;
            return InverseGaussianDistribution.logPdf(childRate, parentRate, shape);
        }
        else {
            throw new RuntimeException ("Parameter for distribution is not recognised");
        }
    }

    double branchRateSample(double parentRate, double time) {

        double var = variance.getParameterValue(0);

        if (!isEpisodic())
            var *= time;

        //if (isLogSpace) {
        //    final double logParentRate = Math.log(parentRate);

        //    return Math.exp(MathUtils.nextGaussian() * Math.sqrt(var) + logParentRate - (var / 2.));
        //} else {
        //    return MathUtils.nextGaussian() * Math.sqrt(var) + parentRate;
        //}

        if(distribution.equals(LOGNORMAL)) {
            final double logParentRate = Math.log(parentRate);

            return Math.exp(MathUtils.nextGaussian() * Math.sqrt(var) + logParentRate - (var / 2.));
        }
        else if(distribution.equals(NORMAL)) {
            return MathUtils.nextGaussian() * Math.sqrt(var) + parentRate;
        }
        else { /* Inverse Gaussian */
            //return Math.random()
            //Random rand = new Random();
            double lambda = (parentRate * parentRate * parentRate) / var;
            return MathUtils.nextInverseGaussian(parentRate, lambda);
        }
    }

    private Parameter variance;
    //boolean isLogSpace = false;
    String distribution;
}
