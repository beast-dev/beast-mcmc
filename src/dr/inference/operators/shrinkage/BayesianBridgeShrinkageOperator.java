package dr.inference.operators.shrinkage;

import dr.inference.distribution.ExponentialTiltedStableDistribution;
import dr.inference.distribution.shrinkage.BayesianBridgeStatisticsProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.distributions.GammaDistribution;

import static dr.inferencexml.operators.shrinkage.BayesianBridgeShrinkageOperatorParser.BAYESIAN_BRIDGE_PARSER;

/**
 * @author Marc A. Suchard
 * @author Akihiko Nishimura
 */

public class BayesianBridgeShrinkageOperator extends SimpleMCMCOperator implements GibbsOperator {

    private final BayesianBridgeStatisticsProvider provider;
    private final Parameter globalScale;
    private final Parameter localScale;
    private final Parameter regressionExponent;
    private final int dim;

    private final GammaDistribution globalScalePrior;


    public BayesianBridgeShrinkageOperator(BayesianBridgeStatisticsProvider bayesianBridge,
                                           GammaDistribution globalScalePrior,
                                           double weight) {
        setWeight(weight);

        this.provider = bayesianBridge;
        this.globalScale = bayesianBridge.getGlobalScale();
        this.localScale = bayesianBridge.getLocalScale();
        this.regressionExponent = bayesianBridge.getExponent();
        this.dim = bayesianBridge.getDimension();

        this.globalScalePrior = globalScalePrior;
    }

    @Override
    public String getOperatorName() {
        return BAYESIAN_BRIDGE_PARSER;
    }

    @Override
    public double doOperation() {

        if (globalScalePrior != null) {
            sampleGlobalScale(); // Order matters
        }

        if (localScale != null) {
            sampleLocalScale();
        }

        return 0;
    }

    private void sampleGlobalScale() {

        double priorShape = globalScalePrior.getShape();
        double priorScale = globalScalePrior.getScale();
        double exponent = regressionExponent.getParameterValue(0);

        double shape = dim / exponent;
        double rate = absSumBeta();

        if (priorShape > 0.0) {
            shape += priorShape;
            rate += 1.0 / priorScale;
        }

        double phi = GammaDistribution.nextGamma(shape, 1.0 / rate);
        double draw = Math.pow(phi, -1.0 / exponent);

        globalScale.setParameterValue(0, draw);

        //        # Conjugate update for phi = 1 / gshrink ** reg_exponent
        //        shape = beta_with_shrinkage.size / reg_exponent
        //        scale = 1 / np.sum(np.abs(beta_with_shrinkage) ** reg_exponent)
        //        phi = self.rg.np_random.gamma(shape, scale=scale)
        //        gshrink = 1 / phi ** (1 / reg_exponent)
        //
        //   To update the global scale parameter τ, we work directly with the exponential-power density, marginalizing out the latent variables {ωj,uj}. This is a crucial source of efficiency in the bridge MCMC, and leads to the favorable mixing evident in Figure 1. From (1), observe that the posterior for ν ≡ τ−α, given β, is conditionally independent of y, and takes the form
        //
        //    p(ν | β) propto νp/α exp(−ν |βj|α) p(ν).
        //    j=1
        //    Therefore if ν has a Gamma(c, d) prior, its conditional posterior will also be a gamma distribution, with hyperparameters c⋆ = c+p/α and d⋆ = d+pj=1 |βj|α. To sample τ, simply draw ν from this gamma distribution, and use the transformation τ = ν−1/α.
    }

    private double absSumBeta() {

        double exponent = regressionExponent.getParameterValue(0);
        double sum = 0.0;
        for (int i = 0; i < dim; ++i) {
            sum += Math.pow(Math.abs(provider.getCoefficient(i)), exponent);
        }

        return sum;
    }

    private void sampleLocalScale() {

        final double exponent = regressionExponent.getParameterValue(0);
        final double global = globalScale.getParameterValue(0);

        for (int i = 0; i < dim; ++i) {
            double draw = ExponentialTiltedStableDistribution.nextTiltedStable(
                    exponent / 2, Math.pow(provider.getCoefficient(i) / global, 2)
            );

            localScale.setParameterValueQuietly(i, Math.sqrt(1 / (2 * draw)));
        }

        localScale.fireParameterChangedEvent();

        //
        //        lshrink_sq = 1 / np.array([
        //            2 * self.rg.tilted_stable(reg_exponent / 2, (beta_j / gshrink) ** 2)
        //            for beta_j in beta_with_shrinkage
        //        ])
        //        lshrink = np.sqrt(lshrink_sq)

        //        if np.any(lshrink == 0):
        //            warn_message_only(
        //                "Local shrinkage parameter under-flowed. Replacing with a small number.")
        //            lshrink[lshrink == 0] = 10e-16
        //        elif np.any(np.isinf(lshrink)):
        //            warn_message_only(
        //                "Local shrinkage parameter under-flowed. Replacing with a large number.")
        //            lshrink[np.isinf(lshrink)] = 2.0 / gshrink
        //
        //        return lshrink
    }
}

//    def slice_sample_global_shrinkage(
//            self, gshrink, beta_with_shrinkage, global_scale, reg_exponent):
//        """ Slice sample phi = 1 / gshrink ** reg_exponent. """
//
//        n_update = 10 # Slice sample for multiple iterations to ensure good mixing.
//
//        # Initialize a gamma distribution object.
//        shape = (beta_with_shrinkage.size + 1) / reg_exponent
//        scale = 1 / np.sum(np.abs(beta_with_shrinkage) ** reg_exponent)
//        gamma_rv = sp.stats.gamma(shape, scale=scale)
//
//        phi = 1 / gshrink
//        for i in range(n_update):
//            u = self.rg.np_random.uniform() \
//                / (1 + (global_scale * phi ** (1 / reg_exponent)) ** 2)
//            upper = (np.sqrt(1 / u - 1) / global_scale) ** reg_exponent
//                # Invert the half-Cauchy density.
//            phi = gamma_rv.ppf(gamma_rv.cdf(upper) * self.rg.np_random.uniform())
//            if np.isnan(phi):
//                # Inverse CDF method can fail if the current conditional
//                # distribution is drastically different from the previous one.
//                # In this case, ignore the slicing variable and just sample from
//                # a Gamma.
//                phi = gamma_rv.rvs()
//        gshrink = 1 / phi ** (1 / reg_exponent)
//
//        return gshrink
