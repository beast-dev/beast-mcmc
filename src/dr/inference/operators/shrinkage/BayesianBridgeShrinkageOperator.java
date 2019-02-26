package dr.inference.operators.shrinkage;

import dr.inference.distribution.shrinkage.BayesianBridgeLikelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.distributions.GammaDistribution;

import static dr.inferencexml.operators.shrinkage.BayesianBridgeShrinkageOperatorParser.BAYESIAN_BRIDGE_PARSER;

/**
 * @author Marc A. Suchard
 */

public class BayesianBridgeShrinkageOperator extends SimpleMCMCOperator implements GibbsOperator {

    private final Parameter globalScale;
    private final Parameter localScale;
    private final Parameter regressionExponent;

    private final GammaDistribution globalScalePrior;


    public BayesianBridgeShrinkageOperator(BayesianBridgeLikelihood bayesianBridge,
                                           GammaDistribution globalScalePrior,
                                           double weight) {
        setWeight(weight);

        this.globalScale = bayesianBridge.getGlobalScale();
        this.localScale = bayesianBridge.getLocalScale();
        this.regressionExponent = bayesianBridge.getExponent();

        this.globalScalePrior = globalScalePrior;

    }

    @Override
    public String getOperatorName() {
        return BAYESIAN_BRIDGE_PARSER;
    }

    @Override
    public double doOperation() {

        sampleGlobalScale(); // Order matters

        if (localScale != null) {
            sampleLocalScale();
        }

        return 0;
    }

    private void sampleGlobalScale() {

        double priorShape = globalScalePrior.getShape();
        double priorScale = globalScalePrior.getScale();
        double exponent = regressionExponent.getParameterValue(0);

        double shape = localScale.getDimension() / exponent;
        double rate = absSumBeta();

        if (priorShape > 0.0) {
            shape += priorShape;
            rate += 1.0 / priorScale;
        }

        double phi = GammaDistribution.nextGamma(shape, 1.0 / rate);
        double draw = Math.pow(phi, -1.0 / exponent);

        globalScale.setParameterValue(0, draw);
    }

    private double absSumBeta() {

        double exponent = regressionExponent.getParameterValue(0);
        double sum = 0.0;
        for (int i = 0; i < localScale.getDimension(); ++i) {
            sum += Math.pow(Math.abs(localScale.getParameterValue(i)), exponent);
        }

        return sum;
    }

    private void sampleLocalScale() {
        // TODO
    }
}

//    if self.prior_type['global_shrinkage'] == 'jeffreys':
//
//        # Conjugate update for phi = 1 / gshrink ** reg_exponent
//        shape = beta_with_shrinkage.size / reg_exponent
//        scale = 1 / np.sum(np.abs(beta_with_shrinkage) ** reg_exponent)
//        phi = self.rg.np_random.gamma(shape, scale=scale)
//        gshrink = 1 / phi ** (1 / reg_exponent)
//
//    elif self.prior_type['global_shrinkage'] == 'half-cauchy':
//
//        gshrink = self.slice_sample_global_shrinkage(
//            gshrink, beta_with_shrinkage, self.prior_param['global_shrinkage']['scale'], reg_exponent
//        )
//



//    To update the global scale parameter τ , we work directly with the exponential-power density, marginalizing out the latent variables {ωj,uj}. This is a crucial source of efficiency in the bridge MCMC, and leads to the favorable mixing evident in Figure 1. From (1), observe that the posterior for ν ≡ τ−α, given β, is conditionally independent of y, and takes the form
//    p
//    p(ν | β) ∝ νp/α exp(−ν 􏰑|βj|α) p(ν).
//    j=1
//    Therefore if ν has a Gamma(c, d) prior, its conditional posterior will also be a gamma dis- tribution, with hyperparameters c⋆ = c+p/α and d⋆ = d+􏰏pj=1 |βj|α. To sample τ, simply draw ν from this gamma distribution, and use the transformation τ = ν−1/α.

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


//       def update_local_shrinkage(self, gshrink, beta_with_shrinkage, reg_exponent):
//
//        lshrink_sq = 1 / np.array([
//            2 * self.rg.tilted_stable(reg_exponent / 2, (beta_j / gshrink) ** 2)
//            for beta_j in beta_with_shrinkage
//        ])
//        lshrink = np.sqrt(lshrink_sq)
//
//        # TODO: Pick the lower and upper bound more carefully.
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