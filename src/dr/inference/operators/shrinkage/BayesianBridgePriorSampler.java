package dr.inference.operators.shrinkage;

import dr.inference.distribution.ExponentialTiltedStableDistribution;
import dr.inference.distribution.shrinkage.BayesianBridgeStatisticsProvider;
import dr.inference.model.Parameter;
import dr.math.distributions.GammaDistribution;
import dr.xml.*;

public class BayesianBridgePriorSampler {
    // todo: find magic number based on input exponent and slabwidth -- sample from global and local scale
    //todo: it's okay for sampleGLobalScale and sampleLocalScale to return doubles, but they should also update the parameter themselves since
    // sampling one  depends on sampling the other! i.e. need to gibbs sample back and forth.
    // todo: make getSD function

    public static final String BAYESIAN_BRIDGE_PRIOR_SAMPLER = "bayesianBridgePriorSampler";
    public static final String STEPS = "steps";

    private final BayesianBridgeStatisticsProvider provider;
    private final Parameter globalScale;
    private final Parameter localScale;
    private final Parameter regressionExponent;
    private final int dim;
    private final GammaDistribution globalScalePrior;


    public BayesianBridgePriorSampler(BayesianBridgeStatisticsProvider dummyBridge, GammaDistribution globalScalePrior, int N){
        this.provider = dummyBridge;
        this.globalScale = dummyBridge.getGlobalScale(); //shouldn't matter
        this.localScale = dummyBridge.getLocalScale(); //shouldn't matter
        this.regressionExponent = dummyBridge.getExponent();
        this.dim = dummyBridge.getDimension();
        this.globalScalePrior = globalScalePrior;

        if(dim != 1) {
            throw new RuntimeException("dimension of dummy bridge should be 1");
        }

        sample(N);
    }

    void sample(int N){
        double globalTotal = 0.0;
        double localTotal = 0.0;
        for (int i = 0; i < N; i++){
            globalTotal += sampleGlobalScale();
            localTotal += sampleLocalScale();
        }

        this.globalScale.setParameterValue(0, globalTotal / N);
        this.localScale.setParameterValue(0, localTotal / N);

    }

    private double absSumBeta() {

        double exponent = regressionExponent.getParameterValue(0);
        double sum = 0.0;
        for (int i = 0; i < dim; ++i) {
                sum += Math.pow(Math.abs(provider.getCoefficient(i)), exponent);
            }

        return sum;
    }

    private double sampleGlobalScale() {

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

//        globalScale.setParameterValue(0, draw);
        globalScale.setParameterValue(0, draw);
        return(draw);

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

    private double sampleLocalScale() {

        final double exponent = regressionExponent.getParameterValue(0);
        final double global = globalScale.getParameterValue(0);


            double draw = ExponentialTiltedStableDistribution.nextTiltedStable(
                        exponent / 2, Math.pow(provider.getCoefficient(0) / global, 2) // what to do with this getcoefficient???
            );

            double drawFinal = Math.sqrt(1 / (2 * draw));
            localScale.setParameterValue(0, drawFinal);
            return(drawFinal);
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

    public double getStandardDeviation() {
        return 1.0;
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return BAYESIAN_BRIDGE_PRIOR_SAMPLER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            BayesianBridgeStatisticsProvider bbStatistics = (BayesianBridgeStatisticsProvider) xo.getChild(BayesianBridgeStatisticsProvider.class);

            GammaDistribution globalScalePrior = (GammaDistribution) xo.getChild(GammaDistribution.class);

            int steps = xo.getIntegerAttribute(STEPS);

            BayesianBridgePriorSampler bbPriorSampler = new BayesianBridgePriorSampler(bbStatistics, globalScalePrior, steps);
            return bbPriorSampler;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        @Override
        public String getParserDescription() {
            return "Classifies increment as 0 or 1 based on arbitrary density cutoff epsilon.";
        }

        @Override
        public Class getReturnType() {
            return BayesianBridgePriorSampler.class;
        }
        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(BayesianBridgeStatisticsProvider.class),
                new ElementRule(GammaDistribution.class),
                AttributeRule.newIntegerRule(STEPS),
        };
    };
}
