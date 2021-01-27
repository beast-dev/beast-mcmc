package dr.inference.operators.shrinkage;

import dr.inference.distribution.ExponentialTiltedStableDistribution;
import dr.inference.distribution.shrinkage.BayesianBridgeDistributionModel;
import dr.inference.distribution.shrinkage.BayesianBridgeStatisticsProvider;
import dr.inference.distribution.shrinkage.JointBayesianBridgeDistributionModel;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.distributions.GammaDistribution;
import dr.math.distributions.NormalDistribution;
import dr.xml.*;

/**
 * @author Alexander Fisher
 */

public class BayesianBridgePriorSampler {
    // todo new plan: empirical cdf by sampling increments from prior
    // todo remove code duplication by accessing relevant functions e.g. jointBayesianBridgeDistribution.getStandardDeviation (this will auto take care of slab/no slab) etc.
    // todo remember to sample from global scale PRIOR
    // inc ~ N(0,1) * sd

    public static final String BAYESIAN_BRIDGE_PRIOR_SAMPLER = "bayesianBridgePriorSampler";
    public static final String STEPS = "steps";

//    private final BayesianBridgeStatisticsProvider provider;
    // todo: make sure these are passed by reference when loaded in the constructor, otherwise need to do the ugly bridge.getGlobalScale().setParameterValue() to use getSD funk
    private final Parameter globalScale;
    private final Parameter localScale;
//    private Parameter slabWidth;
    private final Parameter regressionExponent;
//    private final int dim;
    private final GammaDistribution globalScalePrior;
    private double[] coefficients; // (increments)
    private double sd;
    JointBayesianBridgeDistributionModel bridge;

    public BayesianBridgePriorSampler(JointBayesianBridgeDistributionModel dummyBridge, GammaDistribution globalScalePrior, int N){
//        this.provider = dummyBridge;
        this.globalScale = bridge.getGlobalScale();
        this.localScale = bridge.getLocalScale();
//        this.jointBridge = (bridge instanceof JointBayesianBridgeDistributionModel) ?
//                (JointBayesianBridgeDistributionModel) bridge : null;

        this.bridge = dummyBridge;
        this.regressionExponent = bridge.getExponent();

//        this.globalScale = 1.0;
//        this.localScale = 1.0;
//        this.slabWidth = bridge.getSlabWidth();
//        this.dim = 1;
        this.globalScalePrior = globalScalePrior;

        this.coefficients = new double[N];
        sample(N);
    }

    void sample(int N){
//        double globalTotal = 0.0;
//        double localTotal = 0.0;
        for (int i = 0; i < N; i++){
//            globalTotal += sampleGlobalScale();
//            localTotal += sampleLocalScale();
            sampleGlobalScale();
            coefficients[i] = MathUtils.nextGaussian() * bridge.getStandardDeviation(0);
            sampleLocalScale(i);
        }

//        this.globalScale.setParameterValue(0, globalTotal / N);
//        this.localScale.setParameterValue(0, localTotal / N);
    }

    private void sampleGlobalScale() {

        double priorShape = globalScalePrior.getShape();
        double priorScale = globalScalePrior.getScale();
        double exponent = regressionExponent.getParameterValue(0);

        double phi = GammaDistribution.nextGamma(priorShape, priorScale); // sample

        double draw = Math.pow(phi, -1.0 / exponent); //global scale = phi^(-1/exponent) and phi ~ gamma
        // (phi := nu and global scale := tau in Bayesian Bridge, Polson et al. (2012)

        globalScale.setParameterValue(0, draw);
    }

    private void sampleLocalScale(int i) {

        final double exponent = regressionExponent.getParameterValue(0);
        final double global = globalScale.getParameterValue(0);

        double draw = ExponentialTiltedStableDistribution.nextTiltedStable(
                exponent / 2, Math.pow(coefficients[i] / global, 2)
        );
        localScale.setParameterValueQuietly(0, Math.sqrt(1 / (2 * draw)));
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return BAYESIAN_BRIDGE_PRIOR_SAMPLER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

             JointBayesianBridgeDistributionModel bridge = (JointBayesianBridgeDistributionModel) xo.getChild(JointBayesianBridgeDistributionModel.class);

             if (bridge.getDimension() > 1) {
                 throw new XMLParseException("dim " + bridge.getDimension() + " is not equal to 1. Bayesian Bridge prior sampling not yet implemented for dimensions > 1");
             }

            GammaDistribution globalScalePrior = (GammaDistribution) xo.getChild(GammaDistribution.class);

            int steps = xo.getIntegerAttribute(STEPS);

            BayesianBridgePriorSampler bbPriorSampler = new BayesianBridgePriorSampler(bridge, globalScalePrior, steps);
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
            return "Samples from the prior of a given 'dummy' Bayesian Bridge distribution";
        }

        @Override
        public Class getReturnType() {
            return BayesianBridgePriorSampler.class;
        }
        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(JointBayesianBridgeDistributionModel.class),
                new ElementRule(GammaDistribution.class),
                AttributeRule.newIntegerRule(STEPS),
        };
    };
}
