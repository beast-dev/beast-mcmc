package dr.inference.operators.shrinkage;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.ExponentialTiltedStableDistribution;
import dr.inference.distribution.shrinkage.JointBayesianBridgeDistributionModel;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.distributions.GammaDistribution;
import dr.xml.*;

import java.util.Arrays;

/**
 * @author Alexander Fisher
 */

public class BayesianBridgePriorSampler {
    public static final String BAYESIAN_BRIDGE_PRIOR_SAMPLER = "bayesianBridgePriorSampler";
    public static final String STEPS = "steps";

    private final Parameter globalScale;
    private final Parameter localScale;
    private final Parameter regressionExponent;
    private final GammaDistribution globalScalePrior;
    private double[] coefficients; // (increments)
    private int steps;

    JointBayesianBridgeDistributionModel bridge;

    public BayesianBridgePriorSampler(JointBayesianBridgeDistributionModel dummyBridge, GammaDistribution globalScalePrior, int N) {
        this.localScale = new Parameter.Default(1, 10.0); // initial value does not matter
        this.globalScale = new Parameter.Default(1, dummyBridge.getGlobalScale().getParameterValue(0));
        this.regressionExponent = dummyBridge.getExponent();

        this.bridge = new JointBayesianBridgeDistributionModel(globalScale, localScale, regressionExponent, dummyBridge.getSlabWidth(), 1, false);

        this.globalScalePrior = globalScalePrior;

        this.steps = N;
        this.coefficients = new double[steps];
        sample(steps);
    }

    void sample(int N) {
        for (int i = 0; i < N; i++) {
            sampleGlobalScale();
            coefficients[i] = MathUtils.nextGaussian() * bridge.getStandardDeviation(0);
            sampleLocalScale(i);
        }
        Arrays.sort(coefficients);
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

    public double getEpsilon(double targetProb) {
        double probToFindEpsilon = (targetProb / 2) + 0.5;
        int i = (int) Math.round(this.steps * probToFindEpsilon);
        double epsilon = coefficients[i];
        if (epsilon < 0.0) {
            throw new RuntimeException("Target probability too small or need more samples from the prior.");
        }
        return coefficients[i];
    }

    public int getSteps() {
        return steps;
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

            GammaDistribution globalScalePrior = null;

            DistributionLikelihood prior = (DistributionLikelihood) xo.getChild(DistributionLikelihood.class);
            if (prior != null) {
                if (prior.getDistribution() instanceof GammaDistribution) {
                    globalScalePrior = (GammaDistribution) prior.getDistribution();
                } else {
                    throw new XMLParseException("Currently only gamma prior on global scale implemented.");
                }
            }

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
                new ElementRule(DistributionLikelihood.class),
                AttributeRule.newIntegerRule(STEPS),
        };
    };
}
