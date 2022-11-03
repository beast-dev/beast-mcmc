package dr.inference.operators.shrinkage;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.shrinkage.BayesianBridgeRNG;
import dr.inference.distribution.shrinkage.JointBayesianBridgeDistributionModel;
import dr.inference.model.Statistic;
import dr.math.distributions.GammaDistribution;
import dr.inference.operators.shrinkage.BayesianBridgeShrinkageOperator;
import dr.xml.*;

/**
 * @author Andy Magee
 * @author Yucai Shao
 */

public class BayesianBridgeGlobalScaleEffectivePriorSampler extends Statistic.Abstract implements Reportable {
    public static final String PRIOR_SAMPLER = "bayesianBridgeGlobalScaleEffectivePriorSampler";
    BayesianBridgeShrinkageOperator operator;
    JointBayesianBridgeDistributionModel bridge;

    private final GammaDistribution globalScalePrior;

    public BayesianBridgeGlobalScaleEffectivePriorSampler(JointBayesianBridgeDistributionModel tempBridge, GammaDistribution globalScalePrior) {
        this.bridge = new JointBayesianBridgeDistributionModel(tempBridge.getGlobalScale(), tempBridge.getLocalScale(), tempBridge.getExponent(), tempBridge.getSlabWidth(), tempBridge.getDimension(), false);
        this.globalScalePrior = globalScalePrior;
        this.operator = new BayesianBridgeShrinkageOperator(bridge, globalScalePrior, null, Double.MIN_VALUE);
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return PRIOR_SAMPLER;
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

            BayesianBridgeGlobalScaleEffectivePriorSampler BBEGSPS = new BayesianBridgeGlobalScaleEffectivePriorSampler(bridge, globalScalePrior);
            return BBEGSPS;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        @Override
        public String getParserDescription() {
            return "Samples from the \"effective prior\" of the global scale for a shrunken-shoulder-regularized Bayesian Bridge distribution.";
        }

        @Override
        public Class getReturnType() {
            return BayesianBridgeGlobalScaleEffectivePriorSampler.class;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(JointBayesianBridgeDistributionModel.class),
                new ElementRule(DistributionLikelihood.class),
        };
    };

    //************************************************************************
    // Reportable implementation
    //************************************************************************

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();
        sb.append("Using bridge with following parameters:").append("\n");
        sb.append("  globalScale named ").append(bridge.getGlobalScale().getParameterName()).append(" and current value ").append(bridge.getGlobalScale().getParameterValue(0)).append("\n");
        sb.append("  localScale named ").append(bridge.getLocalScale().getParameterName()).append(" and current value ").append(new dr.math.matrixAlgebra.Vector(bridge.getLocalScale().getParameterValues())).append("\n");
        sb.append("  exponent named ").append(bridge.getExponent().getParameterName()).append(" and current value ").append(bridge.getExponent().getParameterValue(0)).append("\n");
        sb.append("  slabWidth named ").append(bridge.getSlabWidth().getParameterName()).append(" and current value ").append(bridge.getSlabWidth().getParameterValue(0)).append("\n");
        sb.append("Bridge gamma prior has shape with current value ").append(globalScalePrior.getShape()).append("and scale with current value ").append(globalScalePrior.getScale()).append("\n");
        return sb.toString();
    }

    //************************************************************************
    // Statistic.Abstract implementation
    //************************************************************************
    @Override
    public int getDimension() {
        return 1;
    }

    // TODO we could consider making a class for this, e.g. PowerTransformedGamma
    // If we then let this be an additional option for the Gibbs operator, this would let us specify the correct prior to be logged in the jointDensity
    private double sampleGlobalScalePrior() {
        double priorShape = globalScalePrior.getShape();
        double priorScale = globalScalePrior.getScale();
        double exponent = bridge.getExponent().getParameterValue(0);

        double phi = GammaDistribution.nextGamma(priorShape, priorScale);

        double draw = Math.pow(phi, -1.0 / exponent);
        return draw;
    }

    private double[] nextRandom() {
        double[] draws;
        double globalScaleDraw = sampleGlobalScalePrior();
        if ( bridge.getSlabWidth() != null) {
            draws = BayesianBridgeRNG.nextRandom(globalScaleDraw, bridge.getExponent().getParameterValue(0), bridge.getSlabWidth().getParameterValue(0), bridge.getDimension());
        } else {
            draws = BayesianBridgeRNG.nextRandom(globalScaleDraw, bridge.getExponent().getParameterValue(0), bridge.getDimension());
        }
        return draws;
    }

    @Override
    public double getStatisticValue(int dim) {
        double exponent = bridge.getExponent().getParameterValue(0);
        double[] betaDraw = nextRandom();

        double sum = 0.0;
        for (int i = 0; i < bridge.getDimension(); ++i) {
            sum += Math.pow(betaDraw[i], 2);
        }

        return sum;
    }
}
