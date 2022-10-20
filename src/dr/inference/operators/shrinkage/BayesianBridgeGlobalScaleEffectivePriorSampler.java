package dr.inference.operators.shrinkage;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.shrinkage.BayesianBridgeRNG;
import dr.inference.distribution.shrinkage.JointBayesianBridgeDistributionModel;
import dr.inference.model.Statistic;
import dr.math.distributions.GammaDistribution;
import dr.xml.*;

/**
 * @author Andy Magee
 * @author Yucai Shao
 */

public class BayesianBridgeGlobalScaleEffectivePriorSampler extends Statistic.Abstract implements Reportable {
    public static final String PRIOR_SAMPLER = "bayesianBridgeGlobalScaleEffectivePriorSampler";

    JointBayesianBridgeDistributionModel bridge;
    private final GammaDistribution globalScalePrior;

    public BayesianBridgeGlobalScaleEffectivePriorSampler(JointBayesianBridgeDistributionModel tempBridge, GammaDistribution globalScalePrior) {
        this.bridge = new JointBayesianBridgeDistributionModel(tempBridge.getGlobalScale(), tempBridge.getLocalScale(), tempBridge.getExponent(), tempBridge.getSlabWidth(), tempBridge.getDimension(), false);
        this.globalScalePrior = globalScalePrior;
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

    private double sampleGlobalScalePrior() {
        double priorShape = globalScalePrior.getShape();
        double priorScale = globalScalePrior.getScale();
        double exponent = bridge.getExponent().getParameterValue(0);

        double phi = GammaDistribution.nextGamma(priorShape, priorScale); // sample

        double draw = Math.pow(phi, -1.0 / exponent); //global scale = phi^(-1/exponent) and phi ~ gamma
        // (phi := nu and global scale := tau in Bayesian Bridge, Polson et al. (2012)
        // bridge.setGlobalScale(draw);
        return draw;
    }

    private double absSumBeta(double[] betaDraw) {

        double exponent = bridge.getExponent().getParameterValue(0);
        double sum = 0.0;
        for (int i = 0; i < bridge.getDimension(); ++i) {
            if (true) {
                sum += Math.pow(Math.abs(betaDraw[i]), exponent);
            }
        }
        return sum;
    }
    private double sampleGlobalScale(double[] betaDraw) {

        double priorShape = globalScalePrior.getShape();
        double priorScale = globalScalePrior.getScale();
        double exponent = bridge.getExponent().getParameterValue(0);
        double effectiveDim = bridge.getDimension();
        double shape = effectiveDim / exponent;
        double rate = absSumBeta(betaDraw);

        if (priorShape > 0.0) {
            shape += priorShape;
            rate += 1.0 / priorScale;
        }

        double phi = GammaDistribution.nextGamma(shape, 1.0 / rate);
        double globalScaleDraw = Math.pow(phi, -1.0 / exponent);
        return globalScaleDraw;
    }

    @Override
    public double getStatisticValue(int dim) {
        double globalScaleDraw = sampleGlobalScalePrior();
        double[] betaDraw = bridge.nextRandom(globalScaleDraw);
        return sampleGlobalScale(betaDraw);
    }
}
