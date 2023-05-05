package dr.inferencexml.operators.shrinkage;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.shrinkage.BayesianBridgeStatisticsProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.shrinkage.BayesianBridgeShrinkageOperator;
import dr.inference.operators.shrinkage.DimensionMismatchedBayesianBridgeShrinkageOperator;
import dr.math.distributions.GammaDistribution;
import dr.xml.*;

import static dr.evoxml.MaskedPatternsParser.MASK;
import static dr.inference.operators.MCMCOperator.WEIGHT;

public class DimensionMismatchedBayesianBridgeShrinkageOperatorParser extends AbstractXMLObjectParser {

    public final static String BAYESIAN_BRIDGE_PARSER = "dimensionMismatchedBayesianBridgeGibbsOperator";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(WEIGHT);

        BayesianBridgeStatisticsProvider bayesianBridge =
                (BayesianBridgeStatisticsProvider) xo.getChild(BayesianBridgeStatisticsProvider.class);


        GammaDistribution globalScalePrior = null;

        // This prior is actually on phi = globalScale^-exponent
        DistributionLikelihood prior = (DistributionLikelihood) xo.getChild(DistributionLikelihood.class);
        if (prior != null) {
            if (prior.getDistribution() instanceof GammaDistribution) {
                globalScalePrior = (GammaDistribution) prior.getDistribution();
            } else {
                throw new XMLParseException("Gibbs sampler only implemented for a gamma prior on globalScale^(-exponent).");
            }
        }

        Parameter mask = (Parameter) xo.getElementFirstChild(MASK);

        BayesianBridgeShrinkageOperator operator = new BayesianBridgeShrinkageOperator(bayesianBridge, globalScalePrior, mask, Double.MIN_VALUE);

        return new DimensionMismatchedBayesianBridgeShrinkageOperator(bayesianBridge, globalScalePrior, mask, operator, weight);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(WEIGHT),
            new ElementRule(BayesianBridgeStatisticsProvider.class),
            new ElementRule(DistributionLikelihood.class, true),
            new ElementRule(MASK, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class),

            }, false),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return DimensionMismatchedBayesianBridgeShrinkageOperator.class;
    }

    @Override
    public String getParserName() {
        return BAYESIAN_BRIDGE_PARSER;
    }
}
