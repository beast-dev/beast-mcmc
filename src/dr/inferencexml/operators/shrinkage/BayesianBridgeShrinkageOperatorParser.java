package dr.inferencexml.operators.shrinkage;

import dr.evomodel.branchratemodel.AutoCorrelatedBranchRatesDistribution;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.shrinkage.BayesianBridgeStatisticsProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.shrinkage.BayesianBridgeShrinkageOperator;
import dr.math.distributions.GammaDistribution;
import dr.xml.*;

import static dr.evoxml.MaskedPatternsParser.MASK;
import static dr.inference.operators.MCMCOperator.WEIGHT;

public class BayesianBridgeShrinkageOperatorParser extends AbstractXMLObjectParser {

    public final static String BAYESIAN_BRIDGE_PARSER = "bayesianBridgeGibbsOperator";

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

        Parameter mask = null;
        if (xo.hasChildNamed(MASK)) {
            mask = (Parameter) xo.getElementFirstChild(MASK);
        }

        return new BayesianBridgeShrinkageOperator(bayesianBridge, globalScalePrior, mask, weight);
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

            }, true),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return BayesianBridgeShrinkageOperator.class;
    }

    @Override
    public String getParserName() {
        return BAYESIAN_BRIDGE_PARSER;
    }
}
