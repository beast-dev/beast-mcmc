package dr.evomodelxml.coalescent.operators;

import dr.evomodel.coalescent.GMRFSkyrideLikelihood;
import dr.evomodel.coalescent.operators.GMRFSkyrideFixedEffectsGibbsOperator;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.math.distributions.MultivariateDistribution;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.xml.*;

/**
 *
 */
public class GMRFSkyrideFixedEffectsGibbsOperatorParser extends AbstractXMLObjectParser {

    public static final String GMRF_GIBBS_OPERATOR = "gmrfFixedEffectsGibbsOperator";

    public String getParserName() {
        return GMRF_GIBBS_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {


        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        GMRFSkyrideLikelihood gmrfLikelihood = (GMRFSkyrideLikelihood) xo.getChild(GMRFSkyrideLikelihood.class);

        MultivariateDistributionLikelihood likelihood =
                (MultivariateDistributionLikelihood) xo.getChild(MultivariateDistributionLikelihood.class);

        MultivariateDistribution prior = likelihood.getDistribution();

        if (prior.getType().compareTo(MultivariateNormalDistribution.TYPE) != 0)
            throw new XMLParseException("Only a multivariate normal distribution is conjugate for the regression coefficients in a GMRF");

        Parameter param = (Parameter) xo.getChild(Parameter.class);

        return new GMRFSkyrideFixedEffectsGibbsOperator(param,
                gmrfLikelihood, prior, weight);

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a Gibbs operator for regression coefficients in a GMRF.";
    }

    public Class getReturnType() {
        return GMRFSkyrideFixedEffectsGibbsOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
//            return null;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(MultivariateDistributionLikelihood.class),
            new ElementRule(Parameter.class),
            new ElementRule(GMRFSkyrideLikelihood.class)
    };
}
