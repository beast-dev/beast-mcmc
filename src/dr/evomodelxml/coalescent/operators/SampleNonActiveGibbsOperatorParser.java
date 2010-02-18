package dr.evomodelxml.coalescent.operators;

import dr.evomodel.coalescent.operators.SampleNonActiveGibbsOperator;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class SampleNonActiveGibbsOperatorParser extends AbstractXMLObjectParser {
    public static String SAMPLE_NONACTIVE_GIBBS_OPERATOR = "sampleNonActiveOperator";
    public static String DISTRIBUTION = "distribution";

    public static String INDICATOR_PARAMETER = "indicators";
    public static String DATA_PARAMETER = "data";

    public String getParserName() {
        return SAMPLE_NONACTIVE_GIBBS_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        XMLObject cxo = xo.getChild(DISTRIBUTION);
        ParametricDistributionModel distribution =
                (ParametricDistributionModel) cxo.getChild(ParametricDistributionModel.class);

        cxo = xo.getChild(DATA_PARAMETER);
        Parameter data = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(INDICATOR_PARAMETER);
        Parameter indicators = (Parameter) cxo.getChild(Parameter.class);

        return new SampleNonActiveGibbsOperator(distribution, data, indicators, weight);

    }

    // ************************************************************************
    // AbstractXMLObjectParser implementation
    // ************************************************************************

    public String getParserDescription() {
        return "This element returns a Gibbs operator for the joint distribution of population sizes.";
    }

    public Class getReturnType() {
        return SampleNonActiveGibbsOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule("distribution",
                    new XMLSyntaxRule[]{new ElementRule(ParametricDistributionModel.class)}),
            new ElementRule(INDICATOR_PARAMETER,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(DATA_PARAMETER,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
    };

}
