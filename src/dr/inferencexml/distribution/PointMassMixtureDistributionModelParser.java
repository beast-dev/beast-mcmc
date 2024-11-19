package dr.inferencexml.distribution;

import dr.inference.distribution.PointMassMixtureDistributionModel;
import dr.inference.model.Parameter;
import dr.inference.model.CompoundParameter;
import dr.xml.*;

public class PointMassMixtureDistributionModelParser extends AbstractXMLObjectParser {

    public static final String POINT_MASS_MIXTURE_DISTRIBUTION_MODEL = "pointMassMixtureDistributionModel";
    public static final String WEIGHTS = "weights";
    public static final String REALIZED_PARAMETERS = "realizedParameters";
    public static final String WEIGHTS_NORMALIZED = "weightsNormalized";

    public String getParserName() {
        return POINT_MASS_MIXTURE_DISTRIBUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter weights =  (Parameter)xo.getElementFirstChild(WEIGHTS);
        CompoundParameter realizedParameters = (CompoundParameter)xo.getElementFirstChild(REALIZED_PARAMETERS);
        boolean weightsNormalized = xo.getAttribute(WEIGHTS_NORMALIZED, false);

        return new PointMassMixtureDistributionModel(weights, realizedParameters, weightsNormalized);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(WEIGHTS,
                    new XMLSyntaxRule[] { new ElementRule(Parameter.class, false) }),
            new ElementRule(REALIZED_PARAMETERS,
                    new XMLSyntaxRule[] { new ElementRule(CompoundParameter.class, false) }),
            AttributeRule.newBooleanRule(WEIGHTS_NORMALIZED, true),
    };

    public String getParserDescription() {
        return "Describes a mixture of point masses " +
                "that can be used in a distributionLikelihood element";
    }

    public Class getReturnType() {
        return PointMassMixtureDistributionModel.class;
    }

}