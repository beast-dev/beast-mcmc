package dr.inferencexml.distribution;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.distribution.RandomWalkModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 *
 */
public class RandomWalkModelParser extends AbstractXMLObjectParser {

    public static final String RANDOM_WALK = "randomWalk";
    public static final String LOG_SCALE = "logScale";

    public String getParserName() {
        return RANDOM_WALK;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter data = (Parameter) xo.getChild(Parameter.class);
        ParametricDistributionModel distribution = (ParametricDistributionModel) xo.getChild(ParametricDistributionModel.class);

        boolean logScale = false;
        if (xo.hasAttribute(LOG_SCALE))
            logScale = xo.getBooleanAttribute(LOG_SCALE);

        return new RandomWalkModel(distribution, data, false, logScale);

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newBooleanRule(LOG_SCALE, true),
            new ElementRule(Parameter.class),
            new XORRule(
                    new ElementRule(ParametricDistributionModel.class),
                    new ElementRule(DistributionLikelihood.class)
            )
    };

    public String getParserDescription() {
        return "Describes a first-order random walk. No prior is assumed on the first data element";
    }

    public Class getReturnType() {
        return RandomWalkModel.class;
    }

}
