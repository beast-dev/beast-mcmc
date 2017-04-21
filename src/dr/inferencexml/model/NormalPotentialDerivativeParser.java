package dr.inferencexml.model;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.NormalPotentialDerivative;
import dr.inference.model.Parameter;
import dr.math.distributions.Distribution;
import dr.xml.*;

/**
 * @author Max Tolkoff
 */
public class NormalPotentialDerivativeParser extends AbstractXMLObjectParser{
    public static final String NORMAL_POTENTIAL_DERIVATIVE = "normalPotentialDerivative";

    public static final String MEAN = "mean";
    public static final String STDEV = "stdev";

    @Override
    public String getParserName() {
        return NORMAL_POTENTIAL_DERIVATIVE;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        double mean = xo.getDoubleAttribute(MEAN);
        double stdev = xo.getDoubleAttribute(STDEV);


        return new NormalPotentialDerivative(mean, stdev, parameter);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MEAN),
            AttributeRule.newDoubleRule(STDEV),
            new ElementRule(Parameter.class),
    };


    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return NormalPotentialDerivative.class;
    }
}
