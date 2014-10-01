package dr.inferencexml.distribution;

import dr.inference.distribution.UniformDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Reads a normal distribution model from a DOM Document element.
 */
public class UniformDistributionModelParser extends AbstractXMLObjectParser {

    public static final String UNIFORM_DISTRIBUTION_MODEL = "uniformDistributionModel";
    public static final String LOWER = "lower";
    public static final String UPPER = "upper";

    public String getParserName() {
        return UNIFORM_DISTRIBUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter upperParam;
        Parameter lowerParam;

        XMLObject cxo = xo.getChild(LOWER);
        if (cxo.getChild(0) instanceof Parameter) {
            lowerParam = (Parameter) cxo.getChild(Parameter.class);
        } else {
            lowerParam = new Parameter.Default(cxo.getDoubleChild(0));
        }

        cxo = xo.getChild(UPPER);
        if (cxo.getChild(0) instanceof Parameter) {
            upperParam = (Parameter) cxo.getChild(Parameter.class);
        } else {
            upperParam = new Parameter.Default(cxo.getDoubleChild(0));
        }

        return new UniformDistributionModel(lowerParam, upperParam);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(LOWER,
                    new XMLSyntaxRule[]{
                            new XORRule(
                                    new ElementRule(Parameter.class),
                                    new ElementRule(Double.class)
                            )}
            ),
            new ElementRule(UPPER,
                    new XMLSyntaxRule[]{
                            new XORRule(
                                    new ElementRule(Parameter.class),
                                    new ElementRule(Double.class)
                            )}
            )
    };

    public String getParserDescription() {
        return "Describes a uniform distribution with a given lower and upper bounds ";
    }

    public Class getReturnType() {
        return UniformDistributionModel.class;
    }
}
