package dr.inferencexml.model;

import dr.inference.model.MaskedParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 */
public class MaskedParameterParser extends AbstractXMLObjectParser {

    public static final String MASKED_PARAMETER = "maskedParameter";
    public static final String MASKING = "mask";
    public static final String COMPLEMENT = "complement";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        XMLObject cxo = xo.getChild(MASKING);
        Parameter mask = (Parameter) cxo.getChild(Parameter.class);

        if (mask.getDimension() != parameter.getDimension())
            throw new XMLParseException("dim(" + parameter.getId() + ") != dim(" + mask.getId() + ")");

        MaskedParameter maskedParameter = new MaskedParameter(parameter);
        boolean ones = !xo.getAttribute(COMPLEMENT, false);
        maskedParameter.addMask(mask, ones);

        return maskedParameter;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Parameter.class),
            new ElementRule(MASKING,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }),
            AttributeRule.newBooleanRule(COMPLEMENT, true),
    };

    public String getParserDescription() {
        return "A masked parameter.";
    }

    public Class getReturnType() {
        return Parameter.class;
    }

    public String getParserName() {
        return MASKED_PARAMETER;
    }
}
