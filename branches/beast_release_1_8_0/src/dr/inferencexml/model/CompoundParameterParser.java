package dr.inferencexml.model;

import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 */
public class CompoundParameterParser extends AbstractXMLObjectParser {

    public static final String COMPOUND_PARAMETER = "compoundParameter";

    public String getParserName() {
        return COMPOUND_PARAMETER;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        CompoundParameter compoundParameter = new CompoundParameter((String) null);

        for (int i = 0; i < xo.getChildCount(); i++) {
            compoundParameter.addParameter((Parameter) xo.getChild(i));
        }

        return compoundParameter;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A multidimensional parameter constructed from its component parameters.";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules;{
        rules = new XMLSyntaxRule[]{
                new ElementRule(Parameter.class, 1, Integer.MAX_VALUE),
        };
    }

    public Class getReturnType() {
        return CompoundParameter.class;
    }
}
