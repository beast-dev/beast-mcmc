package dr.inferencexml.model;

import dr.inference.model.CompoundParameter;
import dr.inference.model.JointParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Andrew Rambaut
 */
public class JointParameterParser extends AbstractXMLObjectParser {

    public static final String COMPOUND_PARAMETER = "jointParameter";

    public String getParserName() {
        return COMPOUND_PARAMETER;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        JointParameter jointParameter = new JointParameter((String) null);

        for (int i = 0; i < xo.getChildCount(); i++) {
            jointParameter.addParameter((Parameter) xo.getChild(i));
        }

        return jointParameter;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A parameter that synchronises its component parameters.";
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
        return JointParameter.class;
    }
}
