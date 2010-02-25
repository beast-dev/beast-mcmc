package dr.inferencexml.model;

import dr.inference.model.DuplicatedParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 */
public class DuplicatedParameterParser extends AbstractXMLObjectParser {

    public static final String DUPLICATED_PARAMETER = "duplicatedParameter";
    public static final String COPIES = "copies";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        XMLObject cxo = xo.getChild(COPIES);
        Parameter dup = (Parameter) cxo.getChild(Parameter.class);

        DuplicatedParameter duplicatedParameter = new DuplicatedParameter(parameter);
        duplicatedParameter.addDuplicationParameter(dup);

        return duplicatedParameter;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Parameter.class),
            new ElementRule(COPIES,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }),
    };

    public String getParserDescription() {
        return "A duplicated parameter.";
    }

    public Class getReturnType() {
        return Parameter.class;
    }

    public String getParserName() {
        return DUPLICATED_PARAMETER;
    }
}
