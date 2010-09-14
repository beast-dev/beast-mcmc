package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SelectorOperator;
import dr.xml.*;

/**
 *
 */
public class SelectorOperatorParser extends AbstractXMLObjectParser {
    public static String SELECTOR_OPERATOR = "selectorOperator";

    public String getParserName() {
        return SELECTOR_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final SelectorOperator op = new SelectorOperator(parameter);
        op.setWeight(weight);
        return op;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return ".";
    }

    public Class getReturnType() {
        return SelectorOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                new ElementRule(Parameter.class),
        };
    }
}
