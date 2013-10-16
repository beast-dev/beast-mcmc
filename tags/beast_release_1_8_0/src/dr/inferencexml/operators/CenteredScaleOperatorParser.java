package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.CenteredScaleOperator;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class CenteredScaleOperatorParser extends AbstractXMLObjectParser {

    public static final String CENTERED_SCALE = "centeredScale";
    public static final String SCALE_FACTOR = "scaleFactor";

    public String getParserName() {
        return CENTERED_SCALE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        double scale = xo.getDoubleAttribute(SCALE_FACTOR);
        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        CenteredScaleOperator op = new CenteredScaleOperator(parameter);
        op.setWeight(weight);
        op.scaleFactor = scale;
        op.mode = CoercionMode.parseMode(xo);
        return op;
    }

    public String getParserDescription() {
        return "A centered-scale operator. This operator scales the the values of a multi-dimensional parameter so as to perserve the mean. It does this by expanding or conrtacting the parameter values around the mean.";
    }

    public Class getReturnType() {
        return CenteredScaleOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newDoubleRule(SCALE_FACTOR),
            AttributeRule.newBooleanRule(CoercableMCMCOperator.AUTO_OPTIMIZE, true),
            new ElementRule(Parameter.class)
    };

}
