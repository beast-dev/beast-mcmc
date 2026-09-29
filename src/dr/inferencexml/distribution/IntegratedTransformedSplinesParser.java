package dr.inferencexml.distribution;

import dr.inference.model.Parameter;
import dr.math.IntegratedTransformedSplines;
import dr.xml.*;



public class IntegratedTransformedSplinesParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "IntegratedTransformedSplines";
    private static final String COEFFICIENTS = "coefficients";
    private static final String INTERCEPT = "intercept";
    private static final String KNOTS = "knots";
    private static final String DEGREE = "degree";
    private static final String TRANSFORM = "transform";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter coefficients = (Parameter) xo.getElementFirstChild(COEFFICIENTS);
        Parameter intercept = (Parameter) xo.getElementFirstChild(INTERCEPT);
        double[] knots = xo.getDoubleArrayAttribute(KNOTS);
        String transform = xo.getStringAttribute(TRANSFORM);
        int degree = xo.getIntegerAttribute(DEGREE);
        String id = xo.hasId() ? xo.getId() : PARSER_NAME;

        return new IntegratedTransformedSplines(coefficients, intercept, knots, degree, transform);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return IntegratedTransformedSplines.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(COEFFICIENTS,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(INTERCEPT,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            AttributeRule.newDoubleArrayRule(KNOTS),
            AttributeRule.newIntegerRule(DEGREE),
            AttributeRule.newStringRule(TRANSFORM),
    };
}