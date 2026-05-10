package dr.inferencexml.distribution;

import dr.inference.model.Parameter;
import dr.math.IntegratedTransformedSplines;
import dr.xml.*;



public class IntegratedTransformedSplinesParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "IntegratedTransformedSplines";
    private static final String COEFFICIENTS = "coefficients";
    private static final String INTERCEPT = "intercept";
    private static final String KNOTS = "knots";
    private static final String LOWERBOUNDARY = "lowerBoundary";
    private static final String UPPERBOUNDARY = "upperBoundary";
    private static final String DEGREE = "degree";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter coefficients = (Parameter) xo.getElementFirstChild(COEFFICIENTS);
        Parameter intercept = null;
        if (xo.hasChildNamed(INTERCEPT)) {
            intercept = (Parameter) xo.getElementFirstChild(INTERCEPT);
        }
        double[] knots = xo.getDoubleArrayAttribute(KNOTS);
        Double lowerBoundary = xo.getDoubleAttribute(LOWERBOUNDARY);
        Double upperBoundary = xo.getDoubleAttribute(UPPERBOUNDARY);
        int degree = xo.getIntegerAttribute(DEGREE);

        return new IntegratedTransformedSplines(coefficients, intercept, knots, lowerBoundary, upperBoundary, degree);
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
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            AttributeRule.newDoubleArrayRule(KNOTS),
            AttributeRule.newDoubleRule(LOWERBOUNDARY),
            AttributeRule.newDoubleRule(UPPERBOUNDARY),
            AttributeRule.newIntegerRule(DEGREE),
    };
}