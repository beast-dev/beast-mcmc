package dr.inferencexml.distribution;

import dr.inference.model.Parameter;
import dr.math.SquaredCachedSplines;
import dr.xml.*;

/**
 * Parses a SquaredCachedSplines model (squared B-spline rate function with Gram-matrix cache).
 *
 * XML element: squaredCachedSplines
 *
 * Example:
 *   <squaredCachedSplines id="splineRates"
 *       lowerBoundary="0" upperBoundary="60000" degree="3"
 *       knots="5000 10000 15000 20000 25000 30000 35000 40000 45000 50000 55000">
 *       <coefficients><parameter id="spline.coefficients" .../></coefficients>
 *       <intercept><parameter id="spline.intercept" value="0.005"/></intercept>
 *   </squaredCachedSplines>
 */
public class SquaredCachedSplinesParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME    = "squaredCachedSplines";
    private static final String COEFFICIENTS  = "coefficients";
    private static final String INTERCEPT     = "intercept";
    private static final String KNOTS         = "knots";
    private static final String LOWER_BOUNDARY = "lowerBoundary";
    private static final String UPPER_BOUNDARY = "upperBoundary";
    private static final String DEGREE        = "degree";

    @Override
    public String getParserName() { return PARSER_NAME; }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Parameter coefficients  = (Parameter) xo.getElementFirstChild(COEFFICIENTS);
        Parameter intercept     = (Parameter) xo.getElementFirstChild(INTERCEPT);
        double[]  knots         = xo.getDoubleArrayAttribute(KNOTS);
        double    lowerBoundary = xo.getDoubleAttribute(LOWER_BOUNDARY);
        double    upperBoundary = xo.getDoubleAttribute(UPPER_BOUNDARY);
        int       degree        = xo.getIntegerAttribute(DEGREE);

        return new SquaredCachedSplines(coefficients, intercept, knots,
                lowerBoundary, upperBoundary, degree);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    @Override
    public String getParserDescription() {
        return "Squared B-spline rate model with interval-keyed Gram matrix cache.";
    }

    @Override
    public Class getReturnType() { return SquaredCachedSplines.class; }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(COEFFICIENTS, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(INTERCEPT,    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            AttributeRule.newDoubleArrayRule(KNOTS),
            AttributeRule.newDoubleRule(LOWER_BOUNDARY),
            AttributeRule.newDoubleRule(UPPER_BOUNDARY),
            AttributeRule.newIntegerRule(DEGREE),
    };
}
