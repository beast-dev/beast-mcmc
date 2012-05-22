package dr.inferencexml.operators;

import dr.xml.*;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.RandomWalkIntegerOperator;
import dr.inference.operators.RandomWalkIntegerNodeHeightWeightedOperator;
import dr.inference.model.Parameter;

/**
 * @author Chieh-Hsi Wu
 *
 * The parser for random walk integer node height weighted operator.
 */
public class RandomWalkIntegerNodeHeightWeightedOperatorParser extends AbstractXMLObjectParser {

    public static final String RANDOM_WALK_INT_NODE_HEIGHT_WGT_OP = "randomWalkIntegerNodeHeightWeightedOperator";

    public static final String WINDOW_SIZE = "windowSize";
    public static final String INTERNAL_NODE_HEIGHTS = "internalNodeHeights";

    public String getParserName() {
        return RANDOM_WALK_INT_NODE_HEIGHT_WGT_OP;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        double d = xo.getDoubleAttribute(WINDOW_SIZE);
        if (d != Math.floor(d)) {
            throw new XMLParseException("The window size of a " + RANDOM_WALK_INT_NODE_HEIGHT_WGT_OP + " should be an integer");
        }

        int windowSize = (int)d;
        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        Parameter internalNodeHeights = (Parameter)xo.getElementFirstChild(INTERNAL_NODE_HEIGHTS);
        
        return new RandomWalkIntegerNodeHeightWeightedOperator(parameter, windowSize, weight, internalNodeHeights);
    }

    public String getParserDescription() {
        return "This element returns a random walk node height weighted operator on a given parameter.";
    }

    public Class getReturnType() {
        return RandomWalkIntegerNodeHeightWeightedOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(WINDOW_SIZE),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(Parameter.class),
            new ElementRule(INTERNAL_NODE_HEIGHTS, new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
    };
}
