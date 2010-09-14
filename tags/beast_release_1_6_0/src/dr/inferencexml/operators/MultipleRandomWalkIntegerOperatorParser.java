package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.MultipleRandomWalkIntegerOperator;
import dr.xml.*;

/**
 */
public class MultipleRandomWalkIntegerOperatorParser extends AbstractXMLObjectParser {
    public static final String MULTIPLE_RANDOM_WALK_INT_OP = "multipleRandomWalkIntegerOperator";
    public static final String SAMPLE_SIZE = "sampleSize";

    public String getParserName() {
        return MULTIPLE_RANDOM_WALK_INT_OP;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        double w = xo.getDoubleAttribute(RandomWalkIntegerOperatorParser.WINDOW_SIZE);
        if (w != Math.floor(w)) {
            throw new XMLParseException("The window size of a randomWalkIntegerOperator should be an integer");
        }

        double s = xo.getDoubleAttribute(SAMPLE_SIZE);
        if (s != Math.floor(s)) {
            throw new XMLParseException("The window size of a randomWalkIntegerOperator should be an integer");
        }

        int windowSize = (int)w;
        int sampleSize = (int)s;
        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        return new MultipleRandomWalkIntegerOperator(parameter, windowSize, sampleSize, weight);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a random walk operator on a given parameter.";
    }

    public Class getReturnType() {
        return MultipleRandomWalkIntegerOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(RandomWalkIntegerOperatorParser.WINDOW_SIZE),
            AttributeRule.newDoubleRule(SAMPLE_SIZE),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(Parameter.class)
    };

}
