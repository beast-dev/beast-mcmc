package dr.inferencexml.operators;

import dr.xml.*;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.RandomWalkIntegerSetSizeWeightedOperator;
import dr.inference.model.Parameter;
import dr.evomodel.tree.MicrosatelliteSamplerTreeModel;

/**
 * @author Chieh-Hsi Wu
 *
 * This is the parser for the random walk integer set size weighted operator.
 */
public class RandomWalkIntegerSetSizeWeightedOperatorParser extends AbstractXMLObjectParser {

    public static final String RANDOM_WALK_INT_SET_SIZE_WGT_OP = "randomWalkIntegerSetSizeWeightedOperator";

    public static final String WINDOW_SIZE = "windowSize";
    public static final String BASE_SET_SIZE = "baseSetSize";

    public String getParserName() {
        return RANDOM_WALK_INT_SET_SIZE_WGT_OP;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        double d = xo.getDoubleAttribute(WINDOW_SIZE);
        if (d != Math.floor(d)) {
            throw new XMLParseException("The window size of a " + RANDOM_WALK_INT_SET_SIZE_WGT_OP + " should be an integer");
        }

        double baseSetSize = xo.getDoubleAttribute(BASE_SET_SIZE);

        int windowSize = (int)d;
        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        MicrosatelliteSamplerTreeModel msatSampleTreeModel = (MicrosatelliteSamplerTreeModel)xo.getChild(MicrosatelliteSamplerTreeModel.class);

        return new RandomWalkIntegerSetSizeWeightedOperator(parameter, windowSize, weight, msatSampleTreeModel, baseSetSize);
    }

    public String getParserDescription() {
        return "This element returns a random walk set size weighted operator on a given parameter.";
    }

    public Class getReturnType() {
        return RandomWalkIntegerSetSizeWeightedOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(WINDOW_SIZE),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newDoubleRule(BASE_SET_SIZE),
            new ElementRule(Parameter.class),
            new ElementRule(MicrosatelliteSamplerTreeModel.class)
    };
}