package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.BitMoveOperator;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 *
 */
public class BitMoveOperatorParser extends AbstractXMLObjectParser {

    public static final String BIT_MOVE_OPERATOR = "bitMoveOperator";
    public static final String NUM_BITS_TO_MOVE = "numBitsToMove";

    public String getParserName() {
        return BIT_MOVE_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        int numBitsToMove = xo.getIntegerAttribute(NUM_BITS_TO_MOVE);

        Parameter bitsParameter = (Parameter) xo.getElementFirstChild("bits");
        Parameter valuesParameter = null;


        if (xo.hasChildNamed("values")) {
            valuesParameter = (Parameter) xo.getElementFirstChild("values");
        }


        return new BitMoveOperator(bitsParameter, valuesParameter, numBitsToMove, weight);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a bit-move operator on a given parameter.";
    }

    public Class getReturnType() {
        return BitMoveOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newIntegerRule(NUM_BITS_TO_MOVE),
            new ElementRule("bits", Parameter.class),
            new ElementRule("values", Parameter.class, "values parameter", true)
    };

}
