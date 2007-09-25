package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * A generic operator that moves k 1 bits to k zero locations.
 *
 * @author Alexei Drummond
 * @version $Id$
 */
public class BitMoveOperator extends SimpleMCMCOperator {

    public static final String BIT_MOVE_OPERATOR = "bitMoveOperator";
    public static final String NUM_BITS_TO_MOVE = "numBitsToMove";


    public BitMoveOperator(Parameter bitsParameter, Parameter valuesParameter, int numBitsToMove, double weight) {
        this.bitsParameter = bitsParameter;
        this.valuesParameter = valuesParameter;

        if (valuesParameter != null && bitsParameter.getDimension() != valuesParameter.getDimension()) {
            throw new IllegalArgumentException("bits parameter must be same length as values parameter");
        }

        this.numBitsToMove = numBitsToMove;
        setWeight(weight);
    }

    /**
     * Pick a random k ones in the vector and move them to a random k zero positions.
     */
    public final double doOperation() throws OperatorFailedException {
        final int dim = bitsParameter.getDimension();
        List<Integer> ones = new ArrayList<Integer>();
        List<Integer> zeros = new ArrayList<Integer>();

        for (int i = 0; i < dim; i++) {
            if (bitsParameter.getParameterValue(i) == 1.0) {
                ones.add(i);
            } else {
                zeros.add(i);
            }
        }

        if (ones.size() >= numBitsToMove && zeros.size() >= numBitsToMove) {

            for (int i = 0; i < numBitsToMove; i++) {

                int myOne = ones.remove(MathUtils.nextInt(ones.size()));
                int myZero = zeros.remove(MathUtils.nextInt(zeros.size()));

                bitsParameter.setParameterValue(myOne, 0.0);
                bitsParameter.setParameterValue(myZero, 1.0);

                if (valuesParameter != null) {
                    double value1 = valuesParameter.getParameterValue(myOne);
                    double value2 = valuesParameter.getParameterValue(myZero);
                    valuesParameter.setParameterValue(myOne, value2);
                    valuesParameter.setParameterValue(myZero, value1);
                }

            }
        } else throw new OperatorFailedException("Not enough bits to move!");

        return 0.0;
    }

    // Interface MCMCOperator
    public final String getOperatorName() {
        return "bitMove(" + bitsParameter.getParameterName() + ", " + numBitsToMove + ")";
    }

    public final String getPerformanceSuggestion() {

        return "no performance suggestion";
    }

    public String toString() {
        return getOperatorName();
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return BIT_MOVE_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);

            int numBitsToMove = xo.getIntegerAttribute(NUM_BITS_TO_MOVE);

            Parameter bitsParameter = (Parameter) xo.getSocketChild("bits");
            Parameter valuesParameter = null;


            if (xo.hasSocket("values")) {
                valuesParameter = (Parameter) xo.getSocketChild("values");
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
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(WEIGHT),
                AttributeRule.newIntegerRule(NUM_BITS_TO_MOVE),
                new ElementRule("bits", Parameter.class),
                new ElementRule("values", Parameter.class, "values parameter", true)
        };

    };
    // Private instance variables

    private Parameter bitsParameter = null;
    private Parameter valuesParameter = null;
    private int numBitsToMove = 1;
}
