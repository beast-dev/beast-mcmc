package dr.inference.operators;

import dr.xml.*;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

/**
 * @author Marc Suchard
 */
public class RateBitExchangeOperator extends SimpleMCMCOperator {

    public static final String OPERATOR_NAME = "rateBitExchangeOperator";
    public static final String BITS = "bits";
    public static final String RATES = "rates";

    RateBitExchangeOperator(Parameter rateParameter, Parameter bitParameter, double weight) {
        this.rateParameter = rateParameter;
        this.bitParameter = bitParameter;
        setWeight(weight);
    }

    public double doOperation() throws OperatorFailedException {
        int dim = rateParameter.getDimension() / 2;

        // Find a pair-set in which at least one bit is non-zero
        int index = -1;
        do {
            index = MathUtils.nextInt(dim);
        } while( bitParameter.getParameterValue(index) + bitParameter.getParameterValue(index+dim) < 1 );

        // Swap (bit,rate) values
        double tmpBit = bitParameter.getParameterValue(index);
        double tmpRate = rateParameter.getParameterValue(index);
        bitParameter.setParameterValue(index, bitParameter.getParameterValue(index+dim));
        rateParameter.setParameterValue(index, rateParameter.getParameterValue(index+dim));
        bitParameter.setParameterValue(index+dim,tmpBit);
        rateParameter.setParameterValue(index+dim,tmpRate);

        return 0;
    }

        // Interface MCMCOperator
    public final String getOperatorName() {
        return OPERATOR_NAME + "(" + bitParameter.getParameterName() +"," + rateParameter.getParameterName() + ")";
    }

    public final String getPerformanceSuggestion() {
        return "No performance suggestion";
    }

    public String toString() {
        return getOperatorName();
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return OPERATOR_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);

            Parameter ratesParameter = (Parameter) ((XMLObject)xo.getChild(RATES)).getChild(Parameter.class);
            Parameter bitsParameter = (Parameter) ((XMLObject)xo.getChild(BITS)).getChild(Parameter.class);


            return new RateBitExchangeOperator(ratesParameter,bitsParameter, weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a bit-flip operator on a given parameter.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(WEIGHT),
                new ElementRule(BITS, new XMLSyntaxRule[] {
            new ElementRule(Parameter.class)
            }),
                new ElementRule(RATES, new XMLSyntaxRule[] {
                        new ElementRule(Parameter.class)
                }),
        };

    };

    // Private instance variables
    private Parameter bitParameter = null;
    private Parameter rateParameter = null;
    
}
