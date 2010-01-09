package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * @author Chieh-Hsi Wu
 * Operator that performs independent random walks on k elements in a parameter
 */
public class MultipleRandomWalkIntegerOperator extends RandomWalkIntegerOperator{
    public static final String SAMPLE_SIZE = "sampleSize";
    private int sampleSize = 1;
    public MultipleRandomWalkIntegerOperator(Parameter parameter, int windowSize, int sampleSize, double weight) {
        super(parameter,windowSize, weight);
        this.sampleSize = sampleSize;
    }

    public double doOperation() {

        int[] shuffledInd = MathUtils.shuffled(parameter.getDimension());
        int index;
        for(int i = 0; i < sampleSize; i++){

            index = shuffledInd[i];

            // a random non zero integer around old value within windowSize * 2
            int oldValue = (int) parameter.getParameterValue(index);
            int newValue;
            int roll = MathUtils.nextInt(2 * windowSize);
            if (roll >= windowSize) {
                newValue = oldValue + 1 + roll - windowSize;

                if (newValue > parameter.getBounds().getUpperLimit(index))
                    newValue = 2 * (int)(double)parameter.getBounds().getUpperLimit(index) - newValue;
            } else {
                newValue = oldValue - 1 - roll;

                if (newValue < parameter.getBounds().getLowerLimit(index))
                    newValue = 2 * (int)(double)parameter.getBounds().getLowerLimit(index) - newValue;
            }

            parameter.setParameterValue(index, newValue);


        }




        return 0.0;
    }

    public static dr.xml.XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return "multipleRandomWalkIntegerOperator";
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);

            double w = xo.getDoubleAttribute(WINDOW_SIZE);
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
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(WINDOW_SIZE),
                AttributeRule.newDoubleRule(SAMPLE_SIZE),
                AttributeRule.newDoubleRule(WEIGHT),
                new ElementRule(Parameter.class)
        };

    };



}
