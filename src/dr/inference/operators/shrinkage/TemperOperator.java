package dr.inference.operators.shrinkage;

import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.xml.*;

import static dr.inferencexml.operators.shrinkage.BayesianBridgeShrinkageOperatorParser.BAYESIAN_BRIDGE_PARSER;

/**
 * @author Marc A. Suchard
 * @author Alexander Fisher
 */

public class TemperOperator extends SimpleMCMCOperator implements GibbsOperator {

    public static final String TEMPER_OPERATOR = "temperOperator";
    public static final String TARGET_PARAMETER = "targetParameter";
    public static final String NUM_STEPS = "numSteps";
    public static final String RATE = "rate";

    private final Parameter parameter;
    private final Parameter target;
    private final int numSteps = 0;
    private final double rate;

    private double[] startValues;

    public TemperOperator(Parameter parameter,
                          Parameter target,
                          double rate,
                          int numSteps,
                          double weight) {

        this.parameter = parameter;
        this.target = target;
        this.rate = rate;

        setWeight(weight);
    }

    @Override
    public String getOperatorName() {
        return BAYESIAN_BRIDGE_PARSER;
    }

    @Override
    public double doOperation() {

        if (startValues == null) {
            startValues = parameter.getParameterValues();
        }

        double currentCount = getCount() + 1;

        double objective = currentCount * Math.exp(-rate);

//        for (int i = 0; i < parameter.getDimension(); ++i) {
//
//            double currentValue = parameter.getParameterValue(i);
//            double targetValue = target.getParameterValue(i % target.getDimension());
//
//            if (objective > 0.1) {
//                double startValue = startValues[i];
//                double x = startValue + (startValue - targetValue) * Math.exp(-objective);
//            } else if (currentValue != targetValue) {
//                parameter.setParameterValue(i, targetValue);
//            }
//        }

        for (int i = 0; i < parameter.getDimension(); ++i) {

            double currentValue = parameter.getParameterValue(i);
            double targetValue = target.getParameterValue(i % target.getDimension());

            if (objective > 0.1) {
                double startValue = startValues[i];
                double x = startValue + (startValue - targetValue) * Math.exp(-objective);
            } else if (currentValue != targetValue) {
                parameter.setParameterValue(i, targetValue);
            }
        }

        return 0;
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TEMPER_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(SimpleMCMCOperator.WEIGHT);

            Parameter hotParameter = (Parameter) xo.getChild(Parameter.class);

            XMLObject cxo = xo.getChild(TARGET_PARAMETER);
            Parameter targetParameter = (Parameter) cxo.getChild(Parameter.class);

            if (targetParameter.getDimension() > hotParameter.getDimension()) {
                throw new XMLParseException("Target parameter cannot have more dimensions than the tempered parameter.");
            }

            int numSteps = xo.getAttribute(NUM_STEPS, 0);
            double rate = xo.getAttribute(RATE, 0);

            if (numSteps == 0 && rate == 0) {
                throw new XMLParseException("Rate OR number of steps must be provided.");
            } else if (numSteps != 0 && rate != 0) {
                throw new XMLParseException("Cannot specify both tempering rate and the number of steps.");
            }

            TemperOperator temper = new TemperOperator(hotParameter, targetParameter, rate, numSteps, weight);
            return temper;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        @Override
        public String getParserDescription() {
            return "Tempers (anneals) a parameter from starting value to target value.";
        }

        @Override
        public Class getReturnType() {
            return TemperOperator.class;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(Parameter.class),
                new ElementRule(TARGET_PARAMETER, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class),
                }),
                AttributeRule.newStringRule(NUM_STEPS, true),
                AttributeRule.newStringRule(RATE, true),
        };
    };
}
