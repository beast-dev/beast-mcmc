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
    private int numSteps;
    private double rate;

    private double[] startValues;
// todo: cleanup and implement or erase numSteps option

    public TemperOperator(Parameter parameter,
                          Parameter target,
                          double rate,
                          int numSteps,
                          double weight) {

        this.parameter = parameter;
        this.target = target;
        this.rate = rate;
        this.numSteps = numSteps;

        startValues = parameter.getParameterValues();

        if (numSteps !=0) {
            throw new RuntimeException("Number of steps not yet implemented, please specifiy a rate.");
        }

//        if(numSteps != 0 && rate == 0){
//            if (parameter.getDimension() > 1) {
//                throw new RuntimeException("Tempering not yet implemented with for parameters with dimension > 1");
//            }
//            this.rate = setRate();
//        }
//        else if (rate !=0 && numSteps == 0){
//            this.numSteps = setSteps();
//        }
//        else{
//            throw new RuntimeException("You must provide exactly one: rate or steps.");
//        }

        setWeight(weight);
    }

    @Override
    public String getOperatorName() {
        return BAYESIAN_BRIDGE_PARSER;
    }

//    private void setRate(){
//        double start = this.parameter.getValue(0);
//        double end = this.target.getValue(0);
//        this.rate = Math.log(end / start) / this.numSteps;
//    }
//
//    private void setSteps(){
//        double start = this.parameter.getValue(0);
//        double end = this.target.getValue(0);
//        this.numSteps = (int) Math.round(Math.log(end / start) / this.rate);
//    }

    @Override
    public double doOperation() {
        // plan: start + (end - start) * percentComplete
        // where percentComplete = [ exp(-Nsteps * exp(-rate)) - 1 ]  is in (0, 1), and importantly goes from *0* to *1*
        // Note we achieve percentComplete = 1 only when Nsteps * exp(-rate) := 'objective' is large

        // factor out a -1 to read code below:
        // start + (start - end) * -1*percentComplete

//        if (startValues == null) {
//            startValues = parameter.getParameterValues();
//        }

        double currentCount = getCount() + 1;

        // seems unnecessary, could just input exp(-rate) directly as "rate" in XML if desired.
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

            // if less than (approx) 90% of the way there, keep going: exp(-2.3) ~ 0.1
            if (objective < 2.3) {
                double startValue = startValues[i];
                double x = startValue + (startValue - targetValue) * (Math.exp(-objective) - 1);
                parameter.setParameterValue(i, x);

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

            if(rate < 0){
                throw new XMLParseException("Rate cannot be negative");
            }

            else if (numSteps == 0 && rate == 0) {
                throw new XMLParseException("Rate OR number of steps must be provided.");
            } else if (numSteps != 0 && rate != 0) {
                throw new XMLParseException("Cannot specify both the tempering rate and the number of steps.");
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
