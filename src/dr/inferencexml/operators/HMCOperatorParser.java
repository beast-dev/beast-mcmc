package dr.inferencexml.operators;

import dr.inference.model.GradientWrtParameterProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.HMCOperator;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 * @author Max Tolkoff
 */
public class HMCOperatorParser extends AbstractXMLObjectParser {
    public final static String HMC_OPERATOR = "HMCOperator";
    public static final String HMC_OPERATOR2 = "hamiltonianMonteCarloOperator";

    public final static String N_STEPS = "nSteps";
    public final static String STEP_SIZE = "stepSize";
    public final static String DRAW_VARIANCE = "drawVariance";

    @Override
    public String getParserName() {
        return HMC_OPERATOR;
    }

    @Override
    public String[] getParserNames() {
        return new String[] { HMC_OPERATOR, HMC_OPERATOR2 };
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        int nSteps = xo.getIntegerAttribute(N_STEPS);
        double stepSize = xo.getDoubleAttribute(STEP_SIZE);
        double drawVariance = xo.getDoubleAttribute(DRAW_VARIANCE);

        GradientWrtParameterProvider derivative = (GradientWrtParameterProvider) xo.getChild(GradientWrtParameterProvider.class);
        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        if (derivative.getDimension() != parameter.getDimension()) {
            throw new XMLParseException("Gradient (" + derivative.getDimension() +
                    ") must be the same dimensions as the parameter (" + parameter.getDimension() + ")");
        }



        return new HMCOperator(CoercionMode.DEFAULT, weight, derivative, parameter,stepSize, nSteps, drawVariance);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newIntegerRule(N_STEPS),
            AttributeRule.newDoubleRule(STEP_SIZE),
            AttributeRule.newDoubleRule(DRAW_VARIANCE),
            new ElementRule(Parameter.class),
            new ElementRule(GradientWrtParameterProvider.class),
    };

    @Override
    public String getParserDescription() {
        return "Returns a Hamiltonian Monte Carlo transition kernel";
    }

    @Override
    public Class getReturnType() {
        return HMCOperator.class;
    }
}
