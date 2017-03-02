package dr.inferencexml.operators;

import dr.inference.model.Likelihood;
import dr.inference.model.PotentialDerivativeInterface;
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

    public final static String N_STEPS = "nSteps";
    public final static String STEP_SIZE = "stepSize";
    public final static String DRAW_VARIANCE = "drawVariance";

    @Override
    public String getParserName() {
        return HMC_OPERATOR;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        int nSteps = xo.getIntegerAttribute(N_STEPS);
        double stepSize = xo.getDoubleAttribute(STEP_SIZE);
        double drawVariance = xo.getDoubleAttribute(DRAW_VARIANCE);

        PotentialDerivativeInterface derivative = (PotentialDerivativeInterface) xo.getChild(PotentialDerivativeInterface.class);
        Parameter parameter = (Parameter) xo.getChild(Parameter.class);



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
            new ElementRule(PotentialDerivativeInterface.class),
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
