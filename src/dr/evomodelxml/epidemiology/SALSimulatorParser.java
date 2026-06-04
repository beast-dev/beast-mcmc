package dr.evomodelxml.epidemiology;

import dr.evomodel.epidemiology.CompartmentalModel;
import dr.evomodel.epidemiology.SALSimulator;
import dr.evomodel.epidemiology.TauLeapingSimulator;
import dr.xml.*;

public class SALSimulatorParser extends AbstractXMLObjectParser {

    public static final String SAL_SIMULATOR = "salSimulator";
    public static final String EPSILON = "epsilon";
    public static final String CRITICAL_NUMBER = "criticalNumber";

    public String getParserName() {
        return SAL_SIMULATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        // default value in Cao et al. step-size selection algorithm is 0.03
        final double epsilon = xo.getDoubleAttribute(EPSILON, 0.03);
        // a reaction is critical if there are fewer than criticalNumber left in that compartment
        // critical number usually between 2 and 20
        final int criticalNumber = xo.getIntegerAttribute(CRITICAL_NUMBER, 5);

        CompartmentalModel compartmentalModel = (CompartmentalModel) xo.getChild(CompartmentalModel.class);

        SALSimulator simulator = new SALSimulator(compartmentalModel, epsilon, criticalNumber);

        return simulator;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents a stochastic simulator that uses the step-anticipation-leap (SAL) generalization " +
                "(Sehl et al., 2009) of the tau leaping procedure with the step size selection algorithm as outlined " +
                "by Cao et al. (2006)";
    }

    public Class getReturnType() {
        return SALSimulator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(CompartmentalModel.class),
            AttributeRule.newDoubleRule(EPSILON, true),
            AttributeRule.newIntegerRule(CRITICAL_NUMBER, true),
    };

}

