package dr.evomodelxml.epidemiology;

import dr.evomodel.epidemiology.CompartmentalModel;
import dr.evomodel.epidemiology.ExactSimulator;
import dr.xml.*;

public class ExactSimulatorParser extends AbstractXMLObjectParser {

    public static final String EXACT_SIMULATOR = "exactStochasticSimulator";

    public String getParserName() {
        return EXACT_SIMULATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        CompartmentalModel compartmentalModel = (CompartmentalModel) xo.getChild(CompartmentalModel.class);

        ExactSimulator simulator = new ExactSimulator(compartmentalModel);

        return simulator;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an exact stochastic simulator that uses the Gillespie algorithm";
    }

    public Class getReturnType() {
        return ExactSimulator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(CompartmentalModel.class),
    };

}

