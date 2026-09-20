package dr.evomodelxml.epidemiology;

import dr.evomodel.epidemiology.CompartmentalModel;
import dr.evomodel.epidemiology.ExactStochasticSimulator;
import dr.inference.model.Parameter;
import dr.xml.*;

public class ExactStochasticSimulatorParser extends AbstractXMLObjectParser {

    public static final String EXACT_SIMULATOR = "exactStochasticSimulator";

    public String getParserName() {
        return EXACT_SIMULATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        CompartmentalModel compartmentalModel = (CompartmentalModel) xo.getChild(CompartmentalModel.class);

        Parameter elapsedTimeOneTrajectory = (Parameter) xo.getChild(Parameter.class);

        ExactStochasticSimulator simulator;

        if(elapsedTimeOneTrajectory == null){
            simulator = new ExactStochasticSimulator(compartmentalModel);
        }else {
            simulator = new ExactStochasticSimulator(compartmentalModel, elapsedTimeOneTrajectory);
        }

        return simulator;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an exact stochastic simulator that uses the Gillespie algorithm";
    }

    public Class getReturnType() {
        return ExactStochasticSimulator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(CompartmentalModel.class),
            new ElementRule(Parameter.class, true),
    };

}

