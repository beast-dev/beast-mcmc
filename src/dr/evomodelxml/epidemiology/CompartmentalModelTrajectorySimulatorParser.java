package dr.evomodelxml.epidemiology;

import dr.evomodel.epidemiology.CompartmentalModelSimulator;
import dr.inference.model.Parameter;
import dr.xml.*;

public class CompartmentalModelTrajectorySimulatorParser extends AbstractXMLObjectParser{

    public static final String COMPARTMENTAL_MODEL_TRAJECTORY_SIMULATOR = "compartmentalModelTrajectorySimulator";

    public String getParserName() {
        return COMPARTMENTAL_MODEL_TRAJECTORY_SIMULATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        CompartmentalModelSimulator simulator = (CompartmentalModelSimulator) xo.getChild(CompartmentalModelSimulator.class);

        simulator.simulateTrajectory();

        // accept the simulated values so they are not restored to their
        // initial XML values during BEAST's post-parsing initialization
        for (Parameter p : simulator.compartmentalModel.compartmentCounts) {
            p.acceptParameterValues();
        }


        return simulator;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element will simulate a stochastic trajectory using a stochastic simulator";
    }

    public Class getReturnType() {
        return CompartmentalModelSimulator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(CompartmentalModelSimulator.class),
    };
}