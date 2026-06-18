package dr.evomodelxml.epidemiology;

import dr.evomodel.epidemiology.StochasticSimulator;
import dr.xml.*;

public class StochasticTrajectorySimulatorParser extends AbstractXMLObjectParser{

    public static final String STOCHASTIC_TRAJECTORY_SIMULATOR = "stochasticTrajectorySimulator";

    public String getParserName() {
        return STOCHASTIC_TRAJECTORY_SIMULATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        StochasticSimulator stochasticSimulator = (StochasticSimulator) xo.getChild(StochasticSimulator.class);

        stochasticSimulator.simulateTrajectory();

        return stochasticSimulator;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element will simulate a stochastic trajectory using a stochastic simulator";
    }

    public Class getReturnType() {
        return StochasticSimulator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(StochasticSimulator.class),
    };
}