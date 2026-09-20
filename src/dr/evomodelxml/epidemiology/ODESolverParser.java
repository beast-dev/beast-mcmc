package dr.evomodelxml.epidemiology;

import dr.evomodel.epidemiology.ODESolver;
import dr.evomodel.epidemiology.TwoPathogenModel;
import dr.xml.*;

public class ODESolverParser extends AbstractXMLObjectParser {

    public static final String ODE_SOLVER = "odeSolver";
    public static final String ABSOLUTE_TOLERANCE = "absoluteTolerance";
    public static final String RELATIVE_TOLERANCE = "relativeTolerance";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TwoPathogenModel twoPathogenModel =
                (TwoPathogenModel) xo.getChild(TwoPathogenModel.class);

        if (xo.hasAttribute(ABSOLUTE_TOLERANCE) && xo.hasAttribute(RELATIVE_TOLERANCE)) {
            double absoluteTolerance = xo.getDoubleAttribute(ABSOLUTE_TOLERANCE);
            double relativeTolerance = xo.getDoubleAttribute(RELATIVE_TOLERANCE);
            return new ODESolver(twoPathogenModel, absoluteTolerance, relativeTolerance);
        }

        return new ODESolver(twoPathogenModel);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TwoPathogenModel.class),
            AttributeRule.newDoubleRule(ABSOLUTE_TOLERANCE, true),
            AttributeRule.newDoubleRule(RELATIVE_TOLERANCE, true),
    };

    @Override
    public String getParserDescription() {
        return "ODE solver for the two-pathogen compartmental model. Solves the " +
                "deterministic skeleton of the model using an adaptive Dormand-Prince " +
                "RK8(5,3) integrator.";
    }

    @Override
    public Class getReturnType() {
        return ODESolver.class;
    }

    @Override
    public String getParserName() {
        return ODE_SOLVER;
    }

}