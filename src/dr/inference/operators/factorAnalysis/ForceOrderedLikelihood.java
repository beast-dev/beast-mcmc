package dr.inference.operators.factorAnalysis;

import dr.inference.model.*;
import dr.xml.*;

public class ForceOrderedLikelihood extends AbstractModelLikelihood {

    private final Parameter parameter;

    /**
     * @param name Model Name
     */
    public ForceOrderedLikelihood(String name, Parameter parameter) {
        super(name);
        this.parameter = parameter;
        addVariable(parameter);
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // do nothing
    }

    @Override
    protected void storeState() {
        // do nothing
    }

    @Override
    protected void restoreState() {
        // do nothing
    }

    @Override
    protected void acceptState() {
        // do nothing
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // TODO
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        int n = parameter.getDimension();
        for (int i = 1; i < n; i++) {
            if (parameter.getParameterValue(i) > parameter.getParameterValue(i - 1)) {
                return Double.NEGATIVE_INFINITY;
            }
        }
        return 0;
    }

    @Override
    public void makeDirty() {
        // do nothing
    }


    private static final String FORCE_ORDERED = "forceOrderedLikelihood";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            Parameter parameter = (Parameter) xo.getChild(Parameter.class);
            return new ForceOrderedLikelihood(FORCE_ORDERED, parameter);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            };
        }

        @Override
        public String getParserDescription() {
            return "Likelihood that forces parameters to be in descending order";
        }

        @Override
        public Class getReturnType() {
            return ForceOrderedLikelihood.class;
        }

        @Override
        public String getParserName() {
            return FORCE_ORDERED;
        }
    };
}
