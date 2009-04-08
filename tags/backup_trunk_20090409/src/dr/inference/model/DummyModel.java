package dr.inference.model;

import dr.inference.loggers.LogColumn;
import dr.xml.*;

/**
 * @author Marc Suchard
 */
public class DummyModel extends AbstractModel implements Likelihood {

    public static final String DUMMY_MODEL = "dummyModel";

    public DummyModel() {
        super(DUMMY_MODEL);
    }

    public DummyModel(Parameter parameter) {
        super(DUMMY_MODEL);
        addParameter(parameter);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    protected final void handleParameterChangedEvent(Parameter parameter, int index, Parameter.ChangeType type) {

    }

    protected void storeState() {

    }

    protected void restoreState() {

    }

    protected void acceptState() {

    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        return 0;
    }

    public void makeDirty() {

    }

    public LogColumn[] getColumns() {
        return new LogColumn[0];
    }

    /**
     * Reads a distribution likelihood from a DOM Document element.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return DUMMY_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) {

//           Parameter parameter = (Parameter) xo.getChild(Parameter.class);
            DummyModel likelihood = new DummyModel();

            for (int i = 0; i < xo.getChildCount(); i++) {
                Parameter parameter = (Parameter) xo.getChild(i);
                likelihood.addParameter(parameter);
            }

            return likelihood;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A function wraps a component model that would otherwise not be registered with the MCMC. Always returns a log likelihood of zero.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(Parameter.class, 1, Integer.MAX_VALUE)
        };

    };


}


