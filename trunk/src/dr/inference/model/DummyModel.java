package dr.inference.model;

import dr.xml.*;
import dr.inference.loggers.LogColumn;

/**
 * @author Marc Suchard
 */
public class DummyModel extends AbstractModel implements Likelihood {

	public static final String DUMMY_MODEL = "dummyModel";

	public DummyModel(Parameter parameter) {
		super(DUMMY_MODEL);
        addParameter(parameter);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    protected void handleParameterChangedEvent(Parameter parameter, int index) {

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

		public String getParserName() { return DUMMY_MODEL; }

		public Object parseXMLObject(XMLObject xo) {

           Parameter parameter = (Parameter) xo.getChild(Parameter.class);
			DummyModel likelihood = new DummyModel(parameter);

			return likelihood;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "A function wraps a component model that would otherwise not be registered with the MCMC. Always returns a log likelihood of zero.";
		}

		public Class getReturnType() { return Likelihood.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new ElementRule(Parameter.class)
		};

	};



}


