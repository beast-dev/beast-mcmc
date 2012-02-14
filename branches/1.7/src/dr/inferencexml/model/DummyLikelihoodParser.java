package dr.inferencexml.model;

import dr.inference.model.DefaultModel;
import dr.inference.model.DummyLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Reads a distribution likelihood from a DOM Document element.
 */
public class DummyLikelihoodParser extends AbstractXMLObjectParser {

    public static final String DUMMY_LIKELIHOOD = "dummyLikelihood";

    public String getParserName() { return DUMMY_LIKELIHOOD; }

    public Object parseXMLObject(XMLObject xo) {

        Model model = (Model)xo.getChild(Model.class);
        Parameter parameter = (Parameter)xo.getChild(Parameter.class);

        if (model == null) {
            model = new DefaultModel();
        }
        final DummyLikelihood likelihood = new DummyLikelihood(model);

        ((DefaultModel)model).addVariable(parameter);

        return likelihood;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A function wraps a component model that would otherwise not be registered with the MCMC. Always returns a log likelihood of zero.";
    }

    public Class getReturnType() { return DummyLikelihood.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            new XORRule(
                    new ElementRule(Model.class, "A model element"),
                    new ElementRule(Parameter.class, "A parameter")
            )
    };
}
