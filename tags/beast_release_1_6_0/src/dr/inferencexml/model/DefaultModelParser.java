package dr.inferencexml.model;

import dr.inference.model.DefaultModel;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLSyntaxRule;

/**
 *
 */
public class DefaultModelParser extends AbstractXMLObjectParser {

    public static final String DUMMY_MODEL = "dummyModel";

    public String getParserName() {
        return DUMMY_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) {

        DefaultModel likelihood = new DefaultModel();

        for (int i = 0; i < xo.getChildCount(); i++) {
            Parameter parameter = (Parameter) xo.getChild(i);
            likelihood.addVariable(parameter);
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
        return DefaultModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Parameter.class, 1, Integer.MAX_VALUE)
    };

}
