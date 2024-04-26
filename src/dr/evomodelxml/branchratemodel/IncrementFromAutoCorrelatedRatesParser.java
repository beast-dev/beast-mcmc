package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.AutoCorrelatedGradientWrtIncrements;
import dr.inference.model.Parameter;
import dr.xml.*;

public class IncrementFromAutoCorrelatedRatesParser extends AbstractXMLObjectParser {

    private static final String INCREMENTS = "increments";

    public String getParserName() {
        return INCREMENTS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        AutoCorrelatedGradientWrtIncrements gradient = (AutoCorrelatedGradientWrtIncrements) xo.getChild(AutoCorrelatedGradientWrtIncrements.class);

        return gradient.getParameter();
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns increments of an auto-correlated branch rate model as a proxy parameter";
    }

    public Class getReturnType() {
        return Parameter.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(AutoCorrelatedGradientWrtIncrements.class),
    };
}

