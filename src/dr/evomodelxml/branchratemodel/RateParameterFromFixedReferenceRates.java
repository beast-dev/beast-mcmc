package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.FixedReferenceRates;
import dr.inference.model.Parameter;
import dr.xml.*;

public class RateParameterFromFixedReferenceRates extends AbstractXMLObjectParser {

    private static final String RATE_PARAMETER = "rateParameter";

    public String getParserName() {
        return RATE_PARAMETER;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        FixedReferenceRates fixedReferenceRates = (FixedReferenceRates) xo.getChild(FixedReferenceRates.class);

        return fixedReferenceRates.getParameter();
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns the rate parameter of a fixedReferenceRate branch rate wrapper";
    }

    public Class getReturnType() {
        return Parameter.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(FixedReferenceRates.class),
    };
}

