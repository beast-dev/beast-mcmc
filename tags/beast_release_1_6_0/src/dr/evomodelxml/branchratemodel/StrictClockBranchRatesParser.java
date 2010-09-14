package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 */
public class StrictClockBranchRatesParser extends AbstractXMLObjectParser {

    public static final String STRICT_CLOCK_BRANCH_RATES = "strictClockBranchRates";
    public static final String RATE = "rate";

    public String getParserName() {
        return STRICT_CLOCK_BRANCH_RATES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter rateParameter = (Parameter) xo.getElementFirstChild(RATE);

        Logger.getLogger("dr.evomodel").info("Using strict molecular clock model.");

        return new StrictClockBranchRates(rateParameter);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return
                "This element provides a strict clock model. " +
                        "All branches have the same rate of molecular evolution.";
    }

    public Class getReturnType() {
        return StrictClockBranchRates.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(RATE, Parameter.class, "The molecular evolutionary rate parameter", false),
    };
}
