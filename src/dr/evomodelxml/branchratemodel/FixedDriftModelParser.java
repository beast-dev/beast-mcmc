package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.FixedDriftModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * Created by mandevgill on 12/22/14.
 */
public class FixedDriftModelParser extends AbstractXMLObjectParser {
    public static final String FIXED_DRIFT = "fixedDriftModel";
    public static final String RATE_ONE = "rateOne";
    public static final String RATE_TWO = "rateTwo";
    public static final String REMAINING_RATES = "remainingRates";
    public static final String RATE_ONE_ID = "rateOneID";
    public static final String RATE_TWO_ID = "rateTwoID";

    public String getParserName() {
        return FIXED_DRIFT;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String idOne = xo.getStringAttribute(RATE_ONE_ID);

        String idTwo = xo.getStringAttribute(RATE_TWO_ID);

        Parameter rateOne = (Parameter) xo.getElementFirstChild(RATE_ONE);

        Parameter rateTwo = (Parameter) xo.getElementFirstChild(RATE_TWO);

        Parameter remainingRates = (Parameter) xo.getElementFirstChild(REMAINING_RATES);

        // Parameter rateOneNum = (Parameter) xo.getElementFirstChild(RATE_ONE_NUM);

        // Parameter rateTwoNum = (Parameter) xo.getElementFirstChild(RATE_TWO_NUM);


        Logger.getLogger("dr.evomodel").info("Using fixed drift model.");


        return new FixedDriftModel(rateOne, rateTwo, remainingRates,
                idOne, idTwo);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return
                "This element returns a relaxed drift model.";
    }

    public Class getReturnType() {
        return FixedDriftModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newStringRule(RATE_ONE_ID, false),
            AttributeRule.newStringRule(RATE_TWO_ID, false),
            new ElementRule(RATE_ONE, Parameter.class, "rate one parameter", false),
            new ElementRule(RATE_TWO, Parameter.class, "rate two parameter", false),
            new ElementRule(REMAINING_RATES, Parameter.class, "remaining rates parameter", false)
            // new ElementRule(RATE_ONE_NUM, Parameter.class, "rate one node num", false),
            // new ElementRule(RATE_TWO_NUM, Parameter.class, "rate two node num", false)
    };


}
