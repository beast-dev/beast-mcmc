package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.RandomLocalClockModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 */
public class RandomLocalClockModelParser extends AbstractXMLObjectParser {

    public static final String LOCAL_BRANCH_RATES = "randomLocalClockModel";
    public static final String RATE_INDICATORS = "rateIndicator";
    public static final String RATES = "rates";
    public static final String CLOCK_RATE = "clockRate";
    public static final String RATES_ARE_MULTIPLIERS = "ratesAreMultipliers";

    public String getParserName() {
        return LOCAL_BRANCH_RATES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        Parameter rateIndicatorParameter = (Parameter) xo.getElementFirstChild(RATE_INDICATORS);
        Parameter ratesParameter = (Parameter) xo.getElementFirstChild(RATES);
        Parameter meanRateParameter = null;

        if (xo.hasChildNamed(CLOCK_RATE)) {
            meanRateParameter = (Parameter) xo.getElementFirstChild(CLOCK_RATE);
        }

        boolean ratesAreMultipliers = xo.getAttribute(RATES_ARE_MULTIPLIERS, false);

        Logger.getLogger("dr.evomodel").info("Using random local clock (RLC) model.");
        Logger.getLogger("dr.evomodel").info("  rates at change points are parameterized to be " +
                (ratesAreMultipliers ? " relative to parent rates." : "independent of parent rates."));

        return new RandomLocalClockModel(tree, meanRateParameter, rateIndicatorParameter,
                ratesParameter, ratesAreMultipliers);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return
                "This element returns an random local clock (RLC) model." +
                        "Each branch either has a new rate or " +
                        "inherits the rate of the branch above it depending on the indicator vector, " +
                        "which is itself sampled.";
    }

    public Class getReturnType() {
        return RandomLocalClockModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(TreeModel.class),
            new ElementRule(RATE_INDICATORS, Parameter.class, "The rate change indicators parameter", false),
            new ElementRule(RATES, Parameter.class, "The rates parameter", false),
            new ElementRule(CLOCK_RATE, Parameter.class, "The mean rate across all local clocks", true),
            AttributeRule.newBooleanRule(RATES_ARE_MULTIPLIERS, false)
    };
}
