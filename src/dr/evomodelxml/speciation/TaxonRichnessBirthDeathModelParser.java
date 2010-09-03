package dr.evomodelxml.speciation;

import dr.evolution.util.Units;
import dr.evomodel.speciation.TaxonRichnessBirthDeathModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class TaxonRichnessBirthDeathModelParser extends AbstractXMLObjectParser {

    public static final String TAXON_RICHNESS_MODEL = "randomLocalYuleModel";
    public static final String MEAN_RATE = "meanRate";
    public static final String BIRTH_RATE = "birthRates";
    public static final String BIRTH_RATE_INDICATORS = "indicators";
    public static final String RATES_AS_MULTIPLIERS = "ratesAsMultipliers";

    public String getParserName() {
        return TAXON_RICHNESS_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(BIRTH_RATE);
        Parameter brParameter = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(BIRTH_RATE_INDICATORS);
        Parameter indicatorsParameter = (Parameter) cxo.getChild(Parameter.class);

        Parameter meanRate = (Parameter) xo.getElementFirstChild(MEAN_RATE);

        boolean ratesAsMultipliers = xo.getBooleanAttribute(RATES_AS_MULTIPLIERS);

        int dp = xo.getAttribute("dp", 4);

        return new TaxonRichnessBirthDeathModel(brParameter, indicatorsParameter, meanRate, ratesAsMultipliers, units, dp);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A speciation model of a Yule process whose rate can change at random nodes in the tree.";
    }

    public Class getReturnType() {
        return TaxonRichnessBirthDeathModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(BIRTH_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(BIRTH_RATE_INDICATORS,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(MEAN_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            AttributeRule.newBooleanRule(RATES_AS_MULTIPLIERS),
            AttributeRule.newIntegerRule("dp", true),
            XMLUnits.SYNTAX_RULES[0]
    };
}
