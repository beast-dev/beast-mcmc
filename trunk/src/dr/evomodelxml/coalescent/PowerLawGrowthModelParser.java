package dr.evomodelxml.coalescent;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.ExponentialGrowthModel;
import dr.evomodel.coalescent.PowerLawGrowthModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a ExponentialGrowth.
 */
public class PowerLawGrowthModelParser extends AbstractXMLObjectParser {

    public static String N0 = "n0";
    public static String POWER_LAW_GROWTH_MODEL = "powerLawGrowth";

    public static String POWER = "power";


    public String getParserName() {
        return POWER_LAW_GROWTH_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(N0);
        Parameter N0Param = (Parameter) cxo.getChild(Parameter.class);
        Parameter rParam;


        cxo = xo.getChild(POWER);
        rParam = (Parameter) cxo.getChild(Parameter.class);
        return new PowerLawGrowthModel(N0Param, rParam, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A demographic model of growth according to a power law.";
    }

    public Class getReturnType() {
        return PowerLawGrowthModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(N0,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),


            new ElementRule(POWER,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),

            XMLUnits.SYNTAX_RULES[0]
    };


}
