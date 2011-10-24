package dr.evomodel.epidemiology;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.ExponentialGrowthModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a SIR model.
 */
public class SIRModelParser extends AbstractXMLObjectParser {

    public static String SIR_MODEL = "sirModel";

    public static String TRANSMISSION_RATE = "transmissionRate";
    public static String RECOVERY_RATE = "recoveryRate";


    public String getParserName() {
        return SIR_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(TRANSMISSION_RATE);
        Parameter transmissionRateParameter = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(RECOVERY_RATE);
        Parameter recoveryRateParameter = (Parameter) cxo.getChild(Parameter.class);

        return new SIRModel(transmissionRateParameter, recoveryRateParameter, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A demographic model of exponential growth.";
    }

    public Class getReturnType() {
        return ExponentialGrowthModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(TRANSMISSION_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(RECOVERY_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            XMLUnits.SYNTAX_RULES[0]
    };


}
