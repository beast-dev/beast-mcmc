package dr.evomodel.epidemiology;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.ExponentialGrowthModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a SIR model.
 *
 * @author Trevor Bedford
 * @author Andrew Rambaut
 */
public class SIRModelParser extends AbstractXMLObjectParser {

    public static String SIR_MODEL = "sirEpidemiology";
    public static String REPRODUCTIVE_NUMBER = "reproductiveNumber";
    public static String RECOVERY_RATE = "recoveryRate";
    public static String HOST_POPULATION_SIZE = "hostPopulationSize";
    public static String PROPORTIONS = "proportions";

    public String getParserName() {
        return SIR_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(REPRODUCTIVE_NUMBER);
        Parameter reproductiveNumberParameter = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(RECOVERY_RATE);
        Parameter recoveryRateParameter = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(HOST_POPULATION_SIZE);
        Parameter hostPopulationSizeParameter = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(PROPORTIONS);
        Parameter proportionsParameter = (Parameter) cxo.getChild(Parameter.class);

        return new SIRModel(reproductiveNumberParameter, recoveryRateParameter,
                hostPopulationSizeParameter, proportionsParameter, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A demographic model of epidemic spread.";
    }

    public Class getReturnType() {
        return SIRModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(REPRODUCTIVE_NUMBER,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(RECOVERY_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(HOST_POPULATION_SIZE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(PROPORTIONS,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            XMLUnits.SYNTAX_RULES[0]
    };


}