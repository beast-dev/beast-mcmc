package dr.evomodelxml;

import dr.evolution.util.Units;
import dr.evomodel.speciation.BirthDeathGernhard08Model;
import dr.evoxml.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Joseph Heled
 */
public class BirthDeathModelParser extends AbstractXMLObjectParser {


    public static final String BIRTH_DEATH_MODEL = "birthDeathModel";
    public static String BIRTHDIFF_RATE = "birthMinusDeathRate";
    public static String RELATIVE_DEATH_RATE = "relativeDeathRate";

    public static String BIRTHDIFF_RATE_PARAM_NAME = "birthDeath.BminusDRate";
    public static String RELATIVE_DEATH_RATE_PARAM_NAME = "birthDeath.DoverB";

    public String getParserName() {
        return BIRTH_DEATH_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLParser.Utils.getUnitsAttr(xo);

        Parameter birthParameter = (Parameter) xo.getElementFirstChild(BIRTHDIFF_RATE);
        Parameter deathParameter = (Parameter) xo.getElementFirstChild(RELATIVE_DEATH_RATE);

        Logger.getLogger("dr.evomodel").info("Using Gernhard 2008 birth-death model on tree: Gernhard T (2008) J Theor Biol, In press");

        return new BirthDeathGernhard08Model(birthParameter, deathParameter, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Gernhard (2008) model of speciation (equation at bottom of page 19 of draft).";
    }

    public Class getReturnType() {
        return BirthDeathGernhard08Model.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(BIRTHDIFF_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(RELATIVE_DEATH_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            XMLUnits.SYNTAX_RULES[0]
    };
}