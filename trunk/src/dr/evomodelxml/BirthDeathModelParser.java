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
    public static String SAMPLE_RATE = "sampleRate";

    public static String BIRTH_DEATH = "birthDeath";
    public static String BIRTHDIFF_RATE_PARAM_NAME = BIRTH_DEATH + ".BminusDRate";
    public static String RELATIVE_DEATH_RATE_PARAM_NAME = BIRTH_DEATH + ".DoverB";

    public String getParserName() {
        return BIRTH_DEATH_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLParser.Utils.getUnitsAttr(xo);

        Parameter birthParameter = (Parameter) xo.getElementFirstChild(BIRTHDIFF_RATE);
        Parameter deathParameter = (Parameter) xo.getElementFirstChild(RELATIVE_DEATH_RATE);
        Parameter sampleParameter = xo.hasChildNamed(SAMPLE_RATE) ?
                (Parameter) xo.getElementFirstChild(SAMPLE_RATE) : null;

        Logger.getLogger("dr.evomodel").info("Using Gernhard08 birth-death model on tree: Gernhard T (2008) J Theor Biol, In press");

        return new BirthDeathGernhard08Model(birthParameter, deathParameter, sampleParameter, units);
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

    private final XMLSyntaxRule[] rules = {
            new ElementRule(BIRTHDIFF_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(RELATIVE_DEATH_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(SAMPLE_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            XMLUnits.SYNTAX_RULES[0]
    };
}