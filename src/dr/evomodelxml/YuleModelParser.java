package dr.evomodelxml;

import dr.evolution.util.Units;
import dr.evomodel.speciation.BirthDeathGernhard08Model;
import dr.evomodel.speciation.YuleModel;
import dr.evoxml.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 */
public class YuleModelParser extends AbstractXMLObjectParser {

    public static String BIRTH_RATE = "birthRate";

    public String getParserName() {
        return YuleModel.YULE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLParser.Utils.getUnitsAttr(xo);

        XMLObject cxo = (XMLObject) xo.getChild(BIRTH_RATE);
        Parameter brParameter = (Parameter) cxo.getChild(Parameter.class);
        Parameter deathParameter = new Parameter.Default(0.0);

        Logger.getLogger("dr.evomodel").info("Using Yule prior on tree");

        return new BirthDeathGernhard08Model(brParameter, deathParameter, null, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A speciation model of a simple constant rate Birth-death process.";
    }

    public Class getReturnType() {
        return BirthDeathGernhard08Model.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(BIRTH_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            XMLUnits.SYNTAX_RULES[0]
    };
}