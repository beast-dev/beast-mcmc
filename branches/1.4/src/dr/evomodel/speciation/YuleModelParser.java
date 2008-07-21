package dr.evomodel.speciation;

import dr.evomodel.speciation.BirthDeathGernhard08Model;
import dr.evoxml.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Alexei Drummond
 */
public class YuleModelParser extends AbstractXMLObjectParser {

    public static final String YULE_MODEL = "yuleModel";
    public static String BIRTH_RATE = "birthRate";

    public String getParserName() {
        return YuleModelParser.YULE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        int units = XMLParser.Utils.getUnitsAttr(xo);

        XMLObject cxo = (XMLObject) xo.getChild(BIRTH_RATE);
        Parameter brParameter = (Parameter) cxo.getChild(Parameter.class);
        Parameter deathParameter = new Parameter.Default(0.0);

        return new BirthDeathGernhard08Model(brParameter, deathParameter, units);
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

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(BIRTH_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            XMLUnits.SYNTAX_RULES[0]
    };
}