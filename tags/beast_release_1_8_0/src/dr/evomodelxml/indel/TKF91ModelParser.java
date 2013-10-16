package dr.evomodelxml.indel;

import dr.evolution.util.Units;
import dr.evomodel.indel.TKF91Model;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a ConstantPopulation.
 */
public class TKF91ModelParser extends AbstractXMLObjectParser {

    public static final String TKF91_MODEL = "tkf91Model";
    public static final String TKF91_LENGTH_DIST = "lengthDistribution";
    public static final String BIRTH_RATE = "birthRate";
    public static final String DEATH_RATE = "deathRate";

    public String getParserName() {
        return TKF91_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter lengthDistParameter = (Parameter) xo.getElementFirstChild(TKF91_LENGTH_DIST);
        Parameter deathParameter = (Parameter) xo.getElementFirstChild(DEATH_RATE);
        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        return new TKF91Model(lengthDistParameter, deathParameter, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "The TKF91 (Thorne, Kishino & Felsenstein 1991) model of insertion-deletion.";
    }

    public Class getReturnType() {
        return TKF91Model.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TKF91_LENGTH_DIST,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(DEATH_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            XMLUnits.SYNTAX_RULES[0]
    };
}
