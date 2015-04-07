package dr.evomodelxml.coalescent;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.AsymptoticGrowthModel;
import dr.evomodel.coalescent.LogisticGrowthModel;
import dr.evomodel.coalescent.PeakAndDeclineModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an XMLObject into LogisticGrowthModel.
 */
public class AsymptoticGrowthModelParser extends AbstractXMLObjectParser {

    public static String ASYMPTOTE_VALUE = "asymptoteValue";
    public static String ASYMPTOTIC_GROWTH_MODEL = "asymptoticGrowth";

    public static String SHAPE = "shape";

    public String getParserName() {
        return ASYMPTOTIC_GROWTH_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(ASYMPTOTE_VALUE);
        Parameter asymptoticValueParam = (Parameter) cxo.getChild(Parameter.class);


        Parameter rParam;

        cxo = xo.getChild(SHAPE);
        Parameter shapeParam = (Parameter) cxo.getChild(Parameter.class);


        return new AsymptoticGrowthModel(asymptoticValueParam, shapeParam, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Logistic growth demographic model.";
    }

    public Class getReturnType() {
        return LogisticGrowthModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            XMLUnits.SYNTAX_RULES[0],
            new ElementRule(ASYMPTOTE_VALUE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),


                    new ElementRule(SHAPE,
                            new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
    };
}
