package dr.evomodelxml.coalescent;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.coalescent.LinearGrowthModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a ConstantPopulation.
 */
public class LinearGrowthModelParser extends AbstractXMLObjectParser {

    public static String LINEAR_GROWTH_MODEL = "linearGrowth";
    public static String SLOPE = "slope";

    public String getParserName() {
        return LINEAR_GROWTH_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(SLOPE);
        Parameter slopeParam = (Parameter) cxo.getChild(Parameter.class);

        return new ConstantPopulationModel(slopeParam, units);
    }


    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A demographic model representing linear growth starting at time zero";
    }

    public Class getReturnType() {
        return LinearGrowthModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            XMLUnits.UNITS_RULE,
            new ElementRule(SLOPE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
    };

}
