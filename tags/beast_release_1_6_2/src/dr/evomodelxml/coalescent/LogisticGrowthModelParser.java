package dr.evomodelxml.coalescent;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.LogisticGrowthModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an XMLObject into LogisticGrowthModel.
 */
public class LogisticGrowthModelParser extends AbstractXMLObjectParser {

    public static String POPULATION_SIZE = "populationSize";
    public static String LOGISTIC_GROWTH_MODEL = "logisticGrowth";

    public static String GROWTH_RATE = "growthRate";
    public static String DOUBLING_TIME = "doublingTime";
    public static String TIME_50 = "t50";
    public static String ALPHA = "alpha";

    public String getParserName() {
        return LOGISTIC_GROWTH_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(POPULATION_SIZE);
        Parameter N0Param = (Parameter) cxo.getChild(Parameter.class);

        boolean usingGrowthRate = true;
        Parameter rParam;
        if (xo.getChild(GROWTH_RATE) != null) {
            cxo = xo.getChild(GROWTH_RATE);
            rParam = (Parameter) cxo.getChild(Parameter.class);
        } else {
            cxo = xo.getChild(DOUBLING_TIME);
            rParam = (Parameter) cxo.getChild(Parameter.class);
            usingGrowthRate = false;
        }

        cxo = xo.getChild(TIME_50);
        Parameter cParam = (Parameter) cxo.getChild(Parameter.class);

        return new LogisticGrowthModel(N0Param, rParam, cParam, 0.5, units, usingGrowthRate);
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
            new ElementRule(POPULATION_SIZE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                    "This parameter represents the carrying capacity (maximum population size). " +
                            "If the shape is very large then the current day population size will be very close to the carrying capacity."),
            new XORRule(

                    new ElementRule(GROWTH_RATE,
                            new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                            "This parameter determines the rate of growth during the exponential phase. See exponentialGrowth for details."),
                    new ElementRule(DOUBLING_TIME,
                            new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                            "This parameter determines the doubling time at peak growth rate.")
            ),
            new ElementRule(TIME_50,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                    "This parameter represents the time in the past when the population had half of the carrying capacity (population size). " +
                            "It is therefore a positive number with the same units as divergence times. " +
                            "A scale operator is recommended with a starting value near zero. " +
                            "A lower bound of zero should be employed and an upper bound is required!")
    };
}
