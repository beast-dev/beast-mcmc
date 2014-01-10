package dr.evomodel.epidemiology.casetocase;

import dr.evolution.util.Units;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author mhall
 */
public class LogisticGrowthN0N50ModelParser extends AbstractXMLObjectParser {

    public static String LOGISTIC_GROWTH_MODEL = "logisticGrowthN0N50";

    public static String STARTING_POPULATION_SIZE = "startingPopulationSize";
    public static String EXCESS_CARRYING_CAPACITY = "excessCarryingCapacity";
    public static String TIME_OF_ORIGIN = "timeOfOrigin";
    public static String GROWTH_RATE = "growthRate";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(STARTING_POPULATION_SIZE);
        Parameter n50Param = (Parameter)cxo.getChild(Parameter.class);

        cxo = xo.getChild(EXCESS_CARRYING_CAPACITY);
        Parameter n0Minus2N50LimitParam = (Parameter)cxo.getChild(Parameter.class);

        cxo = xo.getChild(TIME_OF_ORIGIN);
        Parameter t50Param = (Parameter)cxo.getChild(Parameter.class);

        cxo = xo.getChild(GROWTH_RATE);
        Parameter rParam = (Parameter)cxo.getChild(Parameter.class);

        return new LogisticGrowthN0N50Model(n50Param, n0Minus2N50LimitParam, rParam, t50Param, units);
    }

    public String getParserDescription() {
        return "Logistic growth demographic model";
    }

    public Class getReturnType() {
        return LogisticGrowthN0N50Model.class;
    }

    public String getParserName() {
        return LOGISTIC_GROWTH_MODEL;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            XMLUnits.SYNTAX_RULES[0],
            new ElementRule(STARTING_POPULATION_SIZE,
                    new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
                    "This parameter represents the population size at the start of the process"),
            new ElementRule(GROWTH_RATE,
                    new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
                    "This parameter determines the rate of growth during the exponential phase."),
            new ElementRule(TIME_OF_ORIGIN,
                    new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
                    "This parameter represents the time at which the process started"),
            new ElementRule(EXCESS_CARRYING_CAPACITY,
                    new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
                    "The minimum carrying capacity in this model is twice the starting population size. This " +
                            "parameter represents how much larger the carrying capacity is than this"),
    };

}
