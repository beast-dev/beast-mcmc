package dr.evomodelxml.coalescent;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.CataclysmicDemographicModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a ExponentialGrowth.
 */
public class CataclysmicDemographicModelParser extends AbstractXMLObjectParser {

    public static final String POPULATION_SIZE = "populationSize";
    public static final String GROWTH_RATE = "growthRate";
    public static final String SPIKE_SIZE = "spikeFactor";
    public static final String TIME_OF_CATACLYSM = "timeOfCataclysm";

    public static final String CATACLYSM_MODEL = "cataclysm";

    public String getParserName() {
        return CATACLYSM_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(POPULATION_SIZE);
        Parameter N0Param = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(GROWTH_RATE);
        Parameter rParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(SPIKE_SIZE);
        Parameter N1Param = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(TIME_OF_CATACLYSM);
        Parameter tParam = (Parameter) cxo.getChild(Parameter.class);

        return new CataclysmicDemographicModel(N0Param, N1Param, rParam, tParam, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A demographic model of exponential growth.";
    }

    public Class getReturnType() {
        return CataclysmicDemographicModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(POPULATION_SIZE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(GROWTH_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                    "The rate of exponential growth before the cataclysmic event."),
            new ElementRule(SPIKE_SIZE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                    "The factor larger the population size was at its height."),
            new ElementRule(TIME_OF_CATACLYSM,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                    "The time of the cataclysmic event that lead to exponential decline."),
            XMLUnits.SYNTAX_RULES[0]
    };
}
