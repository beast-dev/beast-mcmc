package dr.evomodelxml.coalescent;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.ExponentialSawtoothModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a ExponentialGrowth.
 */
public class ExponentialSawtoothModelParser extends AbstractXMLObjectParser {

    public static final String POPULATION_SIZE = "populationSize";
    public static final String GROWTH_RATE = "growthRate";
    public static final String WAVELENGTH = "wavelength";
    public static final String OFFSET = "offset";

    public static final String EXPONENTIAL_SAWTOOTH = "exponentialSawtooth";

    public String getParserName() {
        return EXPONENTIAL_SAWTOOTH;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(POPULATION_SIZE);
        Parameter N0Param = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(GROWTH_RATE);
        Parameter rParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(WAVELENGTH);
        Parameter wavelengthParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(OFFSET);
        Parameter tParam = (Parameter) cxo.getChild(Parameter.class);

        return new ExponentialSawtoothModel(N0Param, rParam, wavelengthParam, tParam, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A demographic model of succesive exponential growth and periodic population crashes.";
    }

    public Class getReturnType() {
        return ExponentialSawtoothModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(POPULATION_SIZE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(GROWTH_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                    "The rate of exponential growth during the growth phase."),
            new ElementRule(WAVELENGTH,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                    "The wavelength between successive population crashes."),
            new ElementRule(OFFSET,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                    "The proportion of the last growth phase that is not achieved at the final sample time."),
            XMLUnits.SYNTAX_RULES[0]
    };
}
