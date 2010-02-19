package dr.evomodelxml.coalescent;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.ConstantExponentialModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a ConstantExponentialModel.
 */
public class ConstantExponentialModelParser extends AbstractXMLObjectParser {

    public static final String CONSTANT_EXPONENTIAL_MODEL = "constantExponential";
    public static final String POPULATION_SIZE = "populationSize";
    public static final String GROWTH_PHASE_START_TIME = "growthPhaseStartTime";

    public static final String GROWTH_RATE = "growthRate";
    public static final String DOUBLING_TIME = "doublingTime";

    public String getParserName() {
        return CONSTANT_EXPONENTIAL_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(POPULATION_SIZE);
        Parameter N0Param = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(GROWTH_PHASE_START_TIME);
        Parameter timeParam = (Parameter) cxo.getChild(Parameter.class);

        Parameter rParam;
        boolean usingGrowthRate = true;

        if (xo.getChild(GROWTH_RATE) != null) {
            cxo = xo.getChild(GROWTH_RATE);
            rParam = (Parameter) cxo.getChild(Parameter.class);
        } else {
            cxo = xo.getChild(DOUBLING_TIME);
            rParam = (Parameter) cxo.getChild(Parameter.class);
            usingGrowthRate = false;
        }

        return new ConstantExponentialModel(N0Param, timeParam, rParam, units, usingGrowthRate);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A demographic model of constant population size followed by exponential growth.";
    }

    public Class getReturnType() {
        return ConstantExponentialModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            XMLUnits.SYNTAX_RULES[0],
            new ElementRule(POPULATION_SIZE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(GROWTH_PHASE_START_TIME,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new XORRule(

                    new ElementRule(GROWTH_RATE,
                            new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                            "A value of zero represents a constant population size, negative values represent decline towards the present, " +
                                    "positive numbers represents exponential growth towards the present. " +
                                    "A random walk operator is recommended for this parameter with a starting value of 0.0 and no upper or lower limits."),
                    new ElementRule(DOUBLING_TIME,
                            new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                            "This parameter determines the doubling time.")
            )
    };

}
