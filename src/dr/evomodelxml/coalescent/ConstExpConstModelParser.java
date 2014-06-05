package dr.evomodelxml.coalescent;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.ConstExpConstModel;
import dr.evomodel.coalescent.ConstantExponentialModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a ConstantExponentialModel.
 */
public class ConstExpConstModelParser extends AbstractXMLObjectParser {

    public static String CONST_EXP_CONST_MODEL = "constExpConst";
    public static String POPULATION_SIZE = "populationSize";
    public static String ANCESTRAL_POPULATION_SIZE = "ancestralPopulationSize";
    public static String FINAL_PHASE_START_TIME = "finalPhaseStartTime";
    public static String GROWTH_PHASE_TIME = "growthPhaseTime";

    public String getParserName() {
        return CONST_EXP_CONST_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(POPULATION_SIZE);
        Parameter N0Param = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(ANCESTRAL_POPULATION_SIZE);
        Parameter N1Param = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(FINAL_PHASE_START_TIME);
        Parameter timeParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(GROWTH_PHASE_TIME);
        Parameter epochParam = (Parameter) cxo.getChild(Parameter.class);

        return new ConstExpConstModel(N0Param, N1Param, timeParam, epochParam, units);
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
            new ElementRule(ANCESTRAL_POPULATION_SIZE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(FINAL_PHASE_START_TIME,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(GROWTH_PHASE_TIME,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
    };


}
