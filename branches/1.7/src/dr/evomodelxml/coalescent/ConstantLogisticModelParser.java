package dr.evomodelxml.coalescent;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.ConstantLogisticModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a ExponentialGrowth.
 */
public class ConstantLogisticModelParser extends AbstractXMLObjectParser {

    public static final String CONSTANT_LOGISTIC_MODEL = "constantLogistic";
    private static final String POPULATION_SIZE = "populationSize";
    private static final String ANCESTRAL_POPULATION_SIZE = "ancestralPopulationSize";

    private static final String GROWTH_RATE = "growthRate";
    private static final String SHAPE = "shape";
    private static final String ALPHA = "alpha";

    public String getParserName() {
        return CONSTANT_LOGISTIC_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(POPULATION_SIZE);
        Parameter N0Param = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(ANCESTRAL_POPULATION_SIZE);
        Parameter N1Param = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(GROWTH_RATE);
        Parameter rParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(SHAPE);
        Parameter cParam = (Parameter) cxo.getChild(Parameter.class);

        double alpha = xo.getDoubleAttribute(ALPHA);

        return new ConstantLogisticModel(N0Param, N1Param, rParam, cParam, alpha, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A demographic model of constant population size followed by logistic growth.";
    }

    public Class getReturnType() {
        return ConstantLogisticModel.class;
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
            new ElementRule(GROWTH_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(SHAPE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            AttributeRule.newDoubleRule(ALPHA)
    };

}
