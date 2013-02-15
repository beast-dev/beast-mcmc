package dr.evomodelxml.coalescent;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.ConstantLogisticModel;
import dr.evomodel.coalescent.ExponentialLogisticModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a ExponentialGrowth.
 */
public class ExponentialLogisticModelParser extends AbstractXMLObjectParser {

    public static final String EXPONENTIAL_LOGISTIC_MODEL = "exponentialLogistic";
    public static final String POPULATION_SIZE = "populationSize";
    public static final String TRANSITION_TIME = "transitionTime";

    public static final String LOGISTIC_GROWTH_RATE = "logisticGrowthRate";
    public static final String LOGISTIC_SHAPE = "logisticShape";
    public static final String EXPONENTIAL_GROWTH_RATE = "exponentialGrowthRate";

    public static final String ALPHA = "alpha";

    public String getParserName() {
        return EXPONENTIAL_LOGISTIC_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(POPULATION_SIZE);
        Parameter N0Param = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(LOGISTIC_GROWTH_RATE);
        Parameter logGrowthParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(LOGISTIC_SHAPE);
        Parameter shapeParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(EXPONENTIAL_GROWTH_RATE);
        Parameter expGrowthParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(TRANSITION_TIME);
        Parameter timeParam = (Parameter) cxo.getChild(Parameter.class);

        double alpha = xo.getDoubleAttribute(ALPHA);

        return new ExponentialLogisticModel(N0Param, logGrowthParam, shapeParam, expGrowthParam, timeParam, alpha, units);
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

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            XMLUnits.SYNTAX_RULES[0],
            new ElementRule(POPULATION_SIZE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(LOGISTIC_GROWTH_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(LOGISTIC_SHAPE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(EXPONENTIAL_GROWTH_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(TRANSITION_TIME,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            AttributeRule.newDoubleRule(ALPHA)
    };
}
