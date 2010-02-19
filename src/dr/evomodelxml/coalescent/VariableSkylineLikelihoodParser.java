package dr.evomodelxml.coalescent;

import dr.evomodel.coalescent.VariableSkylineLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;


/**
 */
public class VariableSkylineLikelihoodParser extends AbstractXMLObjectParser {

    public static final String SKYLINE_LIKELIHOOD = "ovariableSkyLineLikelihood";
    public static final String POPULATION_SIZES = "populationSizes";
    public static final String INDICATOR_PARAMETER = "indicators";
    public static final String LOG_SPACE = "logUnits";

    public static final String TYPE = "type";
    public static final String STEPWISE = "stepwise";
    public static final String LINEAR = "linear";
    public static final String EXPONENTIAL = "exponential";

    public String getParserName() {
        return SKYLINE_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(POPULATION_SIZES);
        Parameter param = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(INDICATOR_PARAMETER);
        Parameter param2 = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(CoalescentLikelihoodParser.POPULATION_TREE);
        TreeModel treeModel = (TreeModel) cxo.getChild(TreeModel.class);

        VariableSkylineLikelihood.Type type = VariableSkylineLikelihood.Type.STEPWISE;
        /* if (xo.hasAttribute(LINEAR) && !xo.getBooleanAttribute(LINEAR))
        {
            type = VariableSkylineLikelihood.Type.STEPWISE;
        }*/

        if (xo.hasAttribute(TYPE)) {
            final String s = xo.getStringAttribute(TYPE);
            if (s.equalsIgnoreCase(STEPWISE)) {
                type = VariableSkylineLikelihood.Type.STEPWISE;
            } else if (s.equalsIgnoreCase(LINEAR)) {
                type = VariableSkylineLikelihood.Type.LINEAR;
            } else if (s.equalsIgnoreCase(EXPONENTIAL)) {
                type = VariableSkylineLikelihood.Type.EXPONENTIAL;
            } else {
                throw new XMLParseException("Unknown Bayesian Skyline type: " + s);
            }
        }

        boolean logSpace = xo.getBooleanAttribute(LOG_SPACE);

        Logger.getLogger("dr.evomodel").info("Variable skyline plot: " + type.toString() + " control points");

        return new VariableSkylineLikelihood(treeModel, param, param2, type, logSpace);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of the tree given the population size vector.";
    }

    public Class getReturnType() {
        return VariableSkylineLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(TYPE, true),
            new ElementRule(POPULATION_SIZES, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(INDICATOR_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(CoalescentLikelihoodParser.POPULATION_TREE, new XMLSyntaxRule[]{
                    new ElementRule(TreeModel.class)
            }),
            AttributeRule.newBooleanRule(LOG_SPACE)
    };

}
