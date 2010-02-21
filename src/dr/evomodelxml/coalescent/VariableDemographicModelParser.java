package dr.evomodelxml.coalescent;

import dr.evomodel.coalescent.VariableDemographicModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.speciation.SpeciesBindingsParser;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 */
public class VariableDemographicModelParser extends AbstractXMLObjectParser {
    public static final String MODEL_NAME = "variableDemographic";
    public static final String POPULATION_SIZES = "populationSizes";
    public static final String INDICATOR_PARAMETER = "indicators";
    public static final String POPULATION_TREES = "trees";
    private static final String PLOIDY = SpeciesBindingsParser.PLOIDY;
    public static final String POP_TREE = "ptree";

    public static final String LOG_SPACE = "logUnits";
    public static final String USE_MIDPOINTS = "useMidpoints";

    public static final String TYPE = "type";
    //public static final String STEPWISE = "stepwise";
    //public static final String LINEAR = "linear";
    //public static final String EXPONENTIAL = "exponential";

    public static final String demoElementName = "demographic";

    public String getParserName() {
        return MODEL_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(POPULATION_SIZES);
        Parameter popParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(INDICATOR_PARAMETER);
        Parameter indicatorParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(POPULATION_TREES);

        final int nc = cxo.getChildCount();
        TreeModel[] treeModels = new TreeModel[nc];
        double[] populationFactor = new double[nc];

        for (int k = 0; k < treeModels.length; ++k) {
            final XMLObject child = (XMLObject) cxo.getChild(k);
            populationFactor[k] = child.hasAttribute(PLOIDY) ? child.getDoubleAttribute(PLOIDY) : 1.0;

            treeModels[k] = (TreeModel) child.getChild(TreeModel.class);
        }

        VariableDemographicModel.Type type = VariableDemographicModel.Type.STEPWISE;

        if (xo.hasAttribute(TYPE)) {
            final String s = xo.getStringAttribute(TYPE);
            if (s.equalsIgnoreCase(VariableDemographicModel.Type.STEPWISE.toString())) {
                type = VariableDemographicModel.Type.STEPWISE;
            } else if (s.equalsIgnoreCase(VariableDemographicModel.Type.LINEAR.toString())) {
                type = VariableDemographicModel.Type.LINEAR;
            } else if (s.equalsIgnoreCase(VariableDemographicModel.Type.EXPONENTIAL.toString())) {
                type = VariableDemographicModel.Type.EXPONENTIAL;
            } else {
                throw new XMLParseException("Unknown Bayesian Skyline type: " + s);
            }
        }

        final boolean logSpace = xo.getAttribute(LOG_SPACE, false) || type == VariableDemographicModel.Type.EXPONENTIAL;
        final boolean useMid = xo.getAttribute(USE_MIDPOINTS, false);

        Logger.getLogger("dr.evomodel").info("Variable demographic: " + type.toString() + " control points");

        return new VariableDemographicModel(treeModels, populationFactor, popParam, indicatorParam, type,
                logSpace, useMid);

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of the tree given the population size vector.";
    }

    public Class getReturnType() {
        return VariableDemographicModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TYPE, true),
            AttributeRule.newBooleanRule(LOG_SPACE, true),
            AttributeRule.newBooleanRule(USE_MIDPOINTS, true),

            new ElementRule(POPULATION_SIZES, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(INDICATOR_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(POPULATION_TREES, new XMLSyntaxRule[]{
                    new ElementRule(POP_TREE, new XMLSyntaxRule[]{
                            AttributeRule.newDoubleRule(PLOIDY, true),
                            new ElementRule(TreeModel.class),
                    }, 1, Integer.MAX_VALUE)
            })
    };
}
