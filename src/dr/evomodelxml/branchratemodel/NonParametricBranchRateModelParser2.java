package dr.evomodelxml.branchratemodel;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.NonParametricBranchRateModel2;
import dr.evomodel.tree.TreeModel;
import dr.inference.distribution.RandomField;
import dr.inference.model.Parameter;
import dr.xml.*;



public class NonParametricBranchRateModelParser2 extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "nonparametricRates2";
    private static final String COEFFICIENTS = "coefficients";
    private static final String KNOTS = "knots";
    private static final String DEGREE = "degree";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        Parameter coefficients = (Parameter) xo.getElementFirstChild(COEFFICIENTS);
        double[] knots = xo.getDoubleArrayAttribute(KNOTS);
        int degree = xo.getIntegerAttribute(DEGREE);
        String id = xo.hasId() ? xo.getId() : PARSER_NAME;

        return new NonParametricBranchRateModel2(id, tree, coefficients, knots, degree);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return NonParametricBranchRateModelParser2.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(COEFFICIENTS,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            AttributeRule.newDoubleArrayRule(KNOTS),
            AttributeRule.newIntegerRule(DEGREE),
    };
}