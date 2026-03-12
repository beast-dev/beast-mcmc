package dr.evomodelxml.branchratemodel;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.LogLinearBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.distribution.RandomField;
import dr.inference.model.Parameter;
import dr.xml.*;



public class LogLinearBranchRateModelParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "logLinearBranchRate";
    private static final String SLOPE = "slope";
    private static final String INTERCEPT = "intercept";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        Parameter slope = (Parameter) xo.getElementFirstChild(SLOPE);
        Parameter intercept = (Parameter) xo.getElementFirstChild(INTERCEPT);
        String id = xo.hasId() ? xo.getId() : PARSER_NAME;

        return new LogLinearBranchRateModel(id, tree, slope, intercept);
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
        return NonParametricBranchRateModelParser.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(SLOPE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(INTERCEPT,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
    };
}