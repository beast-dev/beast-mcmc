package dr.evomodelxml.branchratemodel;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.NonParametricBranchRateModel2;
import dr.evomodel.tree.TreeModel;
import dr.math.IntegratedTransformedSplines;
import dr.inference.model.Parameter;
import dr.xml.*;



public class NonParametricBranchRateModelParser2 extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "nonparametricRates2";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        IntegratedTransformedSplines splines =
                (IntegratedTransformedSplines) xo.getChild(IntegratedTransformedSplines.class);
        String id = xo.hasId() ? xo.getId() : PARSER_NAME;

        return new NonParametricBranchRateModel2(id, tree, splines);
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
            new ElementRule(IntegratedTransformedSplines.class)
    };
}