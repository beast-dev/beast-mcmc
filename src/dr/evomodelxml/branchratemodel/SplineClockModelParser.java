package dr.evomodelxml.branchratemodel;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.SplineClockModel;
import dr.evomodel.tree.TreeModel;
import dr.math.IntegratedTransformedSplines;
import dr.inference.model.Parameter;
import dr.xml.*;



public class SplineClockModelParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "splineClockModel";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        IntegratedTransformedSplines splines =
                (IntegratedTransformedSplines) xo.getChild(IntegratedTransformedSplines.class);
        String id = xo.hasId() ? xo.getId() : PARSER_NAME;

        return new SplineClockModel(id, tree, splines);
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
        return SplineClockModelParser.class;
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