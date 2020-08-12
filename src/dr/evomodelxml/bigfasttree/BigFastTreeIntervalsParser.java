package dr.evomodelxml.bigfasttree;

import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;


public class BigFastTreeIntervalsParser extends AbstractXMLObjectParser {

    public static final String TREE_INTERVALS = "bigFastTreeIntervals";
    public static final String TREE = "tree";

    public String getParserName() {
        return TREE_INTERVALS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        return new BigFastTreeIntervals(tree);

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Parser for bigFastTreeIntervals.";
    }

    public Class getReturnType() {
        return BigFastTreeIntervals.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class)
    };

}