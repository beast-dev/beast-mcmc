package dr.evomodelxml.operators;

import dr.evomodel.operators.TreeUniform;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class TreeUniformParser extends AbstractXMLObjectParser {

    public static final String TREE_UNIFORM = "treeUniform";
    public static final String COUNT = "count";

    public String getParserName() {
        return TREE_UNIFORM;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final int n = xo.getAttribute(COUNT, 2);
        if( ! ( n == 2 || n == 3) ) {
          throw new XMLParseException("Sorry, only moves of 2 or 3 nodes implemented.");
        }
        return new TreeUniform(n, treeModel, weight);
    }

    // ************************************************************************
    // AbstractXMLObjectParser implementation
    // ************************************************************************

    public String getParserDescription() {
        return "Simultanouesly change height of two nodes.";
    }

    public Class getReturnType() {
        return TreeUniform.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newIntegerRule(COUNT, true),
            new ElementRule(TreeModel.class)
    };

}
