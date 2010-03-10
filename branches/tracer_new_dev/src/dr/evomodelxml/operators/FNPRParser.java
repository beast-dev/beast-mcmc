package dr.evomodelxml.operators;

import dr.evomodel.operators.FNPR;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class FNPRParser extends AbstractXMLObjectParser {

    public static final String FNPR = "FixedNodeheightSubtreePruneRegraft";

    public String getParserName() {
        return FNPR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        return new FNPR(treeModel, weight);
    }

    // ************************************************************************
    // AbstractXMLObjectParser implementation
    // ************************************************************************

    public String getParserDescription() {
        return "This element represents a FNPR operator. "
                + "This operator swaps a random subtree with its uncle.";
    }

    public Class getReturnType() {
        return FNPR.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(TreeModel.class)
    };
}
