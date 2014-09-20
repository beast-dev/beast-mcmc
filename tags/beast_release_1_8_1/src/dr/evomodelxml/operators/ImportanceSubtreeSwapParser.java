package dr.evomodelxml.operators;

import dr.evomodel.operators.ImportanceSubtreeSwap;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class ImportanceSubtreeSwapParser extends AbstractXMLObjectParser {

    public static final String IMPORTANCE_SUBTREE_SWAP = "ImportanceSubtreeSwap";

    public String getParserName() {
        return IMPORTANCE_SUBTREE_SWAP;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final int samples = xo.getIntegerAttribute("samples");

        return new ImportanceSubtreeSwap(treeModel, weight, samples);
    }

    // ************************************************************************
    // AbstractXMLObjectParser implementation
    // ************************************************************************

    public String getParserDescription() {
        return "This element represents a importance guided subtree swap operator. "
                + "This operator swaps a random subtree with a second subtree guided by an importance distribution.";
    }

    public Class getReturnType() {
        return ImportanceSubtreeSwap.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules;{
        rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newIntegerRule("samples"),
                new ElementRule(TreeModel.class)
        };
    }

}
