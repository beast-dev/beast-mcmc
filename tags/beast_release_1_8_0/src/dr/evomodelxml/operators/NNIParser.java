package dr.evomodelxml.operators;

import dr.evomodel.operators.NNI;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class NNIParser extends AbstractXMLObjectParser {

    public static final String NNI = "NearestNeighborInterchange";

    public String getParserName() {
        return NNI;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        return new NNI(treeModel, weight);
    }

    // ************************************************************************
    // AbstractXMLObjectParser
    // implementation
    // ************************************************************************

    public String getParserDescription() {
        return "This element represents a NNI operator. "
                + "This operator swaps a random subtree with its uncle.";
    }

    public Class getReturnType() {
        return NNI.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(TreeModel.class)
    };
}
