package dr.evomodelxml.operators;

import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class WideExchangeOperatorParser extends AbstractXMLObjectParser {
    public static final String WIDE_EXCHANGE = "wideExchange";

    public String getParserName() {
        return WIDE_EXCHANGE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        return new ExchangeOperator(ExchangeOperator.WIDE, treeModel, weight);
    }

    // ************************************************************************
    // AbstractXMLObjectParser implementation
    // ************************************************************************

    public String getParserDescription() {
        return "This element represents a wide exchange operator. "
                + "This operator swaps two random subtrees.";
    }

    public Class getReturnType() {
        return ExchangeOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules;{
        rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                new ElementRule(TreeModel.class)};
    }
}
