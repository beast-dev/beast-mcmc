package dr.evomodelxml.operators;

import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class NarrowExchangeOperatorParser extends AbstractXMLObjectParser {

    public static final String NARROW_EXCHANGE = "narrowExchange";

    public String getParserName() {
        return NARROW_EXCHANGE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        return new ExchangeOperator(ExchangeOperator.NARROW, treeModel, weight);
    }

    // ************************************************************************
    // AbstractXMLObjectParser implementation
    // ************************************************************************

    public String getParserDescription() {
        return "This element represents a narrow exchange operator. "
                + "This operator swaps a random subtree with its uncle.";
    }

    public Class getReturnType() {
        return ExchangeOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(TreeModel.class)
    };
}
