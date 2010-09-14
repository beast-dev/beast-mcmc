package dr.evomodelxml.operators;

import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class ExchangeOperatorParser {

    public static final String NARROW_EXCHANGE = "narrowExchange";
    public static final String WIDE_EXCHANGE = "wideExchange";

    public static XMLObjectParser NARROW_EXCHANGE_OPERATOR_PARSER = new AbstractXMLObjectParser() {
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
    };

    public static XMLObjectParser WIDE_EXCHANGE_OPERATOR_PARSER = new AbstractXMLObjectParser() {
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
                    new ElementRule(TreeModel.class)
            };
        }
    };
}
