package dr.evomodelxml.operators;

import dr.evolution.alignment.PatternList;
import dr.evomodel.operators.ImportanceNarrowExchange;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class ImportanceNarrowExchangeParser extends AbstractXMLObjectParser {

    public static final String INS = "ImportanceNarrowExchange";
    public static final String EPSILON = "epsilon";

    public String getParserName() {
        return INS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        final PatternList patterns = (PatternList) xo.getChild(PatternList.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final double epsilon = xo.getAttribute(EPSILON, 0.1);

        try {
            return new ImportanceNarrowExchange(treeModel, patterns, epsilon, weight);
        } catch( Exception e ) {
            throw new XMLParseException(e.getMessage());
        }
    }

    // ************************************************************************
    // AbstractXMLObjectParser
    // implementation
    // ************************************************************************

    public String getParserDescription() {
        return "This element represents a swap operator. "
                + "This operator swaps a random subtree with its uncle.";
    }

    public Class getReturnType() {
        return ImportanceNarrowExchange.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newDoubleRule(EPSILON, true),
            new ElementRule(PatternList.class),
            new ElementRule(TreeModel.class)
    };
}
