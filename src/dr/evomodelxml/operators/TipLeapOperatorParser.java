package dr.evomodelxml.operators;

import dr.evolution.util.Taxa;
import dr.evomodel.operators.TipLeapOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class TipLeapOperatorParser extends AbstractXMLObjectParser {

    public static final String TIP_LEAP = "tipLeap";

    public static final String SIZE = "size";
    public static final String TARGET_ACCEPTANCE = "targetAcceptance";
    public static final String TAXA = "taxa";

    public String getParserName() {
        return TIP_LEAP;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        AdaptationMode mode = AdaptationMode.parseMode(xo);

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        Taxa taxa = (Taxa) xo.getChild(Taxa.class);

        // size attribute is mandatory
        final double size = xo.getAttribute(SIZE, Double.NaN);
        final double targetAcceptance = xo.getAttribute(TARGET_ACCEPTANCE, 0.234);

        if (size <= 0.0) {
            throw new XMLParseException("The TipLeap size attribute must be positive and non-zero.");
        }

        if (targetAcceptance <= 0.0 || targetAcceptance >= 1.0) {
            throw new XMLParseException("Target acceptance probability has to lie in (0, 1)");
        }
        TipLeapOperator operator = new TipLeapOperator(treeModel, taxa, weight, size, targetAcceptance, mode);

        return operator;
    }

    public String getParserDescription() {
        return "An operator that moves a tip a certain patristic distance.";
    }

    public Class getReturnType() {
        return TipLeapOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newDoubleRule(SIZE, false),
            AttributeRule.newDoubleRule(TARGET_ACCEPTANCE, true),
            AttributeRule.newBooleanRule(AdaptableMCMCOperator.AUTO_OPTIMIZE, true),
            new ElementRule(TreeModel.class),
            new ElementRule(Taxa.class)
    };

}
