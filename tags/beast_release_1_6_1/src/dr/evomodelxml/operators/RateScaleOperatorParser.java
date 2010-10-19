package dr.evomodelxml.operators;

import dr.evomodel.operators.RateScaleOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class RateScaleOperatorParser extends AbstractXMLObjectParser {

    public static final String SCALE_OPERATOR = "rateScaleOperator";
    public static final String SCALE_FACTOR = "scaleFactor";
    public static final String NO_ROOT = "noRoot";

    public String getParserName() {
        return SCALE_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        CoercionMode mode = CoercionMode.parseMode(xo);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);

        final boolean noRoot = xo.getBooleanAttribute(NO_ROOT);

        if (scaleFactor <= 0.0 || scaleFactor >= 1.0) {
            throw new XMLParseException("scaleFactor must be between 0.0 and 1.0");
        }

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

        RateScaleOperator operator = new RateScaleOperator(treeModel, scaleFactor, noRoot, mode);
        operator.setWeight(weight);
        return operator;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a rateScale operator on a given parameter.";
    }

    public Class getReturnType() {
        return RateScaleOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(SCALE_FACTOR),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(CoercableMCMCOperator.AUTO_OPTIMIZE, true),
            AttributeRule.newBooleanRule(NO_ROOT, true),
            new ElementRule(TreeModel.class),
    };

}
