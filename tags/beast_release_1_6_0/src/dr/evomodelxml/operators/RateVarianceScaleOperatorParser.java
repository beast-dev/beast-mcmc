package dr.evomodelxml.operators;

import dr.evomodel.operators.RateVarianceScaleOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class RateVarianceScaleOperatorParser extends AbstractXMLObjectParser {

    public static final String SCALE_OPERATOR = "rateVarianceScaleOperator";
    public static final String SCALE_FACTOR = "scaleFactor";

    public String getParserName() {
        return SCALE_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        CoercionMode mode = CoercionMode.parseMode(xo);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);

        if (scaleFactor <= 0.0 || scaleFactor >= 1.0) {
            throw new XMLParseException("scaleFactor must be between 0.0 and 1.0");
        }

        final TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

        final Parameter variance = (Parameter) xo.getChild(Parameter.class);
        if (variance.getDimension() != 1) {
            throw new XMLParseException("dimension of the variance parameter should be 1");
        }

        RateVarianceScaleOperator operator = new RateVarianceScaleOperator(treeModel, variance, scaleFactor, mode);
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
        return RateVarianceScaleOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(SCALE_FACTOR),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(CoercableMCMCOperator.AUTO_OPTIMIZE, true),
            new ElementRule(TreeModel.class),
            new ElementRule(Parameter.class),
    };
}
