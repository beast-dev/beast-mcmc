package dr.evomodelxml.operators;

import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class WilsonBaldingParser extends AbstractXMLObjectParser {

    public static final String WILSON_BALDING = "wilsonBalding";
    public static final String DEMOGRAPHIC_MODEL = "demographicModel";

    public String getParserName() {
        return WILSON_BALDING;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        final TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

        return new WilsonBalding(treeModel, weight);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(TreeModel.class)
    };

    public String getParserDescription() {
        return "An operator which performs the Wilson-Balding move on a tree";
    }

    public Class getReturnType() {
        return WilsonBalding.class;
    }
}
