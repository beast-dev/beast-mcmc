package dr.evomodelxml.operators;

import dr.evomodel.operators.ImportancePruneAndRegraft;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class ImportancePruneAndRegraftParser extends AbstractXMLObjectParser {

    public static final String IMPORTANCE_PRUNE_AND_REGRAFT = "ImportancePruneAndRegraft";

    public String getParserName() {
        return IMPORTANCE_PRUNE_AND_REGRAFT;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        int samples = xo.getIntegerAttribute("samples");

        return new ImportancePruneAndRegraft(treeModel, weight, samples);
    }

    // ************************************************************************
    // AbstractXMLObjectParser implementation
    // ************************************************************************

    public String getParserDescription() {
        return "This element represents a importance guided prune and regraft operator. "
                + "This operator prunes a random subtree and regrafts it below a node chosen by an importance distribution.";
    }

    public Class getReturnType() {
        return ImportancePruneAndRegraft.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newIntegerRule("samples"),
            new ElementRule(TreeModel.class)
    };

}
