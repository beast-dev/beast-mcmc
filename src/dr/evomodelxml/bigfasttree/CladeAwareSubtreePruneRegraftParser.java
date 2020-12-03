package dr.evomodelxml.bigfasttree;

import dr.evomodel.bigfasttree.CladeAwareSubtreePruneRegraft;
import dr.evomodel.bigfasttree.CladeNodeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;


public class CladeAwareSubtreePruneRegraftParser extends AbstractXMLObjectParser {
    public static final String CLADE_AWARE_SUBTREE_PRUNE_REGRAFT = "cladeAwareSubtreePruneRegraft";
    @Override
    public Object parseXMLObject(XMLObject  xo) throws XMLParseException {
        CladeNodeModel cladeModel =  (CladeNodeModel) xo.getChild(CladeNodeModel.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        return new CladeAwareSubtreePruneRegraft(cladeModel,weight);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "An operator that prunes a regrafts within a clade.";
    }

    @Override
    public Class getReturnType() {
        return CladeAwareSubtreePruneRegraft.class;
    }

    @Override
    public String getParserName() {
        return CLADE_AWARE_SUBTREE_PRUNE_REGRAFT;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(CladeNodeModel.class)
    };
}
