package dr.evomodelxml.bigfasttree;

import dr.evomodel.bigfasttree.CladeAwareSubtreePruneRegraft;
import dr.evomodel.bigfasttree.CladeNodeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;


public class CladeAwareSubtreePruneRegraftParser extends AbstractXMLObjectParser {
    public static final String CLADE_AWARE_SUBTREE_PRUNE_REGRAFT = "cladeAwareSubtreePruneRegraft";
    public static final String SPR_PER_CALL = "rearrangementsPerCall";
    @Override
    public Object parseXMLObject(XMLObject  xo) throws XMLParseException {
        CladeNodeModel cladeModel =  (CladeNodeModel) xo.getChild(CladeNodeModel.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        int SPRperCall;
        if (xo.hasAttribute(SPR_PER_CALL)) {
            SPRperCall = xo.getIntegerAttribute(SPR_PER_CALL);
        }else{
            SPRperCall = 1;
        }
        return new CladeAwareSubtreePruneRegraft(cladeModel,weight,SPRperCall);
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
            AttributeRule.newDoubleRule(SPR_PER_CALL, true,"The number of moves per operator call. Each node will be sampled from the first clade without replacement"),
            new ElementRule(CladeNodeModel.class)
    };
}
