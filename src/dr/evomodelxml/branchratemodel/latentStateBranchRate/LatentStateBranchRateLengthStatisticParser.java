package dr.evomodelxml.branchratemodel.latentStateBranchRate;


import dr.evomodel.branchratemodel.SericolaLatentStateBranchRateModel;
import dr.evomodel.branchratemodel.latentStateBranchRate.LatentStateBranchRateLengthStatistic;
import dr.evomodel.tree.TreeModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class LatentStateBranchRateLengthStatisticParser extends AbstractXMLObjectParser {

    private static final String STATE="state";
    
    public String getParserName() {
        return LatentStateBranchRateLengthStatistic.LATENT_STATE_BRANCH_RATE_LENGTH_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        SericolaLatentStateBranchRateModel latenBranchRateModel = (SericolaLatentStateBranchRateModel) xo.getChild(SericolaLatentStateBranchRateModel.class);
        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        LatentStateBranchRateLengthStatistic.STATE state = LatentStateBranchRateLengthStatistic.STATE.valueOf(xo.getAttribute(STATE, LatentStateBranchRateLengthStatistic.STATE.REPLICATING.name()).toUpperCase());
        return new LatentStateBranchRateLengthStatistic(latenBranchRateModel, tree, state);
    }
  public String getParserDescription() {
        return "This element provides a statistic for the length of a tree spent in a replicating or latent state.";
    }

    public Class getReturnType()  {
        return LatentStateBranchRateLengthStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(SericolaLatentStateBranchRateModel.class, "A branch rate model to provide the rates for the non-latent state"),
            new ElementRule(TreeModel.class, "The tree on which this will operate"),
    };

}
