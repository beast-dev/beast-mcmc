package dr.evomodelxml.indel;

import dr.evolution.alignment.Alignment;
import dr.evomodel.indel.TKF91Likelihood;
import dr.evomodel.indel.TKF91Model;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

/**
 *
 */
public class TKF91LikelihoodParser extends AbstractXMLObjectParser {

    public static final String TKF91_LIKELIHOOD = "tkf91Likelihood";
    public static final String TKF91_DEATH = "deathRate";
    //public static final String MU = "mutationRate";

    public String getParserName() {
        return TKF91_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        Alignment alignment = (Alignment) xo.getChild(Alignment.class);
        GammaSiteModel siteModel = (GammaSiteModel) xo.getChild(GammaSiteModel.class);
        TKF91Model tkfModel = (TKF91Model) xo.getChild(TKF91Model.class);
        return new TKF91Likelihood(tree, alignment, siteModel, tkfModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Returns the total likelihood of a single alignment under the TKF91 model, for a given tree. " +
                "In particular all possible ancestral histories of insertions and deletions leading to the " +
                "alignment of sequences at the tips are taken into account.";
    }

    public Class getReturnType() {
        return TKF91Likelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(Alignment.class),
            new ElementRule(GammaSiteModel.class),
            new ElementRule(TKF91Model.class)
    };
}
