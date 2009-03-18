package dr.evomodel.beagle.parsers;

import dr.xml.*;
import dr.evolution.alignment.PatternList;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.beagle.treelikelihood.BeagleTreeLikelihood;
import dr.evomodel.beagle.sitemodel.SiteRateModel;
import dr.evomodel.beagle.sitemodel.BranchSiteModel;
import dr.evomodel.beagle.sitemodel.HomogenousBranchSiteModel;
import dr.evomodel.beagle.sitemodel.GammaSiteRateModel;
import dr.inference.model.Likelihood;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public class TreeLikelihoodParser extends AbstractXMLObjectParser {

    public static final String TREE_LIKELIHOOD = "treeLikelihood";
    public static final String USE_AMBIGUITIES = "useAmbiguities";
    public static final String DEVICE_NUMBER = "deviceNumber";

    public String getParserName() {
        return TREE_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean useAmbiguities = xo.getAttribute(USE_AMBIGUITIES, false);
        int deviceNumber = xo.getAttribute(DEVICE_NUMBER,1) - 1;

        PatternList patternList = (PatternList) xo.getChild(PatternList.class);
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        GammaSiteRateModel siteRateModel = (GammaSiteRateModel) xo.getChild(GammaSiteRateModel.class);

        BranchSiteModel branchSiteModel = new HomogenousBranchSiteModel(
                siteRateModel.getSubstitutionModel(),
                siteRateModel.getSubstitutionModel().getFrequencyModel());

        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        return new BeagleTreeLikelihood(
                patternList,
                treeModel,
                branchSiteModel,
                siteRateModel,
                branchRateModel,
                useAmbiguities,
                deviceNumber
        );
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of a patternlist on a tree given the site model.";
    }

    public Class getReturnType() {
        return Likelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newBooleanRule(USE_AMBIGUITIES, true),
            AttributeRule.newIntegerRule(DEVICE_NUMBER,true),
            new ElementRule(PatternList.class),
            new ElementRule(TreeModel.class),
            new ElementRule(GammaSiteRateModel.class),
            new ElementRule(BranchRateModel.class, true)
    };
};
