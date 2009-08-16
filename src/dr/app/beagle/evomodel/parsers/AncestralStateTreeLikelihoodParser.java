package dr.app.beagle.evomodel.parsers;

import dr.app.beagle.evomodel.sitemodel.BranchSiteModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.sitemodel.HomogenousBranchSiteModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.xml.*;

/**
 * @author Marc Suchard
 */

public class AncestralStateTreeLikelihoodParser extends AbstractXMLObjectParser {

    public static final String RECONSTRUCTING_TREE_LIKELIHOOD = "ancestralTreeLikelihood";
    public static final String RECONSTRUCTION_TAG = "state";
    public static final String TAG_NAME = "tagName";

    public String getParserName() {
        return RECONSTRUCTING_TREE_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean useAmbiguities = xo.getAttribute(TreeLikelihoodParser.USE_AMBIGUITIES, false);

        PatternList patternList = (PatternList) xo.getChild(PatternList.class);
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        GammaSiteRateModel siteRateModel = (GammaSiteRateModel) xo.getChild(GammaSiteRateModel.class);

        BranchSiteModel branchSiteModel = new HomogenousBranchSiteModel(
                siteRateModel.getSubstitutionModel(),
                siteRateModel.getSubstitutionModel().getFrequencyModel());

        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        DataType dataType = ((SubstitutionModel) xo.getChild(SubstitutionModel.class)).getDataType();

        // default tag is RECONSTRUCTION_TAG
        String tag = xo.getAttribute(TAG_NAME, RECONSTRUCTION_TAG);

        return new AncestralStateBeagleTreeLikelihood(
                patternList,
                treeModel,
                branchSiteModel,
                siteRateModel,
                branchRateModel,
                useAmbiguities,
                dataType,
                tag
        );

    }

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
            AttributeRule.newBooleanRule(TreeLikelihoodParser.USE_AMBIGUITIES, true),
            AttributeRule.newStringRule(TAG_NAME, true),
            new ElementRule(PatternList.class),
            new ElementRule(TreeModel.class),
            new ElementRule(GammaSiteRateModel.class),
            new ElementRule(BranchRateModel.class, true)
    };
}
