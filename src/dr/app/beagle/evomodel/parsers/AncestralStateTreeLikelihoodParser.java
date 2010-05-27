package dr.app.beagle.evomodel.parsers;

import dr.app.beagle.evomodel.sitemodel.BranchSiteModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.app.beagle.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.app.beagle.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

/**
 * @author Marc Suchard
 * @author Andrew Rambaut
 */

public class AncestralStateTreeLikelihoodParser extends TreeLikelihoodParser {

    public static final String RECONSTRUCTING_TREE_LIKELIHOOD = "ancestralTreeLikelihood";
    public static final String RECONSTRUCTION_TAG = "state";
    public static final String RECONSTRUCTION_TAG_NAME = "stateTagName";
    public static final String MAP_RECONSTRUCTION = "useMAP";
    public static final String MARGINAL_LIKELIHOOD = "useMarginalLikelihood";

    public String getParserName() {
        return RECONSTRUCTING_TREE_LIKELIHOOD;
    }

    protected BeagleTreeLikelihood createTreeLikelihood(PatternList patternList, TreeModel treeModel,
                                                        BranchSiteModel branchSiteModel, GammaSiteRateModel siteRateModel,
                                                        BranchRateModel branchRateModel,
                                                        boolean useAmbiguities, PartialsRescalingScheme scalingScheme,
                                                        XMLObject xo) throws XMLParseException {

        SubstitutionModel substModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);

        DataType dataType = substModel.getDataType();

        // default tag is RECONSTRUCTION_TAG
        String tag = xo.getAttribute(RECONSTRUCTION_TAG_NAME, RECONSTRUCTION_TAG);

        boolean useMAP = xo.getAttribute(MAP_RECONSTRUCTION, false);
        boolean useMarginalLogLikelihood = xo.getAttribute(MARGINAL_LIKELIHOOD, true);

        return new AncestralStateBeagleTreeLikelihood(  // Current just returns a BeagleTreeLikelihood
                patternList,
                treeModel,
                branchSiteModel,
                siteRateModel,
                branchRateModel,
                useAmbiguities,
                scalingScheme,
                dataType,
                tag,
                substModel,
                useMAP,
                useMarginalLogLikelihood
        );
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
            AttributeRule.newBooleanRule(TreeLikelihoodParser.USE_AMBIGUITIES, true),
            AttributeRule.newStringRule(RECONSTRUCTION_TAG_NAME, true),
            new ElementRule(PatternList.class),
            new ElementRule(TreeModel.class),
            new ElementRule(GammaSiteRateModel.class),
            new ElementRule(BranchRateModel.class, true),
            new ElementRule(SubstitutionModel.class),
            AttributeRule.newStringRule(TreeLikelihoodParser.SCALING_SCHEME,true),
        };
    }
}
