package dr.evomodelxml.treelikelihood;

import dr.evolution.alignment.PatternList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TipPartialsModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.xml.*;

/**
 */
public class TreeLikelihoodParser extends AbstractXMLObjectParser {

    public static final String TREE_LIKELIHOOD = "treeLikelihood";
    public static final String ANCESTRAL_TREE_LIKELIHOOD = "ancestralTreeLikelihood";
    public static final String USE_AMBIGUITIES = "useAmbiguities";
    public static final String ALLOW_MISSING_TAXA = "allowMissingTaxa";
    public static final String STORE_PARTIALS = "storePartials";
    public static final String SCALING_FACTOR = "scalingFactor";
    public static final String SCALING_THRESHOLD = "scalingThreshold";
    public static final String FORCE_JAVA_CORE = "forceJavaCore";
    public static final String FORCE_RESCALING = "forceRescaling";


    public String getParserName() {
        return TREE_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean useAmbiguities = xo.getAttribute(USE_AMBIGUITIES, false);
        boolean allowMissingTaxa = xo.getAttribute(ALLOW_MISSING_TAXA, false);
        boolean storePartials = xo.getAttribute(STORE_PARTIALS, true);
        boolean forceJavaCore = xo.getAttribute(FORCE_JAVA_CORE, false);

        if (Boolean.valueOf(System.getProperty("java_only"))) {
            forceJavaCore = true;
        }

        PatternList patternList = (PatternList) xo.getChild(PatternList.class);
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        SiteModel siteModel = (SiteModel) xo.getChild(SiteModel.class);

        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        TipPartialsModel tipPartialsModel = (TipPartialsModel) xo.getChild(TipPartialsModel.class);
        if (tipPartialsModel != null && tipPartialsModel.getPatternList() != null) {
            throw new XMLParseException("A sequence error model cannot be used with multiple partitions");
        }

        boolean forceRescaling = xo.getAttribute(FORCE_RESCALING, false);

        return new TreeLikelihood(
                patternList,
                treeModel,
                siteModel,
                branchRateModel,
                tipPartialsModel,
                useAmbiguities, allowMissingTaxa, storePartials, forceJavaCore, forceRescaling);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of a patternlist on a tree given the site model.";
    }

    public Class getReturnType() {
        return TreeLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(USE_AMBIGUITIES, true),
            AttributeRule.newBooleanRule(ALLOW_MISSING_TAXA, true),
            AttributeRule.newBooleanRule(STORE_PARTIALS, true),
            AttributeRule.newBooleanRule(FORCE_JAVA_CORE, true),
            AttributeRule.newBooleanRule(FORCE_RESCALING, true),
            new ElementRule(PatternList.class),
            new ElementRule(TreeModel.class),
            new ElementRule(SiteModel.class),
            new ElementRule(BranchRateModel.class, true),
            new ElementRule(TipPartialsModel.class, true)
    };
}
