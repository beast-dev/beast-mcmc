package dr.evomodel.treelikelihood;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.SubstitutionEpochModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.inference.model.Likelihood;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Marc A. Suchard
 */
public class EpochTreeLikelihood extends TreeLikelihood {
    public static final String LIKE_NAME = "epochTreeLikelihood";

	    public EpochTreeLikelihood(PatternList patternList,
                          TreeModel treeModel,
                          SiteModel siteModel,
                          BranchRateModel branchRateModel,
                          TipStatesModel tipStatesModel,
                          boolean useAmbiguities,
                          boolean allowMissingTaxa,
                          boolean storePartials,
                          boolean forceJavaCore) {
		    super(patternList,  treeModel, siteModel, branchRateModel, tipStatesModel,useAmbiguities,allowMissingTaxa,storePartials,forceJavaCore, false);

	    }


	/**
	  * Traverse the tree calculating partial likelihoods.
	  *
	  * @return whether the partials for this node were recalculated.
	  */
	 protected boolean traverse(Tree tree, NodeRef node) {

	     boolean update = false;

	     int nodeNum = node.getNumber();

	     NodeRef parent = tree.getParent(node);

	     // First update the transition probability matrix(ices) for this branch
	     if (parent != null && updateNode[nodeNum]) {

	         final double branchRate = branchRateModel.getBranchRate(tree, node);

	         // Get the operational time of the branch
		     final double parentNodeHeight = tree.getNodeHeight(parent);
		     final double nodeHeight = tree.getNodeHeight(node);
	         final double branchTime = branchRate * (parentNodeHeight - nodeHeight);

	         if (branchTime < 0.0) {
	             throw new RuntimeException("Negative branch length: " + branchTime);
	         }

	         likelihoodCore.setNodeMatrixForUpdate(nodeNum);

	         for (int i = 0; i < categoryCount; i++) {

	             double branchLength = siteModel.getRateForCategory(i) * branchTime;
	             ((SubstitutionEpochModel)siteModel.getSubstitutionModel()).getTransitionProbabilities(nodeHeight, parentNodeHeight,branchLength, probabilities);
	             likelihoodCore.setNodeMatrix(nodeNum, i, probabilities);
	         }

	         update = true;
	     }

	     // If the node is internal, update the partial likelihoods.
	     if (!tree.isExternal(node)) {

	         // Traverse down the two child nodes
	         NodeRef child1 = tree.getChild(node, 0);
	         final boolean update1 = traverse(tree, child1);

	         NodeRef child2 = tree.getChild(node, 1);
	         final boolean update2 = traverse(tree, child2);

	         // If either child node was updated then update this node too
	         if (update1 || update2) {

	             final int childNum1 = child1.getNumber();
	             final int childNum2 = child2.getNumber();

	             likelihoodCore.setNodePartialsForUpdate(nodeNum);

	             if (integrateAcrossCategories) {
	                 likelihoodCore.calculatePartials(childNum1, childNum2, nodeNum);
	             } else {
	                 likelihoodCore.calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
	             }

	             if (parent == null) {
	                 // No parent this is the root of the tree -
	                 // calculate the pattern likelihoods
	                 double[] frequencies = frequencyModel.getFrequencies();

	                 double[] partials = getRootPartials();

	                 likelihoodCore.calculateLogLikelihoods(partials, frequencies, patternLogLikelihoods);
	             }

	             update = true;
	         }
	     }

	     return update;

	 }


    /**
     * The XML parser
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return LIKE_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            boolean useAmbiguities = xo.getAttribute(TreeLikelihoodParser.USE_AMBIGUITIES, false);
            boolean allowMissingTaxa = xo.getAttribute(TreeLikelihoodParser.ALLOW_MISSING_TAXA, false);
            boolean storePartials = xo.getAttribute(TreeLikelihoodParser.STORE_PARTIALS, true);
            boolean forceJavaCore = xo.getAttribute(TreeLikelihoodParser.FORCE_JAVA_CORE, false);

            PatternList patternList = (PatternList) xo.getChild(PatternList.class);
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            SiteModel siteModel = (SiteModel) xo.getChild(SiteModel.class);

	        if (! (siteModel.getSubstitutionModel() instanceof SubstitutionEpochModel)) {
		        throw new XMLParseException("EpochTreeLikelihood only currently works with a SubstitutionEpochModel");
	        }

            BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

            TipStatesModel tipStatesModel = (TipStatesModel) xo.getChild(TipStatesModel.class);

	        Logger.getLogger("dr.evolution").info("\n ---------------------------------\nCreating EpochTreeLikelihood model.");
	          Logger.getLogger("dr.evolution").info("\tIf you publish results using substitution epoch likelihood, please reference" +
			          " Shapiro and Suchard (in preparation).\n---------------------------------\n");


            return new EpochTreeLikelihood(
                    patternList,
                    treeModel,
                    siteModel,
                    branchRateModel,
                    tipStatesModel,
                    useAmbiguities, allowMissingTaxa, storePartials, forceJavaCore);
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
                AttributeRule.newBooleanRule(TreeLikelihoodParser.USE_AMBIGUITIES, true),
                AttributeRule.newBooleanRule(TreeLikelihoodParser.ALLOW_MISSING_TAXA, true),
                AttributeRule.newBooleanRule(TreeLikelihoodParser.STORE_PARTIALS, true),
                AttributeRule.newBooleanRule(TreeLikelihoodParser.FORCE_JAVA_CORE, true),
                new ElementRule(PatternList.class),
                new ElementRule(TreeModel.class),
                new ElementRule(SiteModel.class),
                new ElementRule(BranchRateModel.class, true),
                new ElementRule(TipStatesModel.class, true)
        };
    };

}