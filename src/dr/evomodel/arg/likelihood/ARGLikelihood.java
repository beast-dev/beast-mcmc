/*
 * TreeLikelihood.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */

package dr.evomodel.arg.likelihood;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.arg.ARGModel;
import dr.evomodel.arg.ARGTree;
import dr.evomodel.arg.operators.ARGPartitioningOperator;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.treelikelihood.*;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * ARGLikelihood - implements a Likelihood Function for sequences on an ancestral recombination graph.
 *
 * @author Marc Suchard
 * @version $Id: ARGLikelihood.java,v 1.3 2006/10/23 04:13:41 msuchard Exp $
 */

public class ARGLikelihood extends AbstractARGLikelihood {

	public static final String ARG_LIKELIHOOD = "argTreeLikelihood";
	public static final String USE_AMBIGUITIES = "useAmbiguities";
	public static final String STORE_PARTIALS = "storePartials";
	public static final String USE_SCALING = "useScaling";

	/**
	 * Constructor.
	 */
	public ARGLikelihood(PatternList patternList,
	                     ARGModel treeModel,
	                     SiteModel siteModel,
	                     BranchRateModel branchRateModel,
	                     boolean useAmbiguities,
	                     boolean storePartials,
	                     boolean useScaling) {

		super(ARG_LIKELIHOOD, patternList, treeModel);

		partition = treeModel.addLikelihoodCalculator(this);

		this.storePartials = storePartials;
		this.useAmbiguities = useAmbiguities;

		try {
			this.siteModel = siteModel;
			addModel(siteModel);

			this.frequencyModel = siteModel.getFrequencyModel();
			addModel(frequencyModel);

			integrateAcrossCategories = siteModel.integrateAcrossCategories();

			this.categoryCount = siteModel.getCategoryCount();

			if (integrateAcrossCategories) {
				if (patternList.getDataType() instanceof dr.evolution.datatype.Nucleotides) {

					if (NativeNucleotideLikelihoodCore.isAvailable()) {

						Logger.getLogger("dr.evomodel").info("TreeLikelihood using native nucleotide likelihood core");
						likelihoodCore = new NativeNucleotideLikelihoodCore();
					} else {

						Logger.getLogger("dr.evomodel").info("TreeLikelihood using Java nucleotide likelihood core");
						likelihoodCore = new NucleotideLikelihoodCore();
					}

				} else if (patternList.getDataType() instanceof dr.evolution.datatype.AminoAcids) {

					if (NativeAminoAcidLikelihoodCore.isAvailable()) {
						Logger.getLogger("dr.evomodel").info("TreeLikelihood using native amino acid likelihood core");
						likelihoodCore = new NativeAminoAcidLikelihoodCore();
					} else {
						Logger.getLogger("dr.evomodel").info("TreeLikelihood using java likelihood core");
						likelihoodCore = new AminoAcidLikelihoodCore();
					}
				} else if (patternList.getDataType() instanceof dr.evolution.datatype.Codons) {
					Logger.getLogger("dr.evomodel").info("TreeLikelihood using Java codon likelihood core");
//					likelihoodCore = new CodonLikelihoodCore(patternList.getStateCount());
					useAmbiguities = true;
					throw new RuntimeException("Still need to merge codon likelihood core");
				} else {
					if (patternList.getDataType() instanceof dr.evolution.datatype.HiddenNucleotides &&
							NativeCovarionLikelihoodCore.isAvailable()) {
						Logger.getLogger("dr.evomodel").info("TreeLikelihood using native covarion likelihood core");
						likelihoodCore = new NativeCovarionLikelihoodCore();
					} else {
						Logger.getLogger("dr.evomodel").info("TreeLikelihood using Java general likelihood core");
						likelihoodCore = new GeneralLikelihoodCore(patternList.getStateCount());
					}

				}
			} else {

				Logger.getLogger("dr.evomodel").info("TreeLikelihood using Java general likelihood core");
				likelihoodCore = new GeneralLikelihoodCore(patternList.getStateCount());
			}
			Logger.getLogger("dr.evomodel").info("  " + (useAmbiguities ? "Using" : "Ignoring") + " ambiguities in tree likelihood.");
			Logger.getLogger("dr.evomodel").info("  Partial likelihood scaling " + (useScaling ? "on." : "off."));

			if (branchRateModel != null) {
				this.branchRateModel = branchRateModel;
				Logger.getLogger("dr.evomodel").info("Branch rate model used: " + branchRateModel.getModelName());
			} else {
				this.branchRateModel = new DefaultBranchRateModel();
			}
			addModel(this.branchRateModel);

			probabilities = new double[stateCount * stateCount];

//			likelihoodCore.initialize(nodeCount, patternCount, categoryCount, integrateAcrossCategories, useScaling);
			likelihoodCore.initialize(nodeCount, patternCount, categoryCount, integrateAcrossCategories);

			int extNodeCount = treeModel.getExternalNodeCount();
			int intNodeCount = treeModel.getInternalNodeCount();

			for (int i = 0; i < extNodeCount; i++) {
				// Find the id of tip i in the patternList
				String id = treeModel.getTaxonId(i);
				int index = patternList.getTaxonIndex(id);

//                System.err.println("id = "+id+"  index = "+index);

				if (index == -1) {
					throw new TaxonList.MissingTaxonException("Taxon, " + id + ", in tree, " + treeModel.getId() +
							", is not found in patternList, " + patternList.getId());
				}

				if (useAmbiguities) {
					setPartials(likelihoodCore, patternList, categoryCount, index, i);
				} else {
					setStates(likelihoodCore, patternList, index, i);
				}
			}

//            System.exit(-1);
			for (int i = 0; i < intNodeCount; i++) {
				likelihoodCore.createNodePartials(extNodeCount + i);
			}
		} catch (TaxonList.MissingTaxonException mte) {
			throw new RuntimeException(mte.toString());
		}

	}

	// **************************************************************
	// ModelListener IMPLEMENTATION
	// **************************************************************

	/**
	 * Handles model changed events from the submodels.
	 */
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		// System.err.println("yoyo");
		if (model == treeModel) {
			if (object instanceof ARGModel.TreeChangedEvent) {

				//      if (((ARGModel.TreeChangedEvent)object).isNodeChanged()) {

				//          updateNodeAndChildren(((ARGModel.TreeChangedEvent)object).getNode());

				//      } else {
				updateAllNodes();

				//      }
			}
			if (object instanceof ARGPartitioningOperator.PartitionChangedEvent) {
				ARGPartitioningOperator.PartitionChangedEvent event =
						(ARGPartitioningOperator.PartitionChangedEvent) object;
				Parameter partitioning = event.getParameter();
				// todo update only nodes below node of partitioning;
				updateAllNodes();
			}

		} else if (model == branchRateModel) {
			updateAllNodes();

		} else if (model == frequencyModel) {

			updateAllNodes();

		} else if (model instanceof SiteModel) {

			updateAllNodes();

		} else {

			throw new RuntimeException("Unknown componentChangedEvent");
		}

		super.handleModelChangedEvent(model, object, index);
	}

	// **************************************************************
	// Model IMPLEMENTATION
	// **************************************************************

	/**
	 * Stores the additional state other than model components
	 */
	protected void storeState() {

		if (storePartials) {
			likelihoodCore.storeState();
		}
		super.storeState();

	}

	/**
	 * Restore the additional stored state
	 */
	protected void restoreState() {

//        if (storePartials) {
		//           likelihoodCore.restoreState();
		//       } else {
		updateAllNodes();
		//      }

		super.restoreState();

	}

	// **************************************************************
	// Likelihood IMPLEMENTATION
	// **************************************************************

	/**
	 * Calculate the log likelihood of the current state.
	 *
	 * @return the log likelihood.
	 */
	protected double calculateLogLikelihood() {
//	    sendState(0);
//	    System.exit(0);
//	    System.err.println("here");

//        ARGTree argTree = new ARGTree(treeModel, 0);
//        ARGModel tree = treeModel;

		ARGTree tree = new ARGTree(treeModel, partition);
//        ARGModel argTree = treeModel;

		NodeRef root = tree.getRoot();
//        NodeRef argRoot = argTree.getRoot();

		/*      final int K = argTree.getNodeCount();
				final int J = tree.getNodeCount();
				System.err.printf("%2d %2d",K,J);
				for(int i=0; i<J; i++)
				  System.err.printf(" : %5.3f %5.3f",
						argTree.getNodeHeight(argTree.getNode(i)),
						tree.getNodeHeight(tree.getNode(i))
				 );
				System.err.println();
				System.err.println(argTree.toGraphString());
				System.err.println(tree.toGraphString());*/
//        System.err.println(treeModel.toString())    ;
//        Tree test = new ARGTree(treeModel,partition);

//        System.err.println(Tree.Utils.newick(test));

		final int I = tree.getExternalNodeCount();

		// todo Figure out how to cache mapping between nodes[] and setPartials/setStates internals

		for (int i = 0; i < I; i++) {
			// Find the id of tip i in the patternList
			String id = tree.getTaxonId(i);
			int index = patternList.getTaxonIndex(id);

			if (useAmbiguities) {
				setPartials(likelihoodCore, patternList, categoryCount, index, i);
			} else {
				setStates(likelihoodCore, patternList, index, i);
			}
		}

		/*       System.err.println("ARG:\n"+treeModel.toGraphString());
			  System.err.println("Calc for:"+Tree.Utils.uniqueNewick(tree, root));

			  System.err.println("Starting tree likelihood");*/
		if (rootPartials == null) {
			rootPartials = new double[patternCount * stateCount];
		}

		if (patternLogLikelihoods == null) {
			patternLogLikelihoods = new double[patternCount];
		}

		if (!integrateAcrossCategories) {
			if (siteCategories == null) {
				siteCategories = new int[patternCount];
			}
			for (int i = 0; i < patternCount; i++) {
				siteCategories[i] = siteModel.getCategoryOfSite(i);
			}
		}

		traverse(tree, root);

		//********************************************************************
		// after traverse all nodes and patterns have been updated --
		//so change flags to reflect this.
		for (int i = 0; i < nodeCount; i++) {
			updateNode[i] = false;
		}
		//********************************************************************

		double logL = 0.0;

		for (int i = 0; i < patternCount; i++) {
//        	System.err.printf("Pattern %2d:  %5.4f  %5.4f\n",i,patternLogLikelihoods[i],patternWeights[i]);
			logL += patternLogLikelihoods[i] * patternWeights[i];
		}

		/*      if (Double.isNaN(logL)) {

					System.err.println("ARG:\n"+treeModel.toGraphString());
					System.err.println("ARG Tree:\n"+tree.toGraphString());
					System.err.println("likelihood for partition = "+this.partition);
					throw new RuntimeException("Likelihood NaN");
				}*/

//		likelihoodCore.checkScaling();
//System.err.println(logL);
//        System.err.println("Finished tree likelihood.");
		return logL;
	}

	/**
	 * Traverse the tree calculating partial likelihoods.
	 *
	 * @return whether the partials for this node were recalculated.
	 */
	private final boolean traverse(Tree tree, NodeRef node) {

		boolean update = false;

		int nodeNum = node.getNumber();
//        System.err.println(nodeNum);

		NodeRef parent = tree.getParent(node);

		// First update the transition probability matrix(ices) for this branch
		if (parent != null && updateNode[nodeNum]) {


			double branchRate = branchRateModel.getBranchRate(tree, node);

			// Get the operational time of the branch
			double branchTime = branchRate * (tree.getNodeHeight(parent) - tree.getNodeHeight(node));
			if (branchTime < 0.0) {
				throw new RuntimeException("Negative branch length: " + branchTime);
			}

			for (int i = 0; i < categoryCount; i++) {
                double branchLength = siteModel.getRateForCategory(i) * branchTime;
                siteModel.getSubstitutionModel().getTransitionProbabilities(branchLength, probabilities);
				likelihoodCore.setNodeMatrix(nodeNum, i, probabilities);
			}

			update = true;
		}

		// If the node is internal, update the partial likelihoods.
		if (!tree.isExternal(node)) {

			int nodeCount = tree.getChildCount(node);
			if (nodeCount != 2)
				throw new RuntimeException("binary trees only!");

			// Traverse down the two child nodes
			NodeRef child1 = tree.getChild(node, 0);
			boolean update1 = traverse(tree, child1);

			NodeRef child2 = tree.getChild(node, 1);
			boolean update2 = traverse(tree, child2);

			// If either child node was updated then update this node too
			if (update1 || update2) {

				int childNum1 = child1.getNumber();
				int childNum2 = child2.getNumber();

				if (integrateAcrossCategories) {
					likelihoodCore.calculatePartials(childNum1, childNum2, nodeNum);
				} else {
					likelihoodCore.calculatePartials(childNum1, childNum2, nodeNum,
							siteCategories);
				}

				if (parent == null) {
					// No parent this is the root of the tree -
					// calculate the pattern likelihoods
					double[] frequencies = frequencyModel.getFrequencies();

					if (integrateAcrossCategories) {

						// moved this call to here, because non-integrating siteModels don't need to support it - AD
						double[] proportions = siteModel.getCategoryProportions();
						likelihoodCore.integratePartials(nodeNum, proportions, rootPartials);
					} else {
						likelihoodCore.getPartials(nodeNum, rootPartials);
					}

					likelihoodCore.calculateLogLikelihoods(rootPartials, frequencies, patternLogLikelihoods);
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
			return ARG_LIKELIHOOD;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			boolean useAmbiguities = false;
			boolean storePartials = true;
			boolean useScaling = false;
			if (xo.hasAttribute(USE_AMBIGUITIES)) {
				useAmbiguities = xo.getBooleanAttribute(USE_AMBIGUITIES);
			}
			if (xo.hasAttribute(STORE_PARTIALS)) {
				storePartials = xo.getBooleanAttribute(STORE_PARTIALS);
			}
			if (xo.hasAttribute(USE_SCALING)) {
				useScaling = xo.getBooleanAttribute(USE_SCALING);
			}
			PatternList patternList = (PatternList) xo.getChild(PatternList.class);
			ARGModel treeModel = (ARGModel) xo.getChild(ARGModel.class);
			SiteModel siteModel = (SiteModel) xo.getChild(SiteModel.class);

			BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

			return new ARGLikelihood(patternList, treeModel, siteModel, branchRateModel, useAmbiguities, storePartials, useScaling);
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
				AttributeRule.newBooleanRule(STORE_PARTIALS, true),
				AttributeRule.newBooleanRule(USE_SCALING, true),
				new ElementRule(PatternList.class),
				new ElementRule(ARGModel.class),
				new ElementRule(SiteModel.class),
				new ElementRule(BranchRateModel.class, true)
		};
	};

	/**
	 * XML Serializer for parallelization
	 *
	 */

//	public Element toXML() {
//		Element likelihoodElement
//	}

	// **************************************************************
	// INSTANCE VARIABLES
	// **************************************************************

	/**
	 * the frequency model for these sites
	 */
	protected FrequencyModel frequencyModel = null;

	/**
	 * the site model for these sites
	 */
	protected SiteModel siteModel = null;

	/**
	 * the branch rate model
	 */
	protected BranchRateModel branchRateModel = null;

	private boolean storePartials = false;

	private boolean integrateAcrossCategories = false;

	/**
	 * the categories for each site
	 */
	protected int[] siteCategories = null;

	/**
	 * the root partial likelihoods
	 */
	protected double[] rootPartials = null;

	/**
	 * the pattern likelihoods
	 */
	protected double[] patternLogLikelihoods = null;

	/**
	 * the number of rate categories
	 */
	protected int categoryCount;

	/**
	 * an array used to store transition probabilities
	 */
	protected double[] probabilities;

	/**
	 * the LikelihoodCore
	 */
	protected LikelihoodCore likelihoodCore;

	private boolean useAmbiguities;
}