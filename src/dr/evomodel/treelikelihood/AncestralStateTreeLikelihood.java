package dr.evomodel.treelikelihood;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.NodeAttributeProvider;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class AncestralStateTreeLikelihood extends TreeLikelihood implements NodeAttributeProvider {

	public static final String RECONSTRUCTING_TREE_LIKELIHOOD = "ancestralTreeLikelihood";
	public static final String RECONSTRUCTION_TAG = "state";
	public static final String TAG_NAME = "tagName";

	private DataType dataType;
	private int[][] reconstructedStates;
	private int redrawTime = 0;
	private String tag;


	/**
	 * Constructor.
	 * Now also takes a DataType so that ancestral states are printed using data codes
	 *
	 * @param patternList     -
	 * @param treeModel       -
	 * @param siteModel       -
	 * @param branchRateModel -
	 * @param useAmbiguities  -
	 * @param useScaling      -
	 * @param storePartials   -
	 * @param dataType        - need to provide the data-type, so that corrent data characters can be returned
	 * @param tag             - string label for reconstruction characters in tree log
	 */
	public AncestralStateTreeLikelihood(PatternList patternList, TreeModel treeModel,
	                                    SiteModel siteModel, BranchRateModel branchRateModel,
	                                    boolean useAmbiguities, boolean storePartials,
	                                    boolean useScaling,
	                                    DataType dataType,
	                                    String tag) {
		super(patternList, treeModel, siteModel, branchRateModel, useAmbiguities, storePartials, useScaling);
		this.dataType = dataType;
		this.tag = tag;

	}


	public String getNodeAttributeLabel() {
		return tag;
	}

	public String getAttributeForNode(Tree tree, NodeRef node) {

		TreeModel treeModel = (TreeModel) tree;

		if (redrawTime == 0) {
			traverseSample(treeModel, tree.getRoot(), null);
		}

		// Function gets called once for each node to log tree
		// After one log, prepare to redraw states
		redrawTime++;
		if (redrawTime == tree.getNodeCount())
			redrawTime = 0;

		return formattedState(reconstructedStates[node.getNumber()], dataType);


	}

	private static String formattedState(int[] state, DataType dataType) {
		StringBuffer sb = new StringBuffer();
		sb.append("\"");
		for (int i : state) {
			sb.append(dataType.getChar(i));
		}
		sb.append("\"");
		return sb.toString();
	}


	/**
	 * Traverse (pre-order) the tree sampling the internal node states.
	 *
	 * @param tree        - TreeModel on which to perform sampling
	 * @param node        - current node
	 * @param parentState - character state of the parent node to 'node'
	 */
	public void traverseSample(TreeModel tree, NodeRef node, int[] parentState) {

		if (reconstructedStates == null)
			reconstructedStates = new int[tree.getNodeCount()][patternCount];

		int nodeNum = node.getNumber();

		NodeRef parent = tree.getParent(node);

		// This function assumes that all partial likelihoods have already been calculated
		// This function also assumes that # of patterns = 1

		// If the node is internal, then sample its state given the state of its parent (pre-order traversal).

		double[] conditionalProbabilities = new double[stateCount];
		int[] state = new int[patternCount];

		if (!tree.isExternal(node)) {

			if (parent == null) {

				// This is the root node
				for (int j = 0; j < patternCount; j++) {
					System.arraycopy(rootPartials, j * stateCount, conditionalProbabilities, 0, stateCount);
					state[j] = MathUtils.randomChoicePDF(conditionalProbabilities);
					reconstructedStates[nodeNum][j] = state[j];
				}

			} else {

				// This is an internal node, but not the root
				double[] partialLikelihood = new double[stateCount * patternCount];
				likelihoodCore.getPartials(nodeNum, partialLikelihood);

//				final double branchRate = branchRateModel.getBranchRate(tree, node);
//
//				            // Get the operational time of the branch
//				final double branchTime = branchRate * ( tree.getNodeHeight(parent) - tree.getNodeHeight(node) );
//
//				for (int i = 0; i < categoryCount; i++) {
//
//				                siteModel.getTransitionProbabilitiesForCategory(i, branchTime, probabilities);
//
//				}
//


				if (categoryCount > 1)
					throw new RuntimeException("Reconstruction not implemented for multiple categories yet.");

				((AbstractLikelihoodCore) likelihoodCore).getNodeMatrix(nodeNum, 0, probabilities);


				for (int j = 0; j < patternCount; j++) {

					int parentIndex = parentState[j] * stateCount;

					for (int i = 0; i < stateCount; i++)
						conditionalProbabilities[i] = partialLikelihood[i] * probabilities[parentIndex + i];

					state[j] = MathUtils.randomChoicePDF(conditionalProbabilities);
					reconstructedStates[nodeNum][j] = state[j];

				}
			}

//			int nodeCount = tree.getChildCount(node);

			// Traverse down the two child nodes
			NodeRef child1 = tree.getChild(node, 0);
			traverseSample(tree, child1, state);

			NodeRef child2 = tree.getChild(node, 1);
			traverseSample(tree, child2, state);
		} else {

			// This is an external leaf
			//patternList.
			((AbstractLikelihoodCore) likelihoodCore).getNodeStates(nodeNum, reconstructedStates[nodeNum]);

		}
	}

	/**
	 * The XML parser
	 */
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return RECONSTRUCTING_TREE_LIKELIHOOD;
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
			TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
			SiteModel siteModel = (SiteModel) xo.getChild(SiteModel.class);

			BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

			DataType dataType = ((SubstitutionModel) xo.getChild(SubstitutionModel.class)).getDataType();

			String tag = RECONSTRUCTION_TAG;
			if (xo.hasAttribute(TAG_NAME))
				tag = xo.getStringAttribute(TAG_NAME);

			return new AncestralStateTreeLikelihood(patternList, treeModel, siteModel,
					branchRateModel, useAmbiguities, storePartials, useScaling, dataType, tag);
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
				AttributeRule.newStringRule(TAG_NAME, true),
				new ElementRule(PatternList.class),
				new ElementRule(TreeModel.class),
				new ElementRule(SiteModel.class),
				new ElementRule(BranchRateModel.class, true),
				new ElementRule(SubstitutionModel.class)
		};
	};

}
