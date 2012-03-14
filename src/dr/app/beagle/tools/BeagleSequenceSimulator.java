package dr.app.beagle.tools;

import java.util.ArrayList;
import java.util.List;

import dr.app.beagle.evomodel.sitemodel.EpochBranchSubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.HKY;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

/**
 * Class for performing random sequence generation for a given site model.
 * Sequences for the leave nodes in the tree are returned as an alignment.
 */

public class BeagleSequenceSimulator {

	/** nr of samples to generate **/
	private int nReplications;
	/** tree used for generating samples **/
	private Tree tree;
	/** site model used for generating samples **/
	private GammaSiteRateModel siteModel;
	/** branch rate model used for generating samples **/
	private BranchRateModel branchRateModel;
	/** nr of categories in site model **/
	private int categoryCount;
	/** nr of states in site model **/
	private int stateCount;

	private static boolean has_ancestralSequence = false;
	private Sequence ancestralSequence;

	/** an array used to transfer transition probabilities **/
	private double[][] m_probabilities;

	/**
	 * Constructor
	 * 
	 * @param tree
	 * @param siteModel
	 * @param branchRateModel
	 * @param sequenceLength: nr of sites to generate
	 */
	BeagleSequenceSimulator(Tree tree, GammaSiteRateModel siteModel,
			BranchRateModel branchRateModel, int sequenceLength) {
		this.tree = tree;
		this.siteModel = siteModel;
		this.branchRateModel = branchRateModel;
		this.nReplications = sequenceLength;
		this.stateCount = siteModel.getSubstitutionModel().getDataType().getStateCount();
		this.categoryCount = siteModel.getCategoryCount();
		this.m_probabilities = new double[categoryCount][stateCount * stateCount];
	} // END: Constructor

	/**
	 * Convert integer representation of sequence into a Sequence
	 * 
	 * @param seq: integer representation of the sequence
	 * @param node: used to determine taxon for sequence
	 * @return Sequence
	 */
	Sequence intArray2Sequence(int[] seq, NodeRef node) {
		
		StringBuilder sSeq = new StringBuilder();
		
		for (int i = 0; i < nReplications; i++) {
			sSeq.append(siteModel.getSubstitutionModel().getDataType().getCode(seq[i]));
		}
		
		return new Sequence(tree.getNodeTaxon(node), sSeq.toString());
	} // intArray2Sequence

	void setAncestralSequence(Sequence seq) {
		ancestralSequence = seq;
		has_ancestralSequence = true;
	}

	int[] sequence2intArray(Sequence seq) {

		if (seq.getLength() != nReplications) {
			throw new RuntimeException("Ancestral sequence length has "
					+ seq.getLength() + " characters " + "expecting "
					+ nReplications + " characters");
		}

		int array[] = new int[nReplications];
		for (int i = 0; i < nReplications; i++) {
			array[i] = siteModel.getSubstitutionModel().getDataType().getState(
					seq.getChar(i));
		}
		
		return array;
	}// END: sequence2intArray

	/**
	 * perform the actual sequence generation
	 * 
	 * @return alignment containing randomly generated sequences for the nodes
	 *         in the leaves of the tree
	 */
	public Alignment simulate() {

		NodeRef root = tree.getRoot();

		double[] categoryProbs = siteModel.getCategoryProportions();
		int[] category = new int[nReplications];
		for (int i = 0; i < nReplications; i++) {
			category[i] = MathUtils.randomChoicePDF(categoryProbs);
		}

		int[] seq = new int[nReplications];

		if (has_ancestralSequence) {

			seq = sequence2intArray(ancestralSequence);

		} else {

			FrequencyModel frequencyModel = siteModel.getSubstitutionModel()
					.getFrequencyModel();
			for (int i = 0; i < nReplications; i++) {
				seq[i] = MathUtils.randomChoicePDF(frequencyModel
						.getFrequencies());
			}

		}// END: ancestral sequence check

		SimpleAlignment alignment = new SimpleAlignment();
		alignment.setReportCountStatistics(true);
		alignment.setDataType(siteModel.getSubstitutionModel().getDataType());

		traverse(root, seq, category, alignment);

		return alignment;
	} // END: simulate

	/**
	 * recursively walk through the tree top down, and add sequence to alignment
	 * whenever a leave node is reached.
	 * 
	 * @param node: reference to the current node, for which we visit all children
	 * @param parentSequence: randomly generated sequence of the parent node
	 * @param category: array of categories for each of the sites
	 * @param alignment
	 */
	void traverse(NodeRef node, int[] parentSequence, int[] category,
			SimpleAlignment alignment) {

		for (int iChild = 0; iChild < tree.getChildCount(node); iChild++) {

			NodeRef child = tree.getChild(node, iChild);
			for (int i = 0; i < categoryCount; i++) {
				getTransitionProbabilities(tree, child, i, m_probabilities[i]);
			}

			int[] seq = new int[nReplications];
			double[] cProb = new double[stateCount];
			for (int i = 0; i < nReplications; i++) {
				System.arraycopy(m_probabilities[category[i]],
						parentSequence[i] * stateCount, cProb, 0, stateCount);
				seq[i] = MathUtils.randomChoicePDF(cProb);
			}

			if (tree.getChildCount(child) == 0) {
				alignment.addSequence(intArray2Sequence(seq, child));
			}
			traverse(tree.getChild(node, iChild), seq, category, alignment);
		}
	} // END: traverse

	void getTransitionProbabilities(Tree tree, NodeRef node, int rateCategory,
			double[] probs) {

		NodeRef parent = tree.getParent(node);

		final double branchRate = branchRateModel.getBranchRate(tree, node);

		// Get the operational time of the branch
		final double branchTime = branchRate
				* (tree.getNodeHeight(parent) - tree.getNodeHeight(node));

		if (branchTime < 0.0) {
			throw new RuntimeException("Negative branch length: " + branchTime);
		}

		double branchLength = siteModel.getRateForCategory(rateCategory)
				* branchTime;

//		if (siteModel.getSubstitutionModel() instanceof SubstitutionEpochModel) {
//			((SubstitutionEpochModel) siteModel.getSubstitutionModel())
//					.getTransitionProbabilities(tree.getNodeHeight(node), tree
//							.getNodeHeight(parent), branchLength, probs);
//			return;
//		}

		siteModel.getSubstitutionModel().getTransitionProbabilities(
				branchLength, probs);
	} // END: getTransitionProbabilities

	/** generate simple site model, for testing purposes **/
	static GammaSiteRateModel getDefaultGammaSiteRateModel() {
		
		List<FrequencyModel> frequencyModelList = new ArrayList<FrequencyModel>();
		List<SubstitutionModel> substModelList = new ArrayList<SubstitutionModel>();
		Parameter epochTransitionTimes = new Parameter.Default(1, 20);
		
		Parameter freqs = new Parameter.Default(new double[] { 0.25, 0.25,
				0.25, 0.25 });
		FrequencyModel freqModel = new FrequencyModel(Nucleotides.INSTANCE, freqs);
		
		Parameter kappa1 = new Parameter.Default(1, 1);
		Parameter kappa2 = new Parameter.Default(10, 1);
		
		HKY hky1 = new HKY(kappa1, freqModel);
		HKY hky2 = new HKY(kappa2, freqModel);
		
		substModelList.add((SubstitutionModel) hky1);
		substModelList.add((SubstitutionModel) hky2);
		
		frequencyModelList.add(freqModel);
		
		EpochBranchSubstitutionModel beagleSubstitutionEpochModel =	new EpochBranchSubstitutionModel(substModelList, frequencyModelList, epochTransitionTimes);
		
		//TODO: pass epoch model to site model
		GammaSiteRateModel gsrm = new GammaSiteRateModel("dupa");
		gsrm.setSubstitutionModel(hky1);

		return gsrm;
	} // END: getDefaultGammaSiteRateModel

	public static void main(String[] args) {

		try {

			int nReplications = 10;

			// create tree
			NewickImporter importer = new NewickImporter(
					"((A:1.0,B:1.0)AB:1.0,(C:1.0,D:1.0)CD:1.0)ABCD;");
			Tree tree = importer.importTree(null);

			// create site model
			GammaSiteRateModel siteModel = getDefaultGammaSiteRateModel();

			// create branch rate model
			BranchRateModel branchRateModel = new DefaultBranchRateModel();

			// feed to sequence simulator and generate leaves
			BeagleSequenceSimulator treeSimulator = new BeagleSequenceSimulator(
					tree, siteModel, branchRateModel, nReplications);

			Sequence ancestralSequence = new Sequence();
			ancestralSequence.appendSequenceString("TCAGGTCAAG");
			treeSimulator.setAncestralSequence(ancestralSequence);

			System.out.println(treeSimulator.simulate().toString());

		} catch (Exception e) {
			e.printStackTrace();
		}// END: try-catch block

	} // END: main

} // class SequenceSimulator
