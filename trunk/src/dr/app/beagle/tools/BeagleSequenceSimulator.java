package dr.app.beagle.tools;

import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.HKY;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.inference.markovjumps.StateHistory;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

@SuppressWarnings("serial")
public class BeagleSequenceSimulator extends SimpleAlignment {

	/** number of replications **/
	private int nReplications;
	/** tree used for generating samples **/
	private Tree tree;
	/** site model used for generating samples **/
	private GammaSiteRateModel siteModel;
	/** branch rate model used for generating samples **/
	private BranchRateModel branchRateModel;
	/** nr of categories in site model **/
	 int categoryCount;
	/** nr of states in site model **/
	 int stateCount;

	private boolean has_ancestralSequence = false;
	private Sequence ancestralSequence;

	/**
	 * Constructor
	 * 
	 * @param tree
	 * @param siteModel
	 * @param branchRateModel
	 * @param nReplications: nr of sites to generate
	 */
	public BeagleSequenceSimulator(Tree tree, GammaSiteRateModel siteModel,
			BranchRateModel branchRateModel, int nReplications) {
		this.tree = tree;
		this.siteModel = siteModel;
		this.branchRateModel = branchRateModel;
		this.nReplications = nReplications;
		stateCount = this.siteModel.getSubstitutionModel().getDataType().getStateCount();
		categoryCount = this.siteModel.getCategoryCount();
	} // END: Constructor

	/**
	 * Convert integer representation of sequence into a Sequence
	 * 
	 * @param seq integer representation of the sequence
	 * @param node used to determine taxon for sequence
	 * @return Sequence
	 */
	Sequence intArray2Sequence(int[] seq, NodeRef node) {
		String sSeq = "";
		DataType dataType = siteModel.getSubstitutionModel().getDataType();
		for (int i = 0; i < nReplications; i++) {
			if (dataType instanceof Codons) {
				String s = dataType.getTriplet(seq[i]);
				sSeq += s;
			} else {
				String c = dataType.getCode(seq[i]);
				sSeq += c;
			}
		}
		return new Sequence(tree.getNodeTaxon(node), sSeq);
	}// END: intArray2Sequence

	void setAncestralSequence(Sequence seq) {
		ancestralSequence = seq;
		has_ancestralSequence = true;
	}// END: setAncestralSequence

	int[] sequence2intArray(Sequence seq) {

		int array[] = new int[nReplications];

		if (seq.getLength() != nReplications) {

			throw new RuntimeException("Ancestral sequence length has "
					+ seq.getLength() + " characters " + "expecting "
					+ nReplications + " characters");

		} else {

			for (int i = 0; i < nReplications; i++) {
				array[i] = siteModel.getSubstitutionModel().getDataType()
						.getState(seq.getChar(i));
			}

		}

		return array;
	}// END: sequence2intArray

	/**
	 * perform the actual sequence generation
	 * 
	 * @return alignment containing randomly generated sequences for the nodes
	 *         in the leaves of the tree
	 */
	public void simulate() {

		double[] lambda = new double[stateCount * stateCount];

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

		}

		this.setReportCountStatistics(true);
		setDataType(siteModel.getSubstitutionModel().getDataType());

		traverse(root, seq, category, this, lambda);
	}//END: simulate

	/**
	 * recursively walk through the tree top down, and add sequence to alignment
	 * whenever a leave node is reached.
	 * 
	 * @param node
	 *            reference to the current node, for which we visit all children
	 * @param parentSequence
	 *            randomly generated sequence of the parent node
	 * @param category
	 *            array of categories for each of the sites
	 * @param alignment
	 */
	private void traverse(NodeRef node, int[] parentSequence, int[] category,
			SimpleAlignment alignment, double[] lambda) {

		for (int iChild = 0; iChild < tree.getChildCount(node); iChild++) {

			NodeRef child = tree.getChild(node, iChild);
			int[] seq = new int[nReplications];
			StateHistory[] histories = new StateHistory[nReplications];

			for (int i = 0; i < nReplications; i++) {
				histories[i] = simulateAlongBranch(tree, child, category[i],
						parentSequence[i], lambda);
				seq[i] = histories[i].getEndingState();
			}

			if (tree.getChildCount(child) == 0) {
				alignment.addSequence(intArray2Sequence(seq, child));
			}
			traverse(tree.getChild(node, iChild), seq, category, alignment,
					lambda);
		}
	}// END: traverse

	private StateHistory simulateAlongBranch(Tree tree, NodeRef node,
			int rateCategory, int startingState, double[] lambda) {

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

		return StateHistory.simulateUnconditionalOnEndingState(0.0,
				startingState, branchLength, lambda, stateCount);
	}// END: simulateAlongBranch

	public String toString() {
		StringBuffer sb = new StringBuffer();
		// alignment output
		sb.append("alignment\n");
		sb.append(super.toString());
		sb.append("\n");
		return sb.toString();
	}// END: toString

	/** generate simple site model, for testing purposes **/
	static GammaSiteRateModel getDefaultGammaSiteRateModel() {
		Parameter kappa = new Parameter.Default(1, 2);
		Parameter freqs = new Parameter.Default(new double[] { 0.25, 0.25,
				0.25, 0.25 });
		FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
		HKY hky = new HKY(kappa, f);
		// TODO
		return new GammaSiteRateModel(hky.getModelName());
	} // getDefaultSiteModel

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

			System.out.println(treeSimulator.toString());

		} catch (Exception e) {
			e.printStackTrace();
		}// END: try-catch block

	} // END: main

}// END: class
