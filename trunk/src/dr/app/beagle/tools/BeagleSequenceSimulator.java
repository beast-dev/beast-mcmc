package dr.app.beagle.tools;

import java.util.ArrayList;
import java.util.List;

import beagle.Beagle;
import beagle.BeagleFactory;
import dr.app.beagle.evomodel.sitemodel.BranchSubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.EpochBranchSubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.sitemodel.HomogenousBranchSubstitutionModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.HKY;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.treelikelihood.BufferIndexHelper;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

public class BeagleSequenceSimulator {

	private TreeModel treeModel;
	private GammaSiteRateModel gammaSiteRateModel;
	private BranchSubstitutionModel branchSubstitutionModel;
	private int sequenceLength;
	private FrequencyModel freqModel;
	private int categoryCount;
	private Beagle beagle;
	private BufferIndexHelper eigenBufferHelper;
	private BufferIndexHelper matrixBufferHelper;

	private int stateCount;

	private boolean has_ancestralSequence = false;
	private Sequence ancestralSequence = null;
	private double[][] probabilities;

	// TODO: get freqModel from branchSubstModel
	public BeagleSequenceSimulator(TreeModel treeModel, //
			GammaSiteRateModel gammaSiteRateModel, //
			BranchRateModel branchRateModel, //
			FrequencyModel freqModel, //
			int sequenceLength //
	) {

		this.treeModel = treeModel;
		this.gammaSiteRateModel = gammaSiteRateModel;
		this.sequenceLength = sequenceLength;
		this.freqModel = freqModel; // gammaSiteRateModel.getSubstitutionModel().getFrequencyModel();
		this.branchSubstitutionModel = (BranchSubstitutionModel) gammaSiteRateModel.getModel(0); // branchSubstitutionModel;

		int tipCount = treeModel.getExternalNodeCount();
		int nodeCount = treeModel.getNodeCount();
		int eigenCount = branchSubstitutionModel.getEigenCount();
		int internalNodeCount = treeModel.getInternalNodeCount();
		int scaleBufferCount = internalNodeCount + 1;

		int compactPartialsCount = tipCount;

		int patternCount = sequenceLength;

		int stateCount = freqModel.getDataType().getStateCount();

		this.categoryCount = gammaSiteRateModel.getCategoryCount();
		this.probabilities = new double[categoryCount][stateCount * stateCount];
		this.stateCount = stateCount;

//		print2DArray(probabilities);
		
		// one partials buffer for each tip and two for each internal node (for store restore)
		BufferIndexHelper partialBufferHelper = new BufferIndexHelper(nodeCount, tipCount);

		// two eigen buffers for each decomposition for store and restore.
		eigenBufferHelper = new BufferIndexHelper(eigenCount, 0);

		// two matrices for each node less the root
		matrixBufferHelper = new BufferIndexHelper(nodeCount, 0);

		// one scaling buffer for each internal node plus an extra for the
		// accumulation, then doubled for store/restore
		BufferIndexHelper scaleBufferHelper = new BufferIndexHelper(scaleBufferCount, 0);

		// null implies no restrictions
		int[] resourceList = null;
		long preferenceFlags = 0;
		long requirementFlags = 0;

		beagle = BeagleFactory.loadBeagleInstance(tipCount, //
				partialBufferHelper.getBufferCount(), //
				compactPartialsCount, //
				stateCount, //
				patternCount, //
				eigenBufferHelper.getBufferCount(), // 
				matrixBufferHelper.getBufferCount() + branchSubstitutionModel.getExtraBufferCount(treeModel), //
				categoryCount, //
				scaleBufferHelper.getBufferCount(), // 
				resourceList, //
				preferenceFlags, //
				requirementFlags //
				);

		double[] categoryRates = gammaSiteRateModel.getCategoryRates();
		beagle.setCategoryRates(categoryRates);
		
	}// END: Constructor

	public void setAncestralSequence(Sequence seq) {
		ancestralSequence = seq;
		has_ancestralSequence = true;
	}// END: setAncestralSequence

	int[] sequence2intArray(Sequence seq) {

		if (seq.getLength() != sequenceLength) {

			throw new RuntimeException("Ancestral sequence length has " + seq.getLength() + " characters " + "expecting " + sequenceLength + " characters");

		}

		int array[] = new int[sequenceLength];
		for (int i = 0; i < sequenceLength; i++) {
			array[i] = freqModel.getDataType().getState(seq.getChar(i));
		}

		return array;
	}// END: sequence2intArray

	Sequence intArray2Sequence(int[] seq, NodeRef node) {
		StringBuilder sSeq = new StringBuilder();
		for (int i = 0; i < sequenceLength; i++) {
			sSeq.append(freqModel.getDataType().getCode(seq[i]));
		}

		return new Sequence(treeModel.getNodeTaxon(node), sSeq.toString());
	} // END: intArray2Sequence

	public Alignment simulate() {

		SimpleAlignment alignment = new SimpleAlignment();
		NodeRef root = treeModel.getRoot();
		double[] categoryProbs = gammaSiteRateModel.getCategoryProportions();
		int[] category = new int[sequenceLength];

		for (int i = 0; i < sequenceLength; i++) {
			category[i] = MathUtils.randomChoicePDF(categoryProbs);
		}

//		printArray(categoryProbs);
//		 printArray(category);
		
		int[] seq = new int[sequenceLength];

		if (has_ancestralSequence) {

			seq = sequence2intArray(ancestralSequence);

		} else {

			for (int i = 0; i < sequenceLength; i++) {
				seq[i] = MathUtils.randomChoicePDF(freqModel.getFrequencies());

			}

		}// END: ancestral sequence check

		alignment.setDataType(freqModel.getDataType());
		alignment.setReportCountStatistics(true);

		traverse(root, seq, category, alignment);

		return alignment;
	}// END: simulate

	void traverse(NodeRef node, int[] parentSequence, int[] category,
			SimpleAlignment alignment) {

		for (int iChild = 0; iChild < treeModel.getChildCount(node); iChild++) {

			NodeRef child = treeModel.getChild(node, iChild);
			int[] sequence = new int[sequenceLength];
			double[] cProb = new double[stateCount];

			for (int i = 0; i < categoryCount; i++) {
				
				getTransitionProbabilities(treeModel, child, i, probabilities[i]);
				
			}

			for (int i = 0; i < sequenceLength; i++) {

				System.arraycopy(probabilities[category[i]], parentSequence[i] * stateCount, cProb, 0, stateCount);
				sequence[i] = MathUtils.randomChoicePDF(cProb);

			}

			if (treeModel.getChildCount(child) == 0) {
				alignment.addSequence(intArray2Sequence(sequence, child));
			}

			traverse(treeModel.getChild(node, iChild), sequence, category, alignment);

		}// END: child nodes loop
	}// END: traverse

	// TODO:
	void getTransitionProbabilities(Tree tree, NodeRef node, int rateCategory,
			double[] probabilities) {

		int nodeNum = node.getNumber();
		matrixBufferHelper.flipOffset(nodeNum);
		int bufferIndex = matrixBufferHelper.getOffsetIndex(nodeNum);
		int eigenIndex = branchSubstitutionModel.getBranchIndex(tree, node, bufferIndex);
		double edgeLengths[] = { tree.getBranchLength(node) };
		int probabilityIndices[] = { bufferIndex };
		int count = 1;

        eigenBufferHelper.flipOffset(eigenIndex);

        branchSubstitutionModel.setEigenDecomposition(
                beagle, eigenIndex, eigenBufferHelper, 0
        );

//		 System.out.println("eigenIndex:" + eigenIndex);
//		 System.out.println("update indices:");
//		 System.out.println(probabilityIndices[0]);
//		 System.out.println("weights:");
//		 System.out.println(edgeLengths[0]);
		
		branchSubstitutionModel.updateTransitionMatrices(beagle, //
				eigenIndex, //
				eigenBufferHelper, //
				probabilityIndices, //
				null, //
				null, //
				edgeLengths, //
				count //
				);

		beagle.getTransitionMatrix(bufferIndex, probabilities);

		System.out.println("bufferIndex: " + bufferIndex);
		printArray(probabilities);
		
	}// END: getTransitionProbabilities

	// /////////////////
	// ---DEBUGGING---//
	// /////////////////

	public static void main(String[] args) {

//		simulateEpochModel();
		simulateHKY();

	} // END: main

	static void simulateHKY() {

		try {

			int sequenceLength = 10;

			// create tree
			NewickImporter importer = new NewickImporter("(SimSeq1:73.7468,(SimSeq2:25.256989999999995,SimSeq3:45.256989999999995):18.48981);");
			Tree tree = importer.importTree(null);
			TreeModel treeModel = new TreeModel(tree);

			// create Frequency Model
			Parameter freqs = new Parameter.Default(new double[] { 0.25, 0.25, 0.25, 0.25 });
			FrequencyModel freqModel = new FrequencyModel(Nucleotides.INSTANCE, freqs);

			// create substitution model
			Parameter kappa = new Parameter.Default(1, 1);
			HKY hky = new HKY(kappa, freqModel);
			HomogenousBranchSubstitutionModel substitutionModel = new HomogenousBranchSubstitutionModel(hky, freqModel);
			
			// create site model
			GammaSiteRateModel siteRateModel = new GammaSiteRateModel("siteModel");
			siteRateModel.addModel(substitutionModel);
			
			// create branch rate model
			BranchRateModel branchRateModel = new DefaultBranchRateModel();

			// feed to sequence simulator and generate leaves
			BeagleSequenceSimulator beagleSequenceSimulator = new BeagleSequenceSimulator(
					treeModel, //
					siteRateModel, //
					branchRateModel, //
					freqModel, //
					sequenceLength //
			);

			Sequence ancestralSequence = new Sequence();
			ancestralSequence.appendSequenceString("AAAAAAAAAA");
			beagleSequenceSimulator.setAncestralSequence(ancestralSequence);

			System.out.println(beagleSequenceSimulator.simulate().toString());

		} catch (Exception e) {
			e.printStackTrace();
		}// END: try-catch block

	}// END: simulateHKY
	
	static void simulateEpochModel() {

		try {

			int sequenceLength = 10;

			// create tree
			NewickImporter importer = new NewickImporter("(SimSeq1:73.7468,(SimSeq2:25.256989999999995,SimSeq3:45.256989999999995):18.48981);");
			Tree tree = importer.importTree(null);
			TreeModel treeModel = new TreeModel(tree);

			// create Frequency Model
			Parameter freqs = new Parameter.Default(new double[] { 0.25, 0.25, 0.25, 0.25 });
			FrequencyModel freqModel = new FrequencyModel(Nucleotides.INSTANCE, freqs);
			List<FrequencyModel> frequencyModelList = new ArrayList<FrequencyModel>();
			frequencyModelList.add(freqModel);

			// create Epoch Model
			Parameter kappa1 = new Parameter.Default(1, 1);
			Parameter kappa2 = new Parameter.Default(1, 10);
			HKY hky1 = new HKY(kappa1, freqModel);
			HKY hky2 = new HKY(kappa2, freqModel);
			List<SubstitutionModel> substModelList = new ArrayList<SubstitutionModel>();
			substModelList.add(hky1);
			substModelList.add(hky2);

			Parameter epochTimes = new Parameter.Default(1, 20);
			EpochBranchSubstitutionModel substitutionModel = new EpochBranchSubstitutionModel(
					substModelList, //
					frequencyModelList, //
					epochTimes //
			);

			// create site model
			GammaSiteRateModel siteRateModel = new GammaSiteRateModel("siteModel");
			siteRateModel.addModel(substitutionModel);
			
			// create branch rate model
			BranchRateModel branchRateModel = new DefaultBranchRateModel();

			// feed to sequence simulator and generate leaves
			BeagleSequenceSimulator beagleSequenceSimulator = new BeagleSequenceSimulator(
					treeModel, //
					siteRateModel, //
					branchRateModel, //
					freqModel, //
					sequenceLength //
			);

			Sequence ancestralSequence = new Sequence();
			ancestralSequence.appendSequenceString("TCAGGTCAAG");
			beagleSequenceSimulator.setAncestralSequence(ancestralSequence);

			System.out.println(beagleSequenceSimulator.simulate().toString());

		} catch (Exception e) {
			e.printStackTrace();
		}// END: try-catch block

	}// END : simulateEpochModel
	
	public static void printArray(int[] category) {
		for (int i = 0; i < category.length; i++) {
			System.out.println(category[i]);
		}
	}// END: printArray

	public static void printArray(double[] category) {
		for (int i = 0; i < category.length; i++) {
			System.out.println(category[i]);
		}
		System.out.print("\n");
	}// END: printArray

	public void print2DArray(double[][] array) {
		for (int row = 0; row < array.length; row++) {
			for (int col = 0; col < array[row].length; col++) {
				System.out.print(array[row][col] + " ");
			}
			System.out.print("\n");
		}
	}// END: print2DArray

	public static void print2DArray(int[][] array) {
		for (int row = 0; row < array.length; row++) {
			for (int col = 0; col < array[row].length; col++) {
				System.out.print(array[row][col] + " ");
			}
			System.out.print("\n");
		}
	}// END: print2DArray

} // END: class
