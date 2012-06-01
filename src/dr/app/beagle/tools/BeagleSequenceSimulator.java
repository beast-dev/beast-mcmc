package dr.app.beagle.tools;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

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
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.substmodel.YangCodonModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

public class BeagleSequenceSimulator {

	private TreeModel treeModel;
	private GammaSiteRateModel siteModel;
	private BranchSubstitutionModel branchSubstitutionModel;
	private int replications;
	private FrequencyModel freqModel;
	private DataType dataType;
	private int categoryCount;
	private int eigenCount;
	private Beagle beagle;
	private BufferIndexHelper eigenBufferHelper;
	private BufferIndexHelper matrixBufferHelper;

	private int stateCount;

	private boolean has_ancestralSequence = false;
	private Sequence ancestralSequence = null;
	private double[][] probabilities;

	// TODO: simplify parser
	public BeagleSequenceSimulator(TreeModel treeModel, //
			BranchSubstitutionModel branchSubstitutionModel,
			GammaSiteRateModel siteModel, //
			BranchRateModel branchRateModel, //
			FrequencyModel freqModel, //
			int replications //
	) {

		this.treeModel = treeModel;
		this.siteModel = siteModel;
		this.replications = replications;
		this.freqModel = freqModel;
		this.dataType = freqModel.getDataType();
//		this.freqModel = siteModel.getSubstitutionModel().getFrequencyModel();
		this.branchSubstitutionModel = branchSubstitutionModel;
//		this.branchSubstitutionModel = (BranchSubstitutionModel) siteModel.getModel(0);
		this.eigenCount = branchSubstitutionModel.getEigenCount();
		
		int tipCount = treeModel.getExternalNodeCount();
		int nodeCount = treeModel.getNodeCount();
		int internalNodeCount = treeModel.getInternalNodeCount();
		int scaleBufferCount = internalNodeCount + 1;

		int compactPartialsCount = tipCount;

		int patternCount = replications;

		
		int stateCount = dataType.getStateCount();

		this.categoryCount = siteModel.getCategoryCount();
		this.probabilities = new double[categoryCount][stateCount * stateCount];
		this.stateCount = stateCount;

		// one partials buffer for each tip and two for each internal node (for store restore)
		BufferIndexHelper partialBufferHelper = new BufferIndexHelper(nodeCount, tipCount);

		// two eigen buffers for each decomposition for store and restore.
		eigenBufferHelper = new BufferIndexHelper(eigenCount, 0);

		// two matrices for each node less the root
		matrixBufferHelper = new BufferIndexHelper(nodeCount, 0);

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
				matrixBufferHelper.getBufferCount() + this.branchSubstitutionModel.getExtraBufferCount(treeModel), //
				categoryCount, //
				scaleBufferCount, //
				resourceList, //
				preferenceFlags, //
				requirementFlags //
				);

	}// END: Constructor

	public void setAncestralSequence(Sequence seq) {
		ancestralSequence = seq;
		has_ancestralSequence = true;
	}// END: setAncestralSequence

	// TODO: get triplets
	private int[] sequence2intArray(Sequence seq) {

		if (seq.getLength() != replications) {

			throw new RuntimeException("Ancestral sequence length has "
					+ seq.getLength() + " characters " + "expecting "
					+ replications + " characters");

		}

		int array[] = new int[replications];

		if (dataType instanceof Codons) {
			
//			for (int i = 0; i < replications; i++) {
//				array[i] = dataType.getState(seq.getChar(i));
//			}// END: replications loop
			
		} else {

			for (int i = 0; i < replications; i++) {
				array[i] = dataType.getState(seq.getChar(i));
			}// END: replications loop

		}// END: dataType check

		return array;
	}// END: sequence2intArray

	// TODO
	private Sequence intArray2Sequence(int[] seq, NodeRef node) {

		StringBuilder sSeq = new StringBuilder();

		if (dataType instanceof Codons) {

			for (int i = 0; i < replications; i++) {

				sSeq.append(dataType.getTriplet(seq[i]));

				System.err.println(sSeq);
				
				
			}// END: replications loop

//			System.exit(-1);
			
		} else {

			for (int i = 0; i < replications; i++) {

				sSeq.append(freqModel.getDataType().getCode(seq[i]));

			}// END: replications loop

		}// END: dataType check

		return new Sequence(treeModel.getNodeTaxon(node), sSeq.toString());
	}// END: intArray2Sequence

	public Alignment simulate() {

		SimpleAlignment alignment = new SimpleAlignment();
		NodeRef root = treeModel.getRoot();
		double[] categoryProbs = siteModel.getCategoryProportions();
		int[] category = new int[replications];

		for (int i = 0; i < replications; i++) {
			category[i] = MathUtils.randomChoicePDF(categoryProbs);
		}

		int[] seq = new int[replications];

		if (has_ancestralSequence) {

			seq = sequence2intArray(ancestralSequence);

		} else {

			for (int i = 0; i < replications; i++) {
				seq[i] = MathUtils.randomChoicePDF(freqModel.getFrequencies());

			}

		}// END: ancestral sequence check

		alignment.setDataType(dataType);
		alignment.setReportCountStatistics(false);

        for (int i = 0; i < eigenCount; i++) {
            eigenBufferHelper.flipOffset(i);

            branchSubstitutionModel.setEigenDecomposition(
                    beagle, i, eigenBufferHelper, 0
            );

    }
		
		double[] categoryRates = siteModel.getCategoryRates();
		beagle.setCategoryRates(categoryRates);   
//	    double[] categoryWeights = gammaSiteRateModel.getCategoryProportions();
//	    beagle.setCategoryWeights(0, categoryWeights);
//        double[] frequencies = branchSubstitutionModel.getStateFrequencies(0);
//        beagle.setStateFrequencies(0, frequencies);
		
		traverse(root, seq, category, alignment);

		return alignment;
	}// END: simulate

	private void traverse(NodeRef node, int[] parentSequence, int[] category,
			SimpleAlignment alignment) {

		for (int iChild = 0; iChild < treeModel.getChildCount(node); iChild++) {

			NodeRef child = treeModel.getChild(node, iChild);
			int[] sequence = new int[replications];
			double[] cProb = new double[stateCount];

			getTransitionProbabilities(treeModel, child, probabilities);
			
			for (int i = 0; i < replications; i++) {

				System.arraycopy(probabilities[category[i]], parentSequence[i] * stateCount, cProb, 0, stateCount);
				sequence[i] = MathUtils.randomChoicePDF(cProb);

			}

			if (treeModel.getChildCount(child) == 0) {
				alignment.addSequence(intArray2Sequence(sequence, child));
			}

			traverse(treeModel.getChild(node, iChild), sequence, category, alignment);

		}// END: child nodes loop
	}// END: traverse

	// TODO
	private void getTransitionProbabilities(Tree tree, NodeRef node,
			double[][] probabilities) {

		int nodeNum = node.getNumber();
		matrixBufferHelper.flipOffset(nodeNum);
		int branchIndex = matrixBufferHelper.getOffsetIndex(nodeNum);
		int eigenIndex = branchSubstitutionModel.getBranchIndex(tree, node, branchIndex);
		int count = 1;

		branchSubstitutionModel.updateTransitionMatrices(beagle, //
				eigenIndex, //
				eigenBufferHelper, //
				new int[] { branchIndex }, //
				null, //
				null, //
				new double[] { tree.getBranchLength(node) }, //
				count //
				);

		double transitionMatrix[] = new double[categoryCount * stateCount * stateCount];
		
		beagle.getTransitionMatrix(branchIndex, //
				transitionMatrix //
				);

		for (int i = 0; i < categoryCount; i++) {

			System.arraycopy(transitionMatrix, i * stateCount, probabilities[i], 0, stateCount * stateCount);
		
		}

//		System.out.println("eigenIndex:" + eigenIndex);
//		System.out.println("bufferIndex: " + branchIndex);
//		System.out.println("weight: " + tree.getBranchLength(node));
//		print2DArray(probabilities);
		
	}// END: getTransitionProbabilities
	
	// /////////////////
	// ---DEBUGGING---//
	// /////////////////

	public static void main(String[] args) {

//		simulateEpochModel();
//		simulateHKY();
		simulateCodon();
		
	} // END: main

	static void simulateCodon() {

		try {

			int sequenceLength = 10;

			// create tree
			NewickImporter importer = new NewickImporter("(SimSeq1:73.7468,(SimSeq2:25.256989999999995,SimSeq3:45.256989999999995):18.48981);");
			Tree tree = importer.importTree(null);
			TreeModel treeModel = new TreeModel(tree);

			// create Frequency Model
			Parameter freqs = new Parameter.Default(new double[] { 0.0163936, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 
					0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 
					0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 
					0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 
					0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344 });
			
			Codons codonDataType = Codons.UNIVERSAL;
			
			 dr.evomodel.substmodel.FrequencyModel freqModel = new  dr.evomodel.substmodel.FrequencyModel(codonDataType, freqs);

			// create codon substitution model
			Parameter kappa = new Parameter.Default(1, 10);
			Parameter omega = new Parameter.Default(1, 10);
			
			YangCodonModel yangCodonModel = new YangCodonModel(codonDataType, omega, kappa, freqModel);
			
//			HKY hky = new HKY(kappa, freqModel);
//			HomogenousBranchSubstitutionModel substitutionModel = new HomogenousBranchSubstitutionModel(hky, freqModel);
			
			// create site model
//			GammaSiteRateModel siteRateModel = new GammaSiteRateModel("siteModel");
//			siteRateModel.addModel(substitutionModel);
			
			// create branch rate model
//			BranchRateModel branchRateModel = new DefaultBranchRateModel();

			// feed to sequence simulator and generate leaves
//			BeagleSequenceSimulator beagleSequenceSimulator = new BeagleSequenceSimulator(
//					treeModel, //
//					substitutionModel,//
//					siteRateModel, //
//					branchRateModel, //
//					freqModel, //
//					sequenceLength //
//			);
//
//			Sequence ancestralSequence = new Sequence();
//			ancestralSequence.appendSequenceString("AAAAAAAAAA");
//			beagleSequenceSimulator.setAncestralSequence(ancestralSequence);
//
//			System.out.println(beagleSequenceSimulator.simulate().toString());

		} catch (Exception e) {
			e.printStackTrace();
		}// END: try-catch block

	}// END: simulateCodon
	
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
			Parameter kappa = new Parameter.Default(1, 10);
			HKY hky = new HKY(kappa, freqModel);
			HomogenousBranchSubstitutionModel substitutionModel = new HomogenousBranchSubstitutionModel(hky, freqModel);
			
			// create site model
			GammaSiteRateModel siteRateModel = new GammaSiteRateModel("siteModel");
//			siteRateModel.addModel(substitutionModel);
			
			// create branch rate model
			BranchRateModel branchRateModel = new DefaultBranchRateModel();

			// feed to sequence simulator and generate leaves
			BeagleSequenceSimulator beagleSequenceSimulator = new BeagleSequenceSimulator(
					treeModel, //
					substitutionModel,//
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
			Parameter kappa1 = new Parameter.Default(1, 10);
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
					substitutionModel,//
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

	public static void printArray(double[] matrix) {
		for (int i = 0; i < matrix.length; i++) {
//			System.out.println(matrix[i]);
			System.out.println(String.format(Locale.US, "%.20f", matrix[i]));
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
