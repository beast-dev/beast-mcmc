package dr.app.bss.test;

import java.util.ArrayList;

import dr.app.beagle.evomodel.branchmodel.HomogeneousBranchModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.HKY;
import dr.app.beagle.tools.BeagleSequenceSimulator;
import dr.app.beagle.tools.Partition;
import dr.app.bss.Utils;
import dr.evolution.coalescent.CoalescentSimulator;
import dr.evolution.coalescent.ExponentialGrowth;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

public class BeagleSeqSimTest {

	public static final boolean simulateInPar = true;
	
	public static void main(String[] args) {

//		simulateTopology();
//		 simulateOnePartition();
//		 simulateTwoPartitions();
		 int N = 10000;
		for (int i = 0; i < N; i++) {
			simulateThreePartitions();
		}
		
	} // END: main

	static void simulateTopology() {

		try {

			int sequenceLength = 10;
			ArrayList<Partition> partitionsList = new ArrayList<Partition>();

			MathUtils.setSeed(666);

			// create tree
			NewickImporter importer = new NewickImporter(
					"(SimSeq1:73.7468,(SimSeq2:25.256989999999995,SimSeq3:45.256989999999995):18.48981);");
			Tree tree = importer.importTree(null);

			// set demographic function
			ExponentialGrowth exponentialGrowth = new ExponentialGrowth(
					Units.Type.YEARS);
			exponentialGrowth.setN0(10);
			exponentialGrowth.setGrowthRate(0.5);

			Taxa taxa = new Taxa();
			for (Taxon taxon : tree.asList()) {

				double absoluteHeight = Utils.getAbsoluteTaxonHeight(taxon,
						tree);

				taxon.setAttribute(Utils.ABSOLUTE_HEIGHT, absoluteHeight);

				// taxon.setAttribute("date", new Date(absoluteHeight,
				// Units.Type.YEARS, true));

				taxa.addTaxon(taxon);

			}// END: taxon loop

			CoalescentSimulator topologySimulator = new CoalescentSimulator();

			TreeModel treeModel = new TreeModel(topologySimulator.simulateTree(
					taxa, exponentialGrowth));

			System.out.println(treeModel.toString());

			Parameter freqs = new Parameter.Default(new double[] { 0.25, 0.25,
					0.25, 0.25 });
			FrequencyModel freqModel = new FrequencyModel(Nucleotides.INSTANCE,
					freqs);

			// create substitution model
			Parameter kappa = new Parameter.Default(1, 10);
			HKY hky = new HKY(kappa, freqModel);
			HomogeneousBranchModel substitutionModel = new HomogeneousBranchModel(
					hky);

			// create site model
			GammaSiteRateModel siteRateModel = new GammaSiteRateModel(
					"siteModel");

			// create branch rate model
			BranchRateModel branchRateModel = new DefaultBranchRateModel();

			// create partition
			Partition partition1 = new Partition(treeModel, //
					substitutionModel,//
					siteRateModel, //
					branchRateModel, //
					freqModel, //
					0, // from
					sequenceLength - 1, // to
					1 // every
			);

			partitionsList.add(partition1);

			// feed to sequence simulator and generate data
			BeagleSequenceSimulator simulator = new BeagleSequenceSimulator(
					partitionsList
			);
			
			System.out.println(simulator.simulate(simulateInPar).toString());
			
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		} // END: try-catch block

	}// END: simulate topology

	static void simulateOnePartition() {

		try {

			MathUtils.setSeed(666);

			System.out.println("Test case 1: simulateOnePartition");

			int sequenceLength = 10;
			ArrayList<Partition> partitionsList = new ArrayList<Partition>();

			// create tree
			NewickImporter importer = new NewickImporter(
					"(SimSeq1:73.7468,(SimSeq2:25.256989999999995,SimSeq3:45.256989999999995):18.48981);");
			Tree tree = importer.importTree(null);
			TreeModel treeModel = new TreeModel(tree);

			// create Frequency Model
			Parameter freqs = new Parameter.Default(new double[] { 0.25, 0.25,
					0.25, 0.25 });
			FrequencyModel freqModel = new FrequencyModel(Nucleotides.INSTANCE,
					freqs);

			// create substitution model
			Parameter kappa = new Parameter.Default(1, 10);
			HKY hky = new HKY(kappa, freqModel);
			HomogeneousBranchModel substitutionModel = new HomogeneousBranchModel(
					hky);

			// create site model
			GammaSiteRateModel siteRateModel = new GammaSiteRateModel(
					"siteModel");

			// create branch rate model
			BranchRateModel branchRateModel = new DefaultBranchRateModel();

			// create partition
			Partition partition1 = new Partition(treeModel, //
					substitutionModel,//
					siteRateModel, //
					branchRateModel, //
					freqModel, //
					0, // from
					sequenceLength - 1, // to
					1 // every
			);

			Sequence ancestralSequence = new Sequence();
			ancestralSequence.appendSequenceString("TCAAGTGAGG");
			partition1.setAncestralSequence(ancestralSequence);

			partitionsList.add(partition1);

			// feed to sequence simulator and generate data
			BeagleSequenceSimulator simulator = new BeagleSequenceSimulator(
					partitionsList
			// , sequenceLength
			);

			System.out.println(simulator.simulate(simulateInPar).toString());
			
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		} // END: try-catch block

	}// END: simulateOnePartition

	static void simulateTwoPartitions() {

		try {

			System.out.println("Test case 2: simulateTwoPartitions");

			int sequenceLength = 11;
			ArrayList<Partition> partitionsList = new ArrayList<Partition>();

			// create tree
			NewickImporter importer = new NewickImporter(
					"(SimSeq1:73.7468,(SimSeq2:25.256989999999995,SimSeq3:45.256989999999995):18.48981);");
			Tree tree = importer.importTree(null);
			TreeModel treeModel = new TreeModel(tree);

			// create Frequency Model
			Parameter freqs = new Parameter.Default(new double[] { 0.25, 0.25,
					0.25, 0.25 });
			FrequencyModel freqModel = new FrequencyModel(Nucleotides.INSTANCE,
					freqs);

			// create substitution model
			Parameter kappa = new Parameter.Default(1, 10);
			HKY hky = new HKY(kappa, freqModel);
			HomogeneousBranchModel substitutionModel = new HomogeneousBranchModel(
					hky);

			// create site model
			GammaSiteRateModel siteRateModel = new GammaSiteRateModel(
					"siteModel");

			// create branch rate model
			BranchRateModel branchRateModel = new DefaultBranchRateModel();

			// create partition
			Partition partition1 = new Partition(treeModel, //
					substitutionModel, //
					siteRateModel, //
					branchRateModel, //
					freqModel, //
					0, // from
					3, // to
					1 // every
			);

			// create partition
			Partition partition2 = new Partition(treeModel, //
					substitutionModel,//
					siteRateModel, //
					branchRateModel, //
					freqModel, //
					4, // from
					sequenceLength - 1, // to
					1 // every
			);

			// Sequence ancestralSequence = new Sequence();
			// ancestralSequence.appendSequenceString("TCAAGTGAGG");
			// partition2.setAncestralSequence(ancestralSequence);

			partitionsList.add(partition1);
			partitionsList.add(partition2);

			// feed to sequence simulator and generate data
			BeagleSequenceSimulator simulator = new BeagleSequenceSimulator(
					partitionsList);
			
			System.out.println(simulator.simulate(simulateInPar).toString());
			
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		} // END: try-catch block

	}// END: simulateTwoPartitions

	static void simulateThreePartitions() {

		try {

			
			MathUtils.setSeed(666);
			
			System.out.println("Test case 3: simulateThreePartitions");

			int sequenceLength = 10;
			ArrayList<Partition2> partitionsList = new ArrayList<Partition2>();

			// create tree
			NewickImporter importer = new NewickImporter(
					"(SimSeq1:73.7468,(SimSeq2:25.256989999999995,SimSeq3:45.256989999999995):18.48981);");
			Tree tree = importer.importTree(null);
			TreeModel treeModel = new TreeModel(tree);

			// create Frequency Model
			Parameter freqs = new Parameter.Default(new double[] { 0.25, 0.25,
					0.25, 0.25 });
			FrequencyModel freqModel = new FrequencyModel(Nucleotides.INSTANCE,
					freqs);

			// create substitution model
			Parameter kappa = new Parameter.Default(1, 10);
			HKY hky = new HKY(kappa, freqModel);
			HomogeneousBranchModel substitutionModel = new HomogeneousBranchModel(
					hky);

			// create site model
			GammaSiteRateModel siteRateModel = new GammaSiteRateModel(
					"siteModel");

			// create branch rate model
			BranchRateModel branchRateModel = new DefaultBranchRateModel();

			// create partition
			Partition2 partition1 = new Partition2(treeModel, //
					substitutionModel, //
					siteRateModel, //
					branchRateModel, //
					freqModel, //
					0, // from
					sequenceLength - 1, // to
					3 // every
			);

			// create partition
			Partition2 partition2 = new Partition2(treeModel, //
					substitutionModel,//
					siteRateModel, //
					branchRateModel, //
					freqModel, //
					1, // from
					sequenceLength - 1, // to
					3 // every
			);

			// create partition
			Partition2 partition3 = new Partition2(treeModel, //
					substitutionModel,//
					siteRateModel, //
					branchRateModel, //
					freqModel, //
					2, // from
					sequenceLength - 1, // to
					3 // every
			);

			partitionsList.add(partition1);
			partitionsList.add(partition2);
			partitionsList.add(partition3);

			// feed to sequence simulator and generate data
			BeagleSequenceSimulator2 simulator = new BeagleSequenceSimulator2(
					partitionsList);
			
			System.out.println(simulator.simulate(simulateInPar).toString());
			
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		} // END: try-catch block

	}// END: simulateThreePartitions

//	static void simulateThreePartitions() {
//
//		try {
//
//			
//			MathUtils.setSeed(666);
//			
//			System.out.println("Test case 3: simulateThreePartitions");
//
//			int sequenceLength = 10;
//			ArrayList<Partition> partitionsList = new ArrayList<Partition>();
//
//			// create tree
//			NewickImporter importer = new NewickImporter(
//					"(SimSeq1:73.7468,(SimSeq2:25.256989999999995,SimSeq3:45.256989999999995):18.48981);");
//			Tree tree = importer.importTree(null);
//			TreeModel treeModel = new TreeModel(tree);
//
//			// create Frequency Model
//			Parameter freqs = new Parameter.Default(new double[] { 0.25, 0.25,
//					0.25, 0.25 });
//			FrequencyModel freqModel = new FrequencyModel(Nucleotides.INSTANCE,
//					freqs);
//
//			// create substitution model
//			Parameter kappa = new Parameter.Default(1, 10);
//			HKY hky = new HKY(kappa, freqModel);
//			HomogeneousBranchModel substitutionModel = new HomogeneousBranchModel(
//					hky);
//
//			// create site model
//			GammaSiteRateModel siteRateModel = new GammaSiteRateModel(
//					"siteModel");
//
//			// create branch rate model
//			BranchRateModel branchRateModel = new DefaultBranchRateModel();
//
//			// create partition
//			Partition partition1 = new Partition(treeModel, //
//					substitutionModel, //
//					siteRateModel, //
//					branchRateModel, //
//					freqModel, //
//					0, // from
//					sequenceLength - 1, // to
//					3 // every
//			);
//
//			// create partition
//			Partition partition2 = new Partition(treeModel, //
//					substitutionModel,//
//					siteRateModel, //
//					branchRateModel, //
//					freqModel, //
//					1, // from
//					sequenceLength - 1, // to
//					3 // every
//			);
//
//			// create partition
//			Partition partition3 = new Partition(treeModel, //
//					substitutionModel,//
//					siteRateModel, //
//					branchRateModel, //
//					freqModel, //
//					2, // from
//					sequenceLength - 1, // to
//					3 // every
//			);
//
//			partitionsList.add(partition1);
//			partitionsList.add(partition2);
//			partitionsList.add(partition3);
//
//			// feed to sequence simulator and generate data
//			BeagleSequenceSimulator simulator = new BeagleSequenceSimulator(
//					partitionsList);
//			
//			System.out.println(simulator.simulate(simulateInPar).toString());
//			
//		} catch (Exception e) {
//			e.printStackTrace();
//			System.exit(-1);
//		} // END: try-catch block
//
//	}// END: simulateThreePartitions
	
}// END: class
