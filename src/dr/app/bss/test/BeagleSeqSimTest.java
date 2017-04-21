/*
 * BeagleSeqSimTest.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.bss.test;

import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Locale;

import dr.evolution.tree.TreeUtils;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchmodel.RandomBranchModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.EmpiricalRateMatrix;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.aminoacid.Blosum62;
import dr.evomodel.substmodel.aminoacid.EmpiricalAminoAcidModel;
import dr.evomodel.substmodel.codon.GY94CodonModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.substmodel.codon.MG94CodonModel;
import dr.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.app.beagle.tools.BeagleSequenceSimulator;
import dr.app.beagle.tools.Partition;
import dr.app.bss.Utils;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.ConvertAlignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.coalescent.CoalescentSimulator;
import dr.evolution.coalescent.ExponentialGrowth;
import dr.evolution.datatype.AminoAcids;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneticCode;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.io.NewickImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTraitProvider;
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

		// start timing
		long tic = System.currentTimeMillis();
		int N = 1;
		for (int i = 0; i < N; i++) {
			// simulateTopology();
//			 simulateOnePartition();
			// simulateTwoPartitions();
			// simulateThreePartitions(i, N);
			// simulateAminoAcid();
			// simulateCodon();
//			ancestralSequenceTree();
			simulateRandomBranchAssignment();

		}
		long toc = System.currentTimeMillis();

		long time = toc - tic;
		System.out.println("Time: " + time + "milliseconds");

	} // END: main

	static void simulateRandomBranchAssignment() {
		
		try {
			
			System.out.println("Test case I dunno which: simulate random branch assignments");

			MathUtils.setSeed(666);

			int sequenceLength = 10;
			ArrayList<Partition> partitionsList = new ArrayList<Partition>();
			
	        File treeFile = new File("/home/filip/Dropbox/BeagleSequenceSimulator/SimTree/SimTree.figtree");
	        Tree tree = Utils.importTreeFromFile(treeFile);
	        TreeModel treeModel = new TreeModel(tree);
			
			// create Frequency Model
			Parameter freqs = new Parameter.Default(Utils.UNIFORM_CODON_FREQUENCIES);
			FrequencyModel freqModel = new FrequencyModel(Codons.UNIVERSAL,
					freqs);

			// create base subst model
			Parameter omegaParameter = new Parameter.Default("omega", 1, 1.0);
			Parameter kappaParameter = new Parameter.Default("kappa", 1, 1.0);
			GY94CodonModel baseSubModel = new GY94CodonModel(Codons.UNIVERSAL, omegaParameter, kappaParameter, freqModel);
			
			RandomBranchModel substitutionModel = new RandomBranchModel(treeModel, baseSubModel, 0.25, false, -1);
			
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

//			Sequence ancestralSequence = new Sequence();
//			ancestralSequence.appendSequenceString("TCAAGTGAGG");
//			partition1.setRootSequence(ancestralSequence);

			partitionsList.add(partition1);

			// feed to sequence simulator and generate data
			BeagleSequenceSimulator simulator = new BeagleSequenceSimulator(
					partitionsList);

			SimpleAlignment alignment = simulator.simulate(simulateInPar, true);

			// alignment.setOutputType(SimpleAlignment.OutputType.NEXUS);
			alignment.setOutputType(SimpleAlignment.OutputType.XML);

			System.out.println(alignment.toString());
	        
	        
		} catch (Exception e) {

			e.printStackTrace();

		}// END: try-catch
		
	}//END: simulateRandomBranchAssignment
	
	static void ancestralSequenceTree() {

		try {

			LinkedHashMap<NodeRef, int[]> sequenceMap = new LinkedHashMap<NodeRef, int[]>();
			DataType dataType = Nucleotides.INSTANCE;

			// create tree
			NewickImporter importer = new NewickImporter(
					"(SimSeq1:73.7468,(SimSeq2:25.256989999999995,SimSeq3:45.256989999999995):18.48981);");
			Tree tree = importer.importTree(null);
			TreeModel treeModel = new TreeModel(tree);

			for (NodeRef node : treeModel.getNodes()) {

				if (treeModel.isExternal(node)) {

					int[] seq = new int[] { 1, 1, 1 };
					sequenceMap.put(node, seq);

				} else {

					int[] seq = new int[] { 2, 2, 2 };
					sequenceMap.put(node, seq);

				}

			}// END: nodes loop

			AncestralSequenceTrait ancestralSequence = new AncestralSequenceTrait(sequenceMap, dataType);
			TreeTraitProvider[] treeTraitProviders = new TreeTraitProvider[] { ancestralSequence };

			StringBuffer buffer = new StringBuffer();
			NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);
			boolean useTipLabels = true;

			TreeUtils.newick(treeModel, //
					treeModel.getRoot(), //
					useTipLabels, //
					TreeUtils.BranchLengthType.LENGTHS_AS_TIME, //
					format, //
					null, //
					treeTraitProviders, //
					null, buffer);

			System.out.println(buffer);

		} catch (IOException e) {

			e.printStackTrace();

		} catch (ImportException e) {

			e.printStackTrace();

		}// END: try-catch

	}// END: annotateTree

	static void simulateTopology() {

		try {

			System.out.println("Test case 1: simulateTopology");

			MathUtils.setSeed(666);

			int sequenceLength = 10;
			ArrayList<Partition> partitionsList = new ArrayList<Partition>();

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
					partitionsList);

			System.out.println(simulator.simulate(simulateInPar, false).toString());

		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		} // END: try-catch block

	}// END: simulate topology

	static void simulateOnePartition() {

		try {

			MathUtils.setSeed(666);

			System.out.println("Test case 2: simulateOnePartition");

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
			partition1.setRootSequence(ancestralSequence);

			partitionsList.add(partition1);

			// feed to sequence simulator and generate data
			BeagleSequenceSimulator simulator = new BeagleSequenceSimulator(
					partitionsList);

			SimpleAlignment alignment = simulator.simulate(simulateInPar, false);

			// alignment.setOutputType(SimpleAlignment.OutputType.NEXUS);
			alignment.setOutputType(SimpleAlignment.OutputType.XML);

			System.out.println(alignment.toString());

		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		} // END: try-catch block

	}// END: simulateOnePartition

	static void simulateTwoPartitions() {

		try {

			System.out.println("Test case 3: simulateTwoPartitions");

			MathUtils.setSeed(666);

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
			Partition Partition = new Partition(treeModel, //
					substitutionModel,//
					siteRateModel, //
					branchRateModel, //
					freqModel, //
					4, // from
					sequenceLength - 1, // to
					1 // every
			);

			Sequence ancestralSequence = new Sequence();
			ancestralSequence.appendSequenceString("TCAAGTG");
			Partition.setRootSequence(ancestralSequence);

			partitionsList.add(partition1);
			partitionsList.add(Partition);

			// feed to sequence simulator and generate data
			BeagleSequenceSimulator simulator = new BeagleSequenceSimulator(
					partitionsList);

			System.out.println(simulator.simulate(simulateInPar, false).toString());

		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		} // END: try-catch block

	}// END: simulateTwoPartitions

	static void simulateThreePartitions(int i, int N) {

		try {

			MathUtils.setSeed(666);

			System.out.println("Test case 3: simulateThreePartitions");

			int sequenceLength = 100000;
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
					sequenceLength - 1, // to
					3 // every
			);

			// create partition
			Partition Partition = new Partition(treeModel, //
					substitutionModel,//
					siteRateModel, //
					branchRateModel, //
					freqModel, //
					1, // from
					sequenceLength - 1, // to
					3 // every
			);

			// create partition
			Partition partition3 = new Partition(treeModel, //
					substitutionModel,//
					siteRateModel, //
					branchRateModel, //
					freqModel, //
					2, // from
					sequenceLength - 1, // to
					3 // every
			);

			partitionsList.add(partition1);
			partitionsList.add(Partition);
			partitionsList.add(partition3);

			// feed to sequence simulator and generate data
			BeagleSequenceSimulator simulator = new BeagleSequenceSimulator(
					partitionsList);

			if (i == (N - 1)) {
				System.out
						.println(simulator.simulate(simulateInPar, false).toString());
			} else {
				simulator.simulate(simulateInPar, false);
			}

		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		} // END: try-catch block

	}// END: simulateThreePartitions

	static void simulateAminoAcid() {

		try {

			System.out.println("Test case 4: simulateAminoAcid");

			MathUtils.setSeed(666);

			int sequenceLength = 10;
			ArrayList<Partition> partitionsList = new ArrayList<Partition>();

			// create tree
			NewickImporter importer = new NewickImporter(
					"(SimSeq1:73.7468,(SimSeq2:25.256989999999995,SimSeq3:45.256989999999995):18.48981);");
			Tree tree = importer.importTree(null);
			TreeModel treeModel = new TreeModel(tree);

			// create site model
			GammaSiteRateModel siteRateModel = new GammaSiteRateModel(
					"siteModel");

			// create branch rate model
			BranchRateModel branchRateModel = new DefaultBranchRateModel();

			// create Frequency Model
			Parameter freqs = new Parameter.Default(new double[] { 0.05, 0.05,
					0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
					0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05 });
			FrequencyModel freqModel = new FrequencyModel(AminoAcids.INSTANCE,
					freqs);

			// create substitution model
			EmpiricalRateMatrix rateMatrix = Blosum62.INSTANCE;

			EmpiricalAminoAcidModel empiricalAminoAcidModel = new EmpiricalAminoAcidModel(
					rateMatrix, freqModel);

			HomogeneousBranchModel substitutionModel = new HomogeneousBranchModel(
					empiricalAminoAcidModel);

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
					partitionsList);

			System.out.println(simulator.simulate(simulateInPar, false).toString());

		} catch (Exception e) {

			e.printStackTrace();
			System.exit(-1);

		} // END: try-catch

	}// END: simulateAminoAcid

	static void simulateCodon() {

		try {

			boolean calculateLikelihood = true;

			System.out.println("Test case 6: simulate codons");

			MathUtils.setSeed(666);

			int sequenceLength = 10;
			ArrayList<Partition> partitionsList = new ArrayList<Partition>();

			// create tree
			NewickImporter importer = new NewickImporter(
					"(SimSeq1:73.7468,(SimSeq2:25.256989999999995,SimSeq3:45.256989999999995):18.48981);");
			Tree tree = importer.importTree(null);
			TreeModel treeModel = new TreeModel(tree);

			// create site model
			GammaSiteRateModel siteRateModel = new GammaSiteRateModel(
					"siteModel");

			// create branch rate model
			BranchRateModel branchRateModel = new DefaultBranchRateModel();

			// create Frequency Model
			Parameter freqs = new Parameter.Default(Utils.UNIFORM_CODON_FREQUENCIES);
			FrequencyModel freqModel = new FrequencyModel(Codons.UNIVERSAL,
					freqs);

			// create substitution model
			Parameter alpha = new Parameter.Default(1, 10);
			Parameter beta = new Parameter.Default(1, 5);
//			Parameter kappa = new Parameter.Default(1, 1);
			MG94CodonModel mg94 = new MG94CodonModel(Codons.UNIVERSAL, alpha,
					beta, freqModel);

			HomogeneousBranchModel substitutionModel = new HomogeneousBranchModel(
					mg94);

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
					partitionsList);

			Alignment alignment = simulator.simulate(simulateInPar, false);

			System.out.println(alignment.toString());

			if (calculateLikelihood) {

				// NewBeagleSequenceLikelihood nbtl = new
				// NewBeagleSequenceLikelihood(alignment, treeModel,
				// substitutionModel, (SiteModel) siteRateModel,
				// branchRateModel, null, false,
				// PartialsRescalingScheme.DEFAULT);

				ConvertAlignment convert = new ConvertAlignment(
						Nucleotides.INSTANCE, GeneticCode.UNIVERSAL, alignment);

				BeagleTreeLikelihood nbtl = new BeagleTreeLikelihood(convert, //
						treeModel, //
						substitutionModel, //
						siteRateModel, //
						branchRateModel, //
						null, //
						false, //
						PartialsRescalingScheme.DEFAULT, true);

				System.out.println("likelihood = " + nbtl.getLogLikelihood());

			}

		} catch (Exception e) {

			e.printStackTrace();
			System.exit(-1);

		} // END: try-catch

	}// END: simulateCodon

}// END: class
