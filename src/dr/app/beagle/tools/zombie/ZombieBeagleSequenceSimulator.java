/*
 * ZombieBeagleSequenceSimulator.java
 *
 * Copyright (C) 2002-2012 Alexei Drummond, Andrew Rambaut & Marc A. Suchard
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

package dr.app.beagle.tools.zombie;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;

import beagle.Beagle;
import beagle.BeagleFactory;
import dr.app.beagle.evomodel.sitemodel.BranchSubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.EpochBranchSubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.sitemodel.HomogenousBranchSubstitutionModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.HKY;
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
import dr.evolution.util.Taxon;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

/**
 * @author Filip Bielejec
 * @version $Id$
 * 
 */
public class ZombieBeagleSequenceSimulator {

	private final boolean DEBUG = false;
	
	private ArrayList<ZombiePartition> partitions;
	private int replications;
    private SimpleAlignment simpleAlignment;
    private DataType dataType;
    private int stateCount;
	private ConcurrentHashMap<Taxon, int[]> alignmentMap;
    private boolean fieldsSet = false;
	
    private int gapFlag = Integer.MAX_VALUE;
    
	public ZombieBeagleSequenceSimulator(ArrayList<ZombiePartition> partitions, //
			int replications //
	) {

		this.partitions = partitions;
		this.replications = replications;

		alignmentMap = new ConcurrentHashMap<Taxon, int[]>();
		
	}// END: Constructor

	public Alignment simulate() {

		try {

			int partitionCount = 0;
			for (ZombiePartition partition : partitions) {

//				threadLocalPartition = new ThreadLocal<Partition>();
//				threadLocalPartition.set(partition);

				if (DEBUG) {
					System.out.println("Simulating for partition " + partitionCount);
				}

				simulatePartition(partition);
//				executor.submit(new simulatePartitionRunnable(threadLocalPartition.get()));
//				simulatePartitionCallers.add(new simulatePartitionRunnable(threadLocalPartition.get()));

//				threadLocalPartition.remove();

				partitionCount++;
			}// END: partitions loop

			// Wait until all threads are finished
//			executor.shutdown();
//			while (!executor.isTerminated()) {
//			}

			if (DEBUG) {
				printHashMap(alignmentMap);
			}

			// compile the alignment
			Iterator<Entry<Taxon, int[]>> iterator = alignmentMap.entrySet().iterator();
			while (iterator.hasNext()) {

				Entry<?, ?> pairs = (Entry<?, ?>) iterator.next();
//				Taxon taxon = (Taxon) pairs.getKey();
//				int[] sequence = (int[]) pairs.getValue();
				
				simpleAlignment.addSequence(intArray2Sequence((Taxon) pairs.getKey(), (int[]) pairs.getValue(), gapFlag));
				iterator.remove();

			}// END: while has next

		} catch (Exception e) {
			e.printStackTrace();
		}// END: try-catch block

		return simpleAlignment;
	}// END: simulate

	private void simulatePartition(ZombiePartition partition) {

		TreeModel treeModel = partition.treeModel;
		BranchSubstitutionModel branchSubstitutionModel = partition.branchSubstitutionModel;
		GammaSiteRateModel siteModel = partition.siteModel;
		FrequencyModel freqModel = partition.freqModel;
		int partitionSiteCount = partition.getPartitionSiteCount();
		
		NodeRef root = treeModel.getRoot();

		// do those only once
		if(!fieldsSet) {
			
			dataType = freqModel.getDataType();
			simpleAlignment = new SimpleAlignment();
			simpleAlignment.setDataType(dataType);
			simpleAlignment.setReportCountStatistics(false);
			stateCount = dataType.getStateCount();
			
			fieldsSet = true;
		}//END: partitionCount check
		
		// Buffer index helpers
		int eigenCount = branchSubstitutionModel.getEigenCount();
		BufferIndexHelper eigenBufferHelper = new BufferIndexHelper(eigenCount, 0);
		
		int nodeCount = treeModel.getNodeCount();
		BufferIndexHelper matrixBufferHelper = new BufferIndexHelper(nodeCount, 0);
		
		int tipCount = treeModel.getExternalNodeCount();
		BufferIndexHelper partialBufferHelper = new BufferIndexHelper(nodeCount, tipCount);

		// load beagle
		Beagle beagle = loadBeagleInstance(partition, eigenBufferHelper,
				matrixBufferHelper, partialBufferHelper);

		// gamma category rates
		double[] categoryRates = siteModel.getCategoryRates();
		beagle.setCategoryRates(categoryRates);

		// weights for gamma category rates
		double[] categoryWeights = siteModel.getCategoryProportions();
		beagle.setCategoryWeights(0, categoryWeights);

		// proportion of sites in each category
		double[] categoryProbs = siteModel.getCategoryProportions();
		int[] category = new int[partitionSiteCount];

		for (int i = 0; i < partitionSiteCount; i++) {
			category[i] = MathUtils.randomChoicePDF(categoryProbs);
		}

		int[] parentSequence = new int[partitionSiteCount];

		// set ancestral sequence for partition if it exists
		if (partition.hasAncestralSequence) {

			parentSequence = sequence2intArray(partition.ancestralSequence);

		} else {

			double[] frequencies = freqModel.getFrequencies();
			for (int i = 0; i < partitionSiteCount; i++) {
				parentSequence[i] = MathUtils.randomChoicePDF(frequencies);
			}

		}// END: ancestral sequence check
		
		for (int i = 0; i < eigenCount; i++) {

			eigenBufferHelper.flipOffset(i);

			branchSubstitutionModel.setEigenDecomposition(beagle, //
					i, //
					eigenBufferHelper, //
					0 //
					);

		}// END: i loop

		int categoryCount = siteModel.getCategoryCount();
		traverse(beagle, partition, root, parentSequence, category, categoryCount, matrixBufferHelper, eigenBufferHelper);

	}// END: simulatePartition

	public Beagle loadBeagleInstance(ZombiePartition partition, //
			BufferIndexHelper eigenBufferHelper, //
			BufferIndexHelper matrixBufferHelper, //
			BufferIndexHelper partialBufferHelper //
			) {

		TreeModel treeModel = partition.treeModel;
		BranchSubstitutionModel branchSubstitutionModel = partition.branchSubstitutionModel;
		GammaSiteRateModel siteModel = partition.siteModel;
		int partitionSiteCount = partition.getPartitionSiteCount();
		
		int tipCount = treeModel.getExternalNodeCount();
		int compactPartialsCount = tipCount;
		int patternCount = partitionSiteCount;
		int categoryCount = siteModel.getCategoryCount();
		int internalNodeCount = treeModel.getInternalNodeCount();
		int scaleBufferCount = internalNodeCount + 1;

		int[] resourceList = new int[] { 0 };
		long preferenceFlags = 0;
		long requirementFlags = 0;

		Beagle beagle = BeagleFactory.loadBeagleInstance(
				tipCount, //
				partialBufferHelper.getBufferCount(), //
				compactPartialsCount, //
				stateCount, //
				patternCount, //
				eigenBufferHelper.getBufferCount(), //
				matrixBufferHelper.getBufferCount() + branchSubstitutionModel.getExtraBufferCount(treeModel), //
				categoryCount, //
				scaleBufferCount, //
				resourceList, //
				preferenceFlags, //
				requirementFlags //
				);

		return beagle;
	}// END: loadBeagleInstance
	
	// TODO: traverse
	private void traverse(Beagle beagle, //
			ZombiePartition partition, //
			NodeRef node, //
			int[] parentSequence, //
			int[] category, //
			int categoryCount, //
			BufferIndexHelper matrixBufferHelper, //
			BufferIndexHelper eigenBufferHelper //
	) {

		TreeModel treeModel = partition.treeModel;
		int partitionSiteCount = partition.getPartitionSiteCount();
		
		for (int iChild = 0; iChild < treeModel.getChildCount(node); iChild++) {

			NodeRef child = treeModel.getChild(node, iChild);
			int[] partitionSequence = new int[partitionSiteCount];
			double[] cProb = new double[stateCount];

			double[][] probabilities = getTransitionProbabilities(beagle, partition, child, categoryCount, stateCount, matrixBufferHelper, eigenBufferHelper);

			if(DEBUG) {
				System.out.println("parent sequence:");
				EpochBranchSubstitutionModel.printArray(parentSequence);
			}
			
			for (int i = 0; i < partitionSiteCount; i++) {

				System.arraycopy(probabilities[category[i]], parentSequence[i] * stateCount, cProb, 0, stateCount);
				partitionSequence[i] = MathUtils.randomChoicePDF(cProb);
				
			}
			
			if (DEBUG) {
				System.out.println("partition sequence:");
				EpochBranchSubstitutionModel.printArray(partitionSequence);
			}
			
			if (treeModel.getChildCount(child) == 0) {

				Taxon taxon = treeModel.getNodeTaxon(child);

				if (alignmentMap.containsKey(taxon)) {

					int j = 0;
					for (int i = partition.from; i <= partition.to; i += partition.every) {

						alignmentMap.get(taxon)[i] = partitionSequence[j];
						j++;

					}// END: i loop

				} else {

					int[] sequence = new int[replications];
					// TODO
					Arrays.fill(sequence, gapFlag);

					int j = 0;
					for (int i = partition.from; i <= partition.to; i += partition.every) {

						sequence[i] = partitionSequence[j];
						j++;

					}// END: i loop

					alignmentMap.put(taxon, sequence);
				}// END: key check			
				
			} //END: tip node check

			traverse(beagle, partition, treeModel.getChild(node, iChild), partitionSequence, category, categoryCount, matrixBufferHelper, eigenBufferHelper);

		}// END: child nodes loop

	}// END: traverse

	private double[][] getTransitionProbabilities(Beagle beagle, //
			ZombiePartition partition, //
			NodeRef node, //
			int categoryCount, //
			int stateCount, //
			BufferIndexHelper matrixBufferHelper, //
			BufferIndexHelper eigenBufferHelper //
	) {

		double[][] probabilities = new double[categoryCount][stateCount * stateCount];

		BranchSubstitutionModel branchSubstitutionModel = partition.branchSubstitutionModel;
		TreeModel treeModel = partition.treeModel;
		BranchRateModel branchRateModel = partition.branchRateModel;
		
		int nodeNum = node.getNumber();
		matrixBufferHelper.flipOffset(nodeNum);
		int branchIndex = matrixBufferHelper.getOffsetIndex(nodeNum);
		int eigenIndex = branchSubstitutionModel.getBranchIndex(treeModel, node, branchIndex);
		int count = 1;

		double branchRate = branchRateModel.getBranchRate(treeModel, node);
		double branchTime = treeModel.getBranchLength(node) * branchRate;
		
        if (branchTime < 0.0) {
            throw new RuntimeException("Negative branch length: " + branchTime);
        }
		
		branchSubstitutionModel.updateTransitionMatrices(beagle, //
				eigenIndex, //
				eigenBufferHelper, //
				new int[] { branchIndex }, //
				null, //
				null, //
				new double[] { branchTime }, //
				count //
				);

		double transitionMatrix[] = new double[categoryCount * stateCount * stateCount];
		
		beagle.getTransitionMatrix(branchIndex, //
				transitionMatrix //
				);

		for (int i = 0; i < categoryCount; i++) {

			System.arraycopy(transitionMatrix, i * stateCount, probabilities[i], 0, stateCount * stateCount);
		
		}

		return probabilities;
	}// END: getTransitionProbabilities

	private Sequence intArray2Sequence(Taxon taxon, int[] seq, int gapFlag) {

		StringBuilder sSeq = new StringBuilder();

		if (dataType instanceof Codons) {

			for (int i = 0; i < replications; i++) {

				int state = seq[i];

				if (state == gapFlag) {
					sSeq.append(dataType.getTriplet(dataType.getGapState()));
				} else {
					sSeq.append(dataType.getTriplet(seq[i]));
				}// END: gap check

			}// END: replications loop

		} else {

			for (int i = 0; i < replications; i++) {

				int state = seq[i];

				if (state == gapFlag) {
					sSeq.append(dataType.getCode(dataType.getGapState()));
				} else {
					sSeq.append(dataType.getCode(seq[i]));
				}// END: gap check

			}// END: replications loop

		}// END: dataType check

		return new Sequence(taxon, sSeq.toString());
	}// END: intArray2Sequence
	
	private int[] sequence2intArray(Sequence sequence) {

		int array[] = new int[replications];

		if (dataType instanceof Codons) {

			int k = 0;
			for (int i = 0; i < replications; i++) {
				array[i] = ((Codons) dataType).getState(sequence.getChar(k), sequence.getChar(k + 1), sequence.getChar(k + 2));
				k += 3;
			}// END: replications loop

		} else {

			for (int i = 0; i < replications; i++) {
				array[i] = dataType.getState(sequence.getChar(i));
			}// END: replications loop

		}// END: dataType check

		return array;
	}// END: sequence2intArray

	// /////////////////
	// ---DEBUGGING---//
	// /////////////////

	public void printHashMap(ConcurrentHashMap<?, ?> hashMap) {
		
		Iterator<?> iterator = hashMap.entrySet().iterator();
		while (iterator.hasNext()) {

			   Entry<?, ?> pairs = (Entry<?, ?>) iterator.next();
		       
			   Taxon taxon = (Taxon)pairs.getKey(); 
			   int[] sequence = (int[])pairs.getValue();
			   
			   System.out.println(taxon.toString());
			   EpochBranchSubstitutionModel.printArray(sequence);
			   
		}// END: while has next
		
	}//END: printHashMap
	
	public static void main(String[] args) {

		simulateOnePartition();
		simulateTwoPartitions();
		simulateThreePartitions(); 

	} // END: main

	static void simulateOnePartition() {

		try {

			System.out.println("Test case 1: simulateOnePartition");
			
			int sequenceLength = 10;
			ArrayList<ZombiePartition> partitionsList = new ArrayList<ZombiePartition>();

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
			HomogenousBranchSubstitutionModel substitutionModel = new HomogenousBranchSubstitutionModel(
					hky, freqModel);

			// create site model
			GammaSiteRateModel siteRateModel = new GammaSiteRateModel(
					"siteModel");

			// create branch rate model
			BranchRateModel branchRateModel = new DefaultBranchRateModel();

			// create partition
			ZombiePartition partition1 = new ZombiePartition(treeModel, //
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
			ZombieBeagleSequenceSimulator simulator = new ZombieBeagleSequenceSimulator(partitionsList, sequenceLength);
			System.out.println(simulator.simulate().toString());
			
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		} // END: try-catch block

	}// END: simulateOnePartition

	static void simulateTwoPartitions() {

		try {

			System.out.println("Test case 2: simulateTwoPartitions");
			
			int sequenceLength = 11;
			ArrayList<ZombiePartition> partitionsList = new ArrayList<ZombiePartition>();

			// create tree
			NewickImporter importer = new NewickImporter(
					"(SimSeq1:73.7468,(SimSeq2:25.256989999999995,SimSeq3:45.256989999999995):18.48981);");
			Tree tree = importer.importTree(null);
			TreeModel treeModel = new TreeModel(tree);

			// create Frequency Model
			Parameter freqs = new Parameter.Default(new double[] { 0.25, 0.25, 0.25, 0.25 });
			FrequencyModel freqModel = new FrequencyModel(Nucleotides.INSTANCE,
					freqs);

			// create substitution model
			Parameter kappa = new Parameter.Default(1, 10);
			HKY hky = new HKY(kappa, freqModel);
			HomogenousBranchSubstitutionModel substitutionModel = new HomogenousBranchSubstitutionModel(
					hky, freqModel);

			// create site model
			GammaSiteRateModel siteRateModel = new GammaSiteRateModel(
					"siteModel");

			// create branch rate model
			BranchRateModel branchRateModel = new DefaultBranchRateModel();

			// create partition
			ZombiePartition partition1 = new ZombiePartition(treeModel, //
					substitutionModel, //
					siteRateModel, //
					branchRateModel, //
					freqModel, //
					0, // from
					4, // to
					1 // every
			);

			// create partition
			ZombiePartition partition2 = new ZombiePartition(treeModel, //
					substitutionModel,//
					siteRateModel, //
					branchRateModel, //
					freqModel, //
					5, // from
					sequenceLength - 1, // to
					1 // every
			);
		
//			Sequence ancestralSequence = new Sequence();
//			ancestralSequence.appendSequenceString("TCAAGTGAGG");
//			partition2.setAncestralSequence(ancestralSequence);
			
			partitionsList.add(partition1);
			partitionsList.add(partition2);
			
			// feed to sequence simulator and generate data
			ZombieBeagleSequenceSimulator simulator = new ZombieBeagleSequenceSimulator(partitionsList, sequenceLength);
			System.out.println(simulator.simulate().toString());
			
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		} // END: try-catch block

	}// END: simulateTwoPartitions
	
	static void simulateThreePartitions() {

		try {

			System.out.println("Test case 3: simulateThreePartitions");
			
			int sequenceLength = 10;
			ArrayList<ZombiePartition> partitionsList = new ArrayList<ZombiePartition>();

			// create tree
			NewickImporter importer = new NewickImporter(
					"(SimSeq1:73.7468,(SimSeq2:25.256989999999995,SimSeq3:45.256989999999995):18.48981);");
			Tree tree = importer.importTree(null);
			TreeModel treeModel = new TreeModel(tree);

			// create Frequency Model
			Parameter freqs = new Parameter.Default(new double[] { 0.25, 0.25, 0.25, 0.25 });
			FrequencyModel freqModel = new FrequencyModel(Nucleotides.INSTANCE,
					freqs);

			// create substitution model
			Parameter kappa = new Parameter.Default(1, 10);
			HKY hky = new HKY(kappa, freqModel);
			HomogenousBranchSubstitutionModel substitutionModel = new HomogenousBranchSubstitutionModel(
					hky, freqModel);

			// create site model
			GammaSiteRateModel siteRateModel = new GammaSiteRateModel(
					"siteModel");

			// create branch rate model
			BranchRateModel branchRateModel = new DefaultBranchRateModel();

			// create partition
			ZombiePartition partition1 = new ZombiePartition(treeModel, //
					substitutionModel, //
					siteRateModel, //
					branchRateModel, //
					freqModel, //
					0, // from
					sequenceLength - 1, // to
					3 // every
			);

			// create partition
			ZombiePartition partition2 = new ZombiePartition(treeModel, //
					substitutionModel,//
					siteRateModel, //
					branchRateModel, //
					freqModel, //
					1, // from
					sequenceLength - 1, // to
					3 // every
			);
			
			// create partition
			ZombiePartition partition3 = new ZombiePartition(treeModel, //
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
			ZombieBeagleSequenceSimulator simulator = new ZombieBeagleSequenceSimulator(partitionsList, sequenceLength);
			System.out.println(simulator.simulate().toString());
			
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		} // END: try-catch block

	}// END: simulateThreePartitions
	
} // END: class
