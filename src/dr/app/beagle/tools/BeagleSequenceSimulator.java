/*
 * BeagleSequenceSimulator.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.beagle.tools;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.math.random.MersenneTwister;

import beagle.Beagle;
import beagle.BeagleFactory;
import dr.app.beagle.evomodel.branchmodel.BranchModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.treelikelihood.BufferIndexHelper;
import dr.app.beagle.evomodel.treelikelihood.SubstitutionModelDelegate;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
import dr.evolution.util.Taxon;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.math.MathUtils;

/**
 * @author Filip Bielejec
 * @version $Id$
 * 
 */
public class BeagleSequenceSimulator {

	List<MersenneTwister> randomList = new ArrayList<MersenneTwister>();

	private final boolean DEBUG = true;

	private ArrayList<Partition> partitions;
	private int siteCount;
	private SimpleAlignment simpleAlignment;
	private DataType dataType;
	private int stateCount;
	private ConcurrentHashMap<Taxon, int[]> alignmentMap;
	private boolean fieldsSet = false;
	private SubstitutionModelDelegate substitutionModelDelegate;

	private int gapFlag = Integer.MAX_VALUE;

	public BeagleSequenceSimulator(ArrayList<Partition> partitions) {

		this.partitions = partitions;
		this.alignmentMap = new ConcurrentHashMap<Taxon, int[]>();

		int siteCount = 0;
		int to = 0;
		for (Partition partition : partitions) {

			to = partition.to;
			if (to > siteCount) {
				siteCount = to;
			}
		}

		this.siteCount = siteCount + 1;
	}// END: Constructor

	private MersenneTwister nextMersenneTwister(int paritionNumber) {
		return new MersenneTwister(MathUtils.nextLong());
	}

	public Alignment simulate(boolean parallel) {

		try {

			// Executor for threads
			int NTHREDS = 1;
			if (parallel) {
				NTHREDS = Math.min(partitions.size(), Runtime.getRuntime()
						.availableProcessors());
			}

			for (int i = 0; i < partitions.size(); ++i) {
				randomList.add(nextMersenneTwister(i));
			}

			ExecutorService executor = Executors.newFixedThreadPool(NTHREDS);
			ThreadLocal<Partition> threadLocalPartition;

			ArrayList<Callable<Void>> simulatePartitionCallers = new ArrayList<Callable<Void>>();

			int partitionCount = 0;
			for (Partition partition : partitions) {

				threadLocalPartition = new ThreadLocal<Partition>();
				threadLocalPartition.set(partition);

				partition.setPartitionNumber(partitionCount);

				simulatePartitionCallers.add(new simulatePartitionCallable(
						threadLocalPartition.get()));
				threadLocalPartition.remove();
				partitionCount++;

			}// END: partitions loop

			executor.invokeAll(simulatePartitionCallers);

			// Wait until all threads are finished
			executor.shutdown();
			while (!executor.isTerminated()) {
			}

			// TODO: simplify
			// compile the alignment
			for (Partition partition : partitions) {

				Map<Taxon, int[]> sequenceMap = partition.getSequenceList();
				Iterator<Entry<Taxon, int[]>> iterator = sequenceMap.entrySet()
						.iterator();

				while (iterator.hasNext()) {

					Entry<Taxon, int[]> pairs = (Entry<Taxon, int[]>) iterator
							.next();

					Taxon taxon = pairs.getKey();
					int[] partitionSequence = pairs.getValue();

					if (alignmentMap.containsKey(taxon)) {

						int j = 0;
						for (int i = partition.from; i <= partition.to; i += partition.every) {

							alignmentMap.get(taxon)[i] = partitionSequence[j];
							j++;

						}// END: i loop

					} else {

						int[] sequence = new int[siteCount];
						// dirty solution for gaps when taxa between the tree topologies don't match
						Arrays.fill(sequence, gapFlag);

						int j = 0;
						for (int i = partition.from; i <= partition.to; i += partition.every) {

							sequence[i] = partitionSequence[j];
							j++;

						}// END: i loop

						alignmentMap.put(taxon, sequence);
					}// END: key check

				}// END: iterator loop

			}// END: partitions loop

			Iterator<Entry<Taxon, int[]>> iterator = alignmentMap.entrySet().iterator();
			while (iterator.hasNext()) {

				Entry<Taxon, int[]> pairs = (Entry<Taxon, int[]>) iterator
						.next();
				simpleAlignment.addSequence(intArray2Sequence(
						(Taxon) pairs.getKey(), (int[]) pairs.getValue(),
						gapFlag));

				iterator.remove();

			}// END: while has next

		} catch (Exception e) {
			e.printStackTrace();
		}// END: try-catch block

		return simpleAlignment;
	}// END: simulate

//	public Alignment simulate() {
//
//		try {
//
//			int partitionCount = 0;
//			for (Partition partition : partitions) {
//
//				// String partitionName = "partition " + partitionCount;
//				partition.setPartitionNumber(partitionCount);// (partitionName);
//
//				simulatePartition(partition);
//
//				partitionCount++;
//			}// END: partitions loop
//
//			if (DEBUG) {
//				Utils.printHashMap(alignmentMap);
//			}
//
//			// compile the alignment
//			Iterator<Entry<Taxon, int[]>> iterator = alignmentMap.entrySet()
//					.iterator();
//			while (iterator.hasNext()) {
//
//				Entry<?, ?> pairs = (Entry<?, ?>) iterator.next();
//				Taxon taxon = (Taxon) pairs.getKey();
//				int[] sequence = (int[]) pairs.getValue();
//
//				simpleAlignment.addSequence(intArray2Sequence(taxon, sequence,
//						gapFlag));
//				iterator.remove();
//
//			}// END: while has next
//
//		} catch (Exception e) {
//			e.printStackTrace();
//		}// END: try-catch block
//
//		return simpleAlignment;
//	}// END: simulate

	private class simulatePartitionCallable implements Callable<Void> {

		private Partition partition;

		private simulatePartitionCallable(Partition partition) {
			this.partition = partition;
		}// END: Constructor

		public Void call() {

			try {

				simulatePartition(partition);

			} catch (Exception e) {
				e.printStackTrace();
			}

			return null;
		}// END: call

	}// END: simulatePartitionCallable class

	private void simulatePartition(Partition partition) {

		TreeModel treeModel = partition.getTreeModel();
		BranchModel branchModel = partition.getBranchModel();
		GammaSiteRateModel siteModel = partition.getSiteModel();
		FrequencyModel freqModel = partition.getFreqModel();
		int partitionSiteCount = partition.getPartitionSiteCount();

		NodeRef root = treeModel.getRoot();

		// do those only once
		if (!fieldsSet) {

			dataType = freqModel.getDataType();
			simpleAlignment = new SimpleAlignment();
			simpleAlignment.setDataType(dataType);
			simpleAlignment.setReportCountStatistics(false);
			stateCount = dataType.getStateCount();

			fieldsSet = true;
		}// END: partitionCount check

		// Buffer index helpers
		int eigenCount = branchModel.getModelCount();
		BufferIndexHelper eigenBufferHelper = new BufferIndexHelper(eigenCount,
				0);

		int nodeCount = treeModel.getNodeCount();
		BufferIndexHelper matrixBufferHelper = new BufferIndexHelper(nodeCount,
				0);

		int tipCount = treeModel.getExternalNodeCount();
		BufferIndexHelper partialBufferHelper = new BufferIndexHelper(
				nodeCount, tipCount);

		substitutionModelDelegate = new SubstitutionModelDelegate(treeModel,
				branchModel);

		// load beagle
		Beagle beagle;
		synchronized (this) {
			beagle = loadBeagleInstance(partition, eigenBufferHelper,
					matrixBufferHelper, partialBufferHelper);
		}

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

				category[i] = randomChoicePDF(categoryProbs, partition.getPartitionNumber());
				// category[i] = MathUtils.randomChoicePDF(categoryProbs);

		}

		int[] parentSequence = new int[partitionSiteCount];

		// set ancestral sequence for partition if it exists
		if (partition.hasAncestralSequence) {

			parentSequence = sequence2intArray(partition.getAncestralSequence());

		} else {

			double[] frequencies = freqModel.getFrequencies();
			for (int i = 0; i < partitionSiteCount; i++) {

					parentSequence[i] = randomChoicePDF(frequencies, partition.getPartitionNumber());
					// parentSequence[i] = MathUtils.randomChoicePDF(frequencies);

			}

		}// END: ancestral sequence check

		substitutionModelDelegate.updateSubstitutionModels(beagle);

		int categoryCount = siteModel.getCategoryCount();
		traverse(beagle, partition, root, parentSequence, category,
				categoryCount, matrixBufferHelper, eigenBufferHelper);

		if (DEBUG) {
			synchronized (this) {
				partition.printSequences();
			}

		}

	}// END: simulatePartition

	public Beagle loadBeagleInstance(Partition partition, //
			BufferIndexHelper eigenBufferHelper, //
			BufferIndexHelper matrixBufferHelper, //
			BufferIndexHelper partialBufferHelper //
	) {

		TreeModel treeModel = partition.getTreeModel();
		GammaSiteRateModel siteModel = partition.getSiteModel();
		int partitionSiteCount = partition.getPartitionSiteCount();

		int tipCount = treeModel.getExternalNodeCount();
		int compactPartialsCount = tipCount;
		int patternCount = partitionSiteCount;
		int categoryCount = siteModel.getCategoryCount();
		int internalNodeCount = treeModel.getInternalNodeCount();

		int[] resourceList = new int[] { 0 };
		long preferenceFlags = 0;
		long requirementFlags = 0;

		BufferIndexHelper scaleBufferHelper = new BufferIndexHelper(
				internalNodeCount + 1, 0);

		Beagle beagle = BeagleFactory.loadBeagleInstance(tipCount,
				partialBufferHelper.getBufferCount(), compactPartialsCount,
				stateCount, patternCount,
				substitutionModelDelegate.getEigenBufferCount(),
				substitutionModelDelegate.getMatrixBufferCount(),
				categoryCount, scaleBufferHelper.getBufferCount(),
				resourceList, preferenceFlags, requirementFlags);

		return beagle;
	}// END: loadBeagleInstance

	private void traverse(Beagle beagle, //
			Partition partition, //
			NodeRef node, //
			int[] parentSequence, //
			int[] category, //
			int categoryCount, //
			BufferIndexHelper matrixBufferHelper, //
			BufferIndexHelper eigenBufferHelper //
	) {

		TreeModel treeModel = partition.getTreeModel();
		int partitionSiteCount = partition.getPartitionSiteCount();

		for (int iChild = 0; iChild < treeModel.getChildCount(node); iChild++) {

			NodeRef child = treeModel.getChild(node, iChild);
			int[] partitionSequence = new int[partitionSiteCount];
			double[] cProb = new double[stateCount];

			double[][] probabilities = getTransitionProbabilities(beagle,
					partition, child, categoryCount, stateCount,
					matrixBufferHelper, eigenBufferHelper);

			for (int i = 0; i < partitionSiteCount; i++) {

				System.arraycopy(probabilities[category[i]], parentSequence[i]
						* stateCount, cProb, 0, stateCount);

				// if(BeagleSeqSimTest.unofficialMT) {
				partitionSequence[i] = randomChoicePDF(cProb,
						partition.getPartitionNumber());
				// } else {
				// partitionSequence[i] = MathUtils.randomChoicePDF(cProb);
				// }

			}

			if (treeModel.getChildCount(child) == 0) {

				Taxon taxon = treeModel.getNodeTaxon(child);
				partition.getSequenceList().put(taxon, partitionSequence);

			} // END: tip node check

			traverse(beagle, partition, treeModel.getChild(node, iChild),
					partitionSequence, category, categoryCount,
					matrixBufferHelper, eigenBufferHelper);

		}// END: child nodes loop

	}// END: traverse

	private double[][] getTransitionProbabilities(Beagle beagle, //
			Partition partition, //
			NodeRef node, //
			int categoryCount, //
			int stateCount, //
			BufferIndexHelper matrixBufferHelper, //
			BufferIndexHelper eigenBufferHelper //
	) {

		double[][] probabilities = new double[categoryCount][stateCount
				* stateCount];

		TreeModel treeModel = partition.getTreeModel();
		BranchRateModel branchRateModel = partition.getBranchRateModel();

		int nodeNum = node.getNumber();
		matrixBufferHelper.flipOffset(nodeNum);
		int branchIndex = nodeNum;

		double branchRate = branchRateModel.getBranchRate(treeModel, node);
		double branchTime = treeModel.getBranchLength(node) * branchRate;

		if (branchTime < 0.0) {
			throw new RuntimeException("Negative branch length: " + branchTime);
		}

		int count = 1;
		substitutionModelDelegate.updateTransitionMatrices(beagle,
				new int[] { branchIndex }, new double[] { branchTime }, count);

		double transitionMatrix[] = new double[categoryCount * stateCount
				* stateCount];

		beagle.getTransitionMatrix(branchIndex, //
				transitionMatrix //
		);

		for (int i = 0; i < categoryCount; i++) {

			System.arraycopy(transitionMatrix, i * stateCount,
					probabilities[i], 0, stateCount * stateCount);

		}

		return probabilities;
	}// END: getTransitionProbabilities

	private Sequence intArray2Sequence(Taxon taxon, int[] seq, int gapFlag) {

		StringBuilder sSeq = new StringBuilder();

		if (dataType instanceof Codons) {

			for (int i = 0; i < siteCount; i++) {

				int state = seq[i];

				if (state == gapFlag) {
					sSeq.append(dataType.getTriplet(dataType.getGapState()));
				} else {
					sSeq.append(dataType.getTriplet(seq[i]));
				}// END: gap check

			}// END: replications loop

		} else {

			for (int i = 0; i < siteCount; i++) {

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

		int array[] = new int[siteCount];

		if (dataType instanceof Codons) {

			int k = 0;
			for (int i = 0; i < siteCount; i++) {
				array[i] = ((Codons) dataType).getState(sequence.getChar(k),
						sequence.getChar(k + 1), sequence.getChar(k + 2));
				k += 3;
			}// END: replications loop

		} else {

			for (int i = 0; i < siteCount; i++) {
				array[i] = dataType.getState(sequence.getChar(i));
			}// END: replications loop

		}// END: dataType check

		return array;
	}// END: sequence2intArray

	// //////////////////////////////
	// ---BEGIN: DUPLICATED CODE---//
	// //////////////////////////////

	private double getTotal(double[] array, int start, int end) {
		double total = 0.0;
		for (int i = start; i < end; i++) {
			total += array[i];
		}
		return total;
	}

	private int randomChoicePDF(double[] pdf, int partitionNumber) {

		MersenneTwister random = randomList.get(partitionNumber);

		double U = random.nextDouble() * getTotal(pdf, 0, pdf.length);
		for (int i = 0; i < pdf.length; i++) {

			U -= pdf[i];
			if (U < 0.0) {
				return i;
			}

		}
		for (int i = 0; i < pdf.length; i++) {
			System.out.println(i + "\t" + pdf[i]);
		}
		
		System.exit(-1);
		throw new RuntimeException(
				"randomChoiceUnnormalized falls through -- negative components in input distribution?");
	}

	// ////////////////////////////
	// ---END: DUPLICATED CODE---//
	// ////////////////////////////

} // END: class
