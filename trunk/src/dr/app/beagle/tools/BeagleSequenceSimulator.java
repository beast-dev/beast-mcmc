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
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;

import beagle.Beagle;
import beagle.BeagleFactory;
import dr.app.beagle.evomodel.branchmodel.BranchModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.treelikelihood.BufferIndexHelper;
import dr.app.beagle.evomodel.treelikelihood.SubstitutionModelDelegate;
import dr.app.bss.Utils;
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

	private final boolean DEBUG = false;
	
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
			// siteCount += partition.getPartitionSiteCount();
			to = partition.to;
			if (to > siteCount) {
				siteCount = to;
			}
		}
		this.siteCount = siteCount + 1;

	}// END: Constructor

	public Alignment simulate() {

		try {

			int partitionCount = 0;
			for (Partition partition : partitions) {

				if (DEBUG) {
					System.out.println("Simulating for partition " + partitionCount);
				}

				simulatePartition(partition);

				partitionCount++;
			}// END: partitions loop

			if (DEBUG) {
				Utils.printHashMap(alignmentMap);
			}

			// compile the alignment
			Iterator<Entry<Taxon, int[]>> iterator = alignmentMap.entrySet().iterator();
			while (iterator.hasNext()) {

				Entry<?, ?> pairs = (Entry<?, ?>) iterator.next();
				Taxon taxon = (Taxon) pairs.getKey();
				int[] sequence = (int[]) pairs.getValue();
				
				simpleAlignment.addSequence(intArray2Sequence(taxon, sequence, gapFlag));
				iterator.remove();

			}// END: while has next

		} catch (Exception e) {
			e.printStackTrace();
		}// END: try-catch block

		return simpleAlignment;
	}// END: simulate

	private void simulatePartition(Partition partition) {

		TreeModel treeModel = partition.treeModel;
		BranchModel branchModel = partition.getBranchModel();
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
		int eigenCount = branchModel.getModelCount();
		BufferIndexHelper eigenBufferHelper = new BufferIndexHelper(eigenCount, 0);
		
		int nodeCount = treeModel.getNodeCount();
		BufferIndexHelper matrixBufferHelper = new BufferIndexHelper(nodeCount, 0);
		
		int tipCount = treeModel.getExternalNodeCount();
		BufferIndexHelper partialBufferHelper = new BufferIndexHelper(nodeCount, tipCount);

		substitutionModelDelegate = new SubstitutionModelDelegate(treeModel, branchModel);
		
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
		
		substitutionModelDelegate.updateSubstitutionModels(beagle);
//		for (int i = 0; i < eigenCount; i++) {
//
//			eigenBufferHelper.flipOffset(i);
//
//			branchModel.setEigenDecomposition(beagle, //
//					i, //
//					eigenBufferHelper, //
//					0 //
//					);
//
//		}// END: i loop

		int categoryCount = siteModel.getCategoryCount();
		traverse(beagle, partition, root, parentSequence, category, categoryCount, matrixBufferHelper, eigenBufferHelper);

	}// END: simulatePartition

	public Beagle loadBeagleInstance(Partition partition, //
			BufferIndexHelper eigenBufferHelper, //
			BufferIndexHelper matrixBufferHelper, //
			BufferIndexHelper partialBufferHelper //
			) {

		TreeModel treeModel = partition.treeModel;
		GammaSiteRateModel siteModel = partition.siteModel;
		int partitionSiteCount = partition.getPartitionSiteCount();
		
		int tipCount = treeModel.getExternalNodeCount();
		int compactPartialsCount = tipCount;
		int patternCount = partitionSiteCount;
		int categoryCount = siteModel.getCategoryCount();
		int internalNodeCount = treeModel.getInternalNodeCount();

		int[] resourceList = new int[] { 0 };
		long preferenceFlags = 0;
		long requirementFlags = 0;

        // one scaling buffer for each internal node plus an extra for the accumulation, then doubled for store/restore
		BufferIndexHelper  scaleBufferHelper = new BufferIndexHelper(internalNodeCount + 1, 0);
		
		Beagle beagle = BeagleFactory.loadBeagleInstance(
                tipCount,
                partialBufferHelper.getBufferCount(),
                compactPartialsCount,
                stateCount,
                patternCount,
                substitutionModelDelegate.getEigenBufferCount(),
                substitutionModelDelegate.getMatrixBufferCount(),
                categoryCount,
                scaleBufferHelper.getBufferCount(), // Always allocate; they may become necessary
                resourceList,
                preferenceFlags,
                requirementFlags
        );
		
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

		TreeModel treeModel = partition.treeModel;
		int partitionSiteCount = partition.getPartitionSiteCount();
		
		for (int iChild = 0; iChild < treeModel.getChildCount(node); iChild++) {

			NodeRef child = treeModel.getChild(node, iChild);
			int[] partitionSequence = new int[partitionSiteCount];
			double[] cProb = new double[stateCount];

			double[][] probabilities = getTransitionProbabilities(beagle, partition, child, categoryCount, stateCount, matrixBufferHelper, eigenBufferHelper);

			if(DEBUG) {
				System.out.println("parent sequence:");
				Utils.printArray(parentSequence);
			}
			
			for (int i = 0; i < partitionSiteCount; i++) {

				System.arraycopy(probabilities[category[i]], parentSequence[i] * stateCount, cProb, 0, stateCount);
				partitionSequence[i] = MathUtils.randomChoicePDF(cProb);
				
			}
			
			if (DEBUG) {
				System.out.println("partition sequence:");
				Utils.printArray(partitionSequence);
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

					int[] sequence = new int[siteCount];
					// TODO: dirty solution for gaps when taxa between the tree topologies don't match
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
			Partition partition, //
			NodeRef node, //
			int categoryCount, //
			int stateCount, //
			BufferIndexHelper matrixBufferHelper, //
			BufferIndexHelper eigenBufferHelper //
	) {

		double[][] probabilities = new double[categoryCount][stateCount * stateCount];

		TreeModel treeModel = partition.treeModel;
		BranchRateModel branchRateModel = partition.branchRateModel;
		
		int nodeNum = node.getNumber();
		matrixBufferHelper.flipOffset(nodeNum);
		int branchIndex = nodeNum;//matrixBufferHelper.getOffsetIndex(nodeNum);
		
		double branchRate = branchRateModel.getBranchRate(treeModel, node);
		double branchTime = treeModel.getBranchLength(node) * branchRate;
		
        if (branchTime < 0.0) {
            throw new RuntimeException("Negative branch length: " + branchTime);
        }
		
        int count = 1;
        substitutionModelDelegate.updateTransitionMatrices(beagle, new int[] { branchIndex }, new double[] { branchTime }, count);
        
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
				array[i] = ((Codons) dataType).getState(sequence.getChar(k), sequence.getChar(k + 1), sequence.getChar(k + 2));
				k += 3;
			}// END: replications loop

		} else {

			for (int i = 0; i < siteCount; i++) {
				array[i] = dataType.getState(sequence.getChar(i));
			}// END: replications loop

		}// END: dataType check

		return array;
	}// END: sequence2intArray

} // END: class
