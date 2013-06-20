/*
 * Partition.java
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

import beagle.Beagle;
import beagle.BeagleFactory;
import dr.app.beagle.evomodel.branchmodel.BranchModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.treelikelihood.BufferIndexHelper;
import dr.app.beagle.evomodel.treelikelihood.SubstitutionModelDelegate;
import dr.app.bss.Utils;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
import dr.evolution.util.Taxon;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.math.MathUtils;
import org.apache.commons.math.random.MersenneTwister;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class Partition {

	private static final boolean DEBUG = false;

	// Constructor fields
	public int from;
	public int to;
	public int every;
	private BranchModel branchModel;
	private TreeModel treeModel;
	private GammaSiteRateModel siteModel;
	private BranchRateModel branchRateModel;
	private FrequencyModel freqModel;

	// Buffer helpers
	private BufferIndexHelper partialBufferHelper;
	private BufferIndexHelper scaleBufferHelper;
	private BufferIndexHelper matrixBufferHelper;

	// Beagle stuff
	private Beagle beagle;
	private SubstitutionModelDelegate substitutionModelDelegate;

	// int fields
	private Integer partitionNumber;
	private int partitionSiteCount;
	private int nodeCount;
	private int tipCount;
	private int internalNodeCount;
	private int stateCount;
	private int compactPartialsCount;
	private int patternCount;
	private int categoryCount;

	// Sequence fields
	private Map<Taxon, int[]> sequenceList;
	private DataType dataType;
	private boolean hasAncestralSequence = false;
	private Sequence ancestralSequence = null;

	private MersenneTwister random;

	public Partition(TreeModel treeModel, //
			BranchModel branchModel, //
			GammaSiteRateModel siteModel, //
			BranchRateModel branchRateModel, //
			FrequencyModel freqModel, //
			int from, //
			int to, //
			int every //
	) {

		this.treeModel = treeModel;
		this.siteModel = siteModel;
		this.freqModel = freqModel;
		this.branchModel = branchModel;
		this.branchRateModel = branchRateModel;

		this.from = from;
		this.to = to;
		this.every = every;

		dataType = freqModel.getDataType();
		partitionSiteCount = getPartitionSiteCount();

		setBufferHelpers();
		setSubstitutionModelDelegate();
		loadBeagleInstance();

		sequenceList = new HashMap<Taxon, int[]>();
		random = new MersenneTwister(MathUtils.nextLong());

	}// END: Constructor

	private void setSubstitutionModelDelegate() {
		substitutionModelDelegate = new SubstitutionModelDelegate(treeModel,
				branchModel);
	}// END: setSubstitutionModelDelegate

	private void setBufferHelpers() {

		nodeCount = treeModel.getNodeCount();
		matrixBufferHelper = new BufferIndexHelper(nodeCount, 0);

		tipCount = treeModel.getExternalNodeCount();
		internalNodeCount = treeModel.getInternalNodeCount();

		partialBufferHelper = new BufferIndexHelper(nodeCount, tipCount);
		scaleBufferHelper = new BufferIndexHelper(internalNodeCount + 1, 0);

	}// END: setBufferHelpers

	public void loadBeagleInstance() {

		compactPartialsCount = tipCount;
		stateCount = dataType.getStateCount();
		patternCount = partitionSiteCount;
		categoryCount = siteModel.getCategoryCount();

		int[] resourceList = new int[] { 0 };
		long preferenceFlags = 0;
		long requirementFlags = 0;

		beagle = BeagleFactory.loadBeagleInstance(tipCount, //
				partialBufferHelper.getBufferCount(), //
				compactPartialsCount, //
				stateCount, //
				patternCount, //
				substitutionModelDelegate.getEigenBufferCount(), //
				substitutionModelDelegate.getMatrixBufferCount(), //
				categoryCount, //
				scaleBufferHelper.getBufferCount(), //
				resourceList, //
				preferenceFlags, //
				requirementFlags);

	}// END: loadBeagleInstance

	public void simulatePartition() {

		try {

			NodeRef root = treeModel.getRoot();

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

				category[i] = randomChoicePDF(categoryProbs, partitionNumber,
						"categories");

			}

			int[] parentSequence = new int[partitionSiteCount];

			// set ancestral sequence for partition if it exists
			if (hasAncestralSequence) {

				if (partitionSiteCount == ancestralSequence.getLength()) {

					parentSequence = sequence2intArray(ancestralSequence);

				} else {

					throw new RuntimeException("Ancestral sequence length of "
							+ ancestralSequence.getLength()
							+ " does not match partition site count of "
							+ partitionSiteCount + ".");

				}

			} else {

				double[] frequencies = freqModel.getFrequencies();
				for (int i = 0; i < partitionSiteCount; i++) {

					parentSequence[i] = randomChoicePDF(frequencies,
							partitionNumber, "root");

				}

			}// END:ancestralSequence check

			if (DEBUG) {
				synchronized (this) {

					System.out.println("parent Sequence:");
					Utils.printArray(parentSequence);
				}
			}

			substitutionModelDelegate.updateSubstitutionModels(beagle);

			int categoryCount = siteModel.getCategoryCount();
			traverse(root, parentSequence, category, categoryCount);

			if (DEBUG) {
				synchronized (this) {
					printSequences();
				}

			}

			beagle.finalize();
//			System.gc();
			
		} catch (Exception e) {
			e.printStackTrace();
		} catch (Throwable e) {
			System.err.println("BeagleException: " + e.getMessage());
			System.exit(-1);
		}

	}// END: simulatePartition

	private void traverse(NodeRef node, //
			int[] parentSequence, //
			int[] category, //
			int categoryCount //
	) {

		for (int iChild = 0; iChild < treeModel.getChildCount(node); iChild++) {

			NodeRef child = treeModel.getChild(node, iChild);
			int[] partitionSequence = new int[partitionSiteCount];
			double[] cProb = new double[stateCount];

			double[][] probabilities = getTransitionProbabilities(child,
					categoryCount, stateCount);

			// synchronized (this) {
			// Utils.print2DArray(probabilities);
			// }

			for (int i = 0; i < partitionSiteCount; i++) {

				System.arraycopy(probabilities[category[i]], parentSequence[i]
						* stateCount, cProb, 0, stateCount);

				// synchronized (this) {
				// Utils.printArray(cProb);
				// }

				partitionSequence[i] = randomChoicePDF(cProb, partitionNumber,
						"seq");

			}

			if (treeModel.getChildCount(child) == 0) {

				Taxon taxon = treeModel.getNodeTaxon(child);
				sequenceList.put(taxon, partitionSequence);

			} // END: tip node check

			traverse(treeModel.getChild(node, iChild), partitionSequence,
					category, categoryCount);

		}// END: child nodes loop

	}// END: traverse

	private int[] sequence2intArray(Sequence sequence) {

		int array[] = new int[partitionSiteCount];

		if (dataType instanceof Codons) {

			int k = 0;
			for (int i = 0; i < partitionSiteCount; i++) {
				array[i] = ((Codons) dataType).getState(sequence.getChar(k),
						sequence.getChar(k + 1), sequence.getChar(k + 2));
				k += 3;
			}// END: replications loop

		} else {

			for (int i = 0; i < partitionSiteCount; i++) {
				array[i] = dataType.getState(sequence.getChar(i));
			}// END: replications loop

		}// END: dataType check

		return array;
	}// END: sequence2intArray

	private double[][] getTransitionProbabilities(NodeRef node, //
			int categoryCount, //
			int stateCount //
	) {

		double[][] probabilities = new double[categoryCount][stateCount
				* stateCount];

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

	private double getTotal(double[] array, int start, int end) {
		double total = 0.0;
		for (int i = start; i < end; i++) {
			total += array[i];
		}
		return total;
	}// END: getTotal

	private int randomChoicePDF(double[] pdf, int partitionNumber, String error) {

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

		System.out.println(error);
		System.exit(-1);
		throw new RuntimeException(
				"randomChoiceUnnormalized falls through -- negative components in input distribution?");

	}// END: randomChoicePDF

	// /////////////
	// --SETTERS--//
	// /////////////

	public void setPartitionNumber(Integer partitionNumber) {
		this.partitionNumber = partitionNumber;
	}

	public void setAncestralSequence(Sequence ancestralSequence) {
		this.ancestralSequence = ancestralSequence;
		this.hasAncestralSequence = true;
	}// END: setAncestralSequence

	// /////////////
	// --GETTERS--//
	// /////////////

	public int getPartitionSiteCount() {
		return ((to - from) / every) + 1;
	}// END: getPartitionSiteCount

	public BranchModel getBranchModel() {
		return this.branchModel;
	}// END: getBranchModelic

	public FrequencyModel getFreqModel() {
		return freqModel;
	}// END: getFreqModel

	public Integer getPartitionNumber() {
		return partitionNumber;
	}// END: getPartitionNumber

	public DataType getDataType() {
		return dataType;
	}// END: getDataType

	public Map<Taxon, int[]> getSequencesMap() {
		return sequenceList;
	}// END: getSequenceList

	public Sequence getAncestralSequence() {
		return ancestralSequence;
	}
	
	// ///////////////
	// --DEBUGGING--//
	// ///////////////

	public void printSequences() {
		System.out.println("partition " + partitionNumber);
		Utils.printMap(sequenceList);
	}// END: printSequences

}// END: class
