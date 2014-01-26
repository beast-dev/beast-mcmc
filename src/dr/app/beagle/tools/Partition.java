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
	private int siteRateCategoryCount;

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
		siteRateCategoryCount = siteModel.getCategoryCount();

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
				siteRateCategoryCount, //
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

				if (ancestralSequence.getLength() == partitionSiteCount) {

					parentSequence = sequence2intArray(ancestralSequence);

				} else if (dataType instanceof Codons && ancestralSequence.getLength() == 3 * partitionSiteCount) {	
					
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

			traverse(root, parentSequence, category);

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
			int[] category //
	) {

		for (int iChild = 0; iChild < treeModel.getChildCount(node); iChild++) {

			NodeRef child = treeModel.getChild(node, iChild);
			int[] partitionSequence = new int[partitionSiteCount];
			double[] cProb = new double[stateCount];

			//TODO: make global, write for all categories
			double[][] probabilities = getTransitionProbabilities(child);

			if (DEBUG) {
				synchronized (this) {
					Utils.print2DArray(probabilities);
				}
			}// END: if DEBUG
			
			for (int i = 0; i < partitionSiteCount; i++) {

				System.arraycopy(probabilities[category[i]], parentSequence[i]
						* stateCount, cProb, 0, stateCount);

//				if (DEBUG) {
//					synchronized (this) {
//						Utils.printArray(cProb);
//					}
//				}// END: if DEBUG
				
				partitionSequence[i] = randomChoicePDF(cProb, partitionNumber,
						"seq");

			}

			if (treeModel.getChildCount(child) == 0) {

				Taxon taxon = treeModel.getNodeTaxon(child);
				sequenceList.put(taxon, partitionSequence);

			} // END: tip node check

			traverse(treeModel.getChild(node, iChild), partitionSequence,
					category);

		}// END: child nodes loop

	}// END: traverse

	private double[][] getTransitionProbabilities(NodeRef node //
	) {

		double[][] probabilities = new double[siteRateCategoryCount][stateCount
				* stateCount];

		int nodeNum = node.getNumber();
		matrixBufferHelper.flipOffset(nodeNum);
		int branchIndex = nodeNum;

		double branchRate = branchRateModel.getBranchRate(treeModel, node);
		double branchLength = treeModel.getBranchLength(node);

//		if (branchLength * branchRate < 0.0) {
//			throw new RuntimeException("Negative branch length: " + branchLength * branchRate);
//		}

		for (int siteRateCat = 0; siteRateCat < siteRateCategoryCount; siteRateCat++) {

			double siteRate = siteModel.getRateForCategory(siteRateCat);
            double branchTime = branchLength * branchRate * siteRate;			
			
			int count = 1;
			substitutionModelDelegate.updateTransitionMatrices(beagle,
					new int[] { branchIndex }, new double[] { branchTime },
					count);

			double transitionMatrix[] = new double[siteRateCategoryCount * stateCount * stateCount];

			beagle.getTransitionMatrix(branchIndex, //
					transitionMatrix //
			);

			System.arraycopy(transitionMatrix, siteRateCat * stateCount,
					probabilities[siteRateCat], 0, stateCount * stateCount);

		}// END: i loop

		return probabilities;
	}// END: getTransitionProbabilities
	
	// ///////////////////////////
	// ---START: EXPERIMENTAL---//
	// ///////////////////////////
	
//	void traverse(NodeRef node, //
//			int [] parentSequence, // 
//			int [] category //
//			) {
//		
//		for (int iChild = 0; iChild < treeModel.getChildCount(node); iChild++) {
//			
//			NodeRef child = treeModel.getChild(node, iChild);
//			
//            for (int i = 0; i < categoryCount; i++) {
//            	getTransitionProbabilities(child, i, probabilities[i]);
//            }
//
////            if(DEBUG){
////            Utils.print2DArray(probabilities);
////            }
//            
//        	int [] seq = new int[partitionSiteCount];
//    		double [] cProb = new double[stateCount];
//    		
//        	for (int i  = 0; i < partitionSiteCount; i++) {
//        		System.arraycopy(probabilities[category[i]], parentSequence[i] * stateCount, cProb, 0, stateCount);
//            	seq[i] = MathUtils.randomChoicePDF(cProb);
//        	}
//
//            if (treeModel.getChildCount(child) == 0) {
//            	Taxon taxon = treeModel.getNodeTaxon(child);
//            	sequenceList.put(taxon, seq);
//            }
//            
//			traverse(treeModel.getChild(node, iChild), seq, category);
//		
//		}//END: child loop
//	} // traverse
//	
//    void getTransitionProbabilities(
//NodeRef node, int rateCategory, double[] probs) {
//
//        NodeRef parent = treeModel.getParent(node);
//
//         double branchRate = branchRateModel.getBranchRate(treeModel, node);
//
//        // Get the operational time of the branch
//         double branchTime = branchRate * (treeModel.getNodeHeight(parent) - treeModel.getNodeHeight(node));
//
//        if (branchTime < 0.0) {
//            throw new RuntimeException("Negative branch length: " + branchTime);
//        }
//
//        double branchLength = siteModel.getRateForCategory(rateCategory) * branchTime;
//
//        siteModel.getSubstitutionModel().getTransitionProbabilities(branchLength, probs);
//    } // getTransitionProbabilities
	
	// /////////////////////////
	// ---END: EXPERIMENTAL---//
	// /////////////////////////
	
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

	private int randomChoicePDF(double[] pdf, int partitionNumber, String error) {

		int samplePos = -Integer.MAX_VALUE;
		double cumProb = 0.0;
		double u = random.nextDouble();
		
		for (int i = 0; i < pdf.length; i++) {
			
			cumProb += pdf[i];
			
			if (u < cumProb) {
				samplePos = i;
				break;
			}
		}

		return samplePos;
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

	public TreeModel getTreeModel() {
		return treeModel;
	}
	
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