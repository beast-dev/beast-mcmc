/*
 * Partition.java
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

package dr.app.beagle.tools;

import java.util.LinkedHashMap;
import java.util.Map;

import dr.evomodel.treedatalikelihood.BufferIndexHelper;
import org.apache.commons.math.random.MersenneTwister;

import beagle.Beagle;
import beagle.BeagleFactory;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.treelikelihood.SubstitutionModelDelegate;
import dr.app.bss.Utils;
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
 */
public class Partition {

	private static final boolean DEBUG = false;

	// Constructor fields
	public int from;
	public int to;
	public int every;
	private BranchModel branchModel;
	private TreeModel treeModel;
	private GammaSiteRateModel siteRateModel;
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
	private LinkedHashMap<Taxon, int[]> alignmentMap;
	//	private LinkedHashMap<NodeRef, int[]> sequencesMap = new LinkedHashMap<NodeRef, int[]>();
	private DataType dataType;
	private boolean hasRootSequence = false;
	private Sequence rootSequence = null;
	private boolean outputAncestralSequences = false;

	// Random number generation
	private MersenneTwister random;

	// Annotating trees
//	private boolean annotateTree = true;

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
		this.siteRateModel = siteModel;
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

		alignmentMap = new LinkedHashMap<Taxon, int[]>();
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
		siteRateCategoryCount = siteRateModel.getCategoryCount();

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
			double[] categoryRates = siteRateModel.getCategoryRates();
			beagle.setCategoryRates(categoryRates);

			// probabilities for gamma category rates
			double[] categoryProbs = siteRateModel.getCategoryProportions();
//			beagle.setCategoryWeights(0, categoryProbs);

//			Utils.printArray(categoryRates);
//			Utils.printArray(categoryProbs);

			int[] category = new int[partitionSiteCount];
			for (int i = 0; i < partitionSiteCount; i++) {

				category[i] = randomChoicePDF(categoryProbs, partitionNumber,
						"categories");

			}

//			category = new int[] {1, 0, 0, 0, 0, 1, 0, 1, 0, 0 };

			if(DEBUG){
				System.out.println("category for each site:");
				Utils.printArray(category);
			}//END: DEBUG

			int[] parentSequence = new int[partitionSiteCount];

			// set ancestral sequence for partition if it exists
			if (hasRootSequence) {

				if (rootSequence.getLength() == partitionSiteCount) {

					parentSequence = sequence2intArray(rootSequence);

				} else if (dataType instanceof Codons && rootSequence.getLength() == 3 * partitionSiteCount) {

					parentSequence = sequence2intArray(rootSequence);

				} else {

					throw new RuntimeException("Ancestral sequence length of "
							+ rootSequence.getLength()
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
					System.out.println();
					System.out.println("root Sequence:");
					Utils.printArray(parentSequence);
				}
			}//END: DEBUG

			substitutionModelDelegate.updateSubstitutionModels(beagle);

			traverse(root, parentSequence, category);

			if (DEBUG) {
				synchronized (this) {
					System.out.println("Simulated alignment:");
					printSequences();
				}
			}//END: DEBUG

			beagle.finalize();

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

		if (DEBUG) {
			synchronized (this) {
				System.out.println();
				System.out.println("I'm at: " + node.toString());
				System.out.println();
			}
		}//END: DEBUG

		for (int iChild = 0; iChild < treeModel.getChildCount(node); iChild++) {

			NodeRef child = treeModel.getChild(node, iChild);
			int[] partitionSequence = new int[partitionSiteCount];
			double[] cProb = new double[stateCount];

			double[][] probabilities = getTransitionProbabilities(child);

			if (DEBUG) {
				synchronized (this) {
					System.out.println("Going to child " + iChild + ": " + child.toString());
					System.out.println("Child finite transition probs matrix:");
					Utils.print2DArray(probabilities, stateCount);
					System.out.println();
				}
			}// END: DEBUG

			for (int i = 0; i < partitionSiteCount; i++) {

                if (DEBUG) {
                    System.out.println("length of copy: " + stateCount);
                    System.out.println("source length: " + probabilities[category[i]].length);
                    System.out.println("destination length: " + cProb.length);
                }

				System.arraycopy(probabilities[category[i]], parentSequence[i] * stateCount, cProb, 0, stateCount);

				if (DEBUG) {
					synchronized (this) {
						System.out.println("site:" + i);
						System.out.println("site probs:");
						Utils.printArray(cProb);
					}
				}// END: DEBUG

				partitionSequence[i] = randomChoicePDF(cProb, partitionNumber,
						"seq");

			}// END: i loop

			if (DEBUG) {
				synchronized (this) {
//					partitionSequence = new int[]{1, 3, 2, 3, 0, 1, 0, 1, 0, 2, 2, 0, 1, 3, 3, 3, 0, 1, 2, 1, 3, 1, 1, 1, 1, 3, 0, 0, 3, 2, 3, 2, 3, 2, 1, 2, 1, 3, 2, 3, 3, 0, 2, 2, 3, 2, 3, 2, 3, 1, 2, 0, 2, 1, 3, 2, 3, 1, 1, 1, 1, 0, 2, 3, 1, 0, 2, 1, 2, 1, 3, 0, 0, 0, 0, 0, 2, 0, 2, 3, 1, 0, 1, 3, 0, 2, 1, 2, 1, 3, 0, 0, 3, 2, 2, 0, 1, 0, 0, 3  };
					System.out.println("Simulated sequence:");
					Utils.printArray(partitionSequence);
				}
			}// END: if DEBUG

//			if(annotateTree) {
//				
//				sequencesMap.put(child, partitionSequence);
//				
//			}

			if (treeModel.getChildCount(child) == 0) {

				Taxon taxon = treeModel.getNodeTaxon(child);
				alignmentMap.put(taxon, partitionSequence);

				if (DEBUG) {
					synchronized (this) {
						System.out.println("Simulated sequence (translated):");
						System.out.println(Utils.intArray2Sequence(taxon, partitionSequence, BeagleSequenceSimulator.gapFlag, dataType).getSequenceString());
					}
				}// END: DEBUG

			} else {

				if(outputAncestralSequences) {

					alignmentMap.put(new Taxon("internalNodeHeight" + treeModel.getNodeHeight(child)), partitionSequence);

				}

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
		double branchTime = branchLength * branchRate;// * siteRate;

		if (DEBUG) {
			synchronized (this) {
				System.out.println("Branch length: " + branchLength
						+ " branch rate: " + branchRate + " branch time: "
						+ branchTime);// + " site rate: " + siteRate);
			}
		}// END: DEBUG

		int count = 1;
		substitutionModelDelegate.updateTransitionMatrices(beagle,
				new int[] { branchIndex }, new double[] { branchTime }, count);

		double transitionMatrix[] = new double[siteRateCategoryCount
				* stateCount * stateCount];

		beagle.getTransitionMatrix(branchIndex, //
				transitionMatrix //
		);

		for (int siteRateCat = 0; siteRateCat < siteRateCategoryCount; siteRateCat++) {

			System.arraycopy(transitionMatrix, siteRateCat * stateCount * stateCount,
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

	public void setRootSequence(Sequence rootSequence) {
		this.rootSequence = rootSequence;
		this.hasRootSequence = true;
	}// END: setAncestralSequence

	public void setOutputAncestralSequences(boolean outputAncestralSequences) {
		this.outputAncestralSequences = outputAncestralSequences;
	}

	// /////////////
	// --GETTERS--//
	// /////////////

//	public boolean isOutputAncestralSequences() {
//		return outputAncestralSequences;
//	}

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

	public Map<Taxon, int[]> getTaxonSequencesMap() {
		return alignmentMap;
	}// END: getSequenceList

//	public LinkedHashMap<NodeRef, int[]> getSequenceMap() {
//		return sequencesMap;
//	}

	public Sequence getRootSequence() {
		return rootSequence;
	}

	// ///////////////
	// --DEBUGGING--//
	// ///////////////

	public void printSequences() {
		System.out.println("partition " + partitionNumber);
		Utils.printMap(alignmentMap);
	}// END: printSequences

}// END: class