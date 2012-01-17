/*
 * EpochBranchSubstitutionModel.java
 *
 * Copyright (c) 2002-2012 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

///////////////////
//---SCRAPBOOK---//
///////////////////

//		timeSlices = new double[nModels + 1];
//		timeSlices[0] = nodeHeight;
//		for (int k = 1; k < nModels; k++) {
//
//			timeSlices[k] = transitionTimes[k - 1];
//
//		}
//		timeSlices[nModels] = parentHeight;
//
//		printArray(timeSlices);

//		timeSlices = new double[nModels + 1];
//		timeSlices[0] = parentHeight;
//		int k = 1;
//		for (int i = transitionTimes.length - 1; i >= 0; i--) {
//		
//			timeSlices[k] = transitionTimes[i];
//			k++;
//			
//		}
//		timeSlices[nModels] = nodeHeight;
//		
//		for (int i = transitionTimes.length - 0; i >= 0; i--) {
//		
//			System.out.println(timeSlices[i] - timeSlices[i+1]);
//		
//		}
//		System.out.println();

// int nodeNum = node.getNumber();
//
// switch (nodeNum) {
// case 0:
// returnValue = 2;
// weights = new double[] { epochTimes.getParameterValue(0),
// 73.7468 - epochTimes.getParameterValue(0) , 0.0};
// convolutionMatricesMap.put(bufferIndex, weights);
// // System.err.println("Case 0 gets index = " + bufferIndex);
// break;
// case 1:
// returnValue = 1;
// break;
// case 2:
// returnValue = 2;
// weights = new double[] { epochTimes.getParameterValue(0) - 10.0083,
// 45.2487 - epochTimes.getParameterValue(0) + 10.0083, 0.0 };
// convolutionMatricesMap.put(bufferIndex, weights);
// // System.err.println("Case 2 gets index = " + bufferIndex);
// break;
// case 3:
// returnValue = 1;
// break;
// }//END: hardcoded (sw/b)itch
//
// System.out.println("============");
// printArray(weights);
// System.out.println("============");

//		int returnValue = 0;
//		for (int i = 0; i < transitionTimes.length - 0; i++) {
//
//			double transitionTime = transitionTimes[i];
//			
//			if (nodeHeight < transitionTime && transitionTime <= parentHeight) {
//
//				// TODO gets proportional lengths
////				weights[i] = (transitionTimes[i] - nodeHeight);
////				weights[i + 1] = (parentHeight - transitionTimes[i]);
//				
//				returnValue = nModels;
//
//				System.out.println("branch of length " + branchLength
//						+ " intersects transition time " + i);
////				System.out.println("\t gets weights of " + weights[i]);
////				System.out.println("\t and " + weights[i + 1]);
//				System.out.println("\t returnValue = " + returnValue
//						+ " (branch needs convolution)");
//				System.out.println();
//
//		
//			} else if (nodeHeight > transitionTime) {
//				// TODO gets full length
//
////				weights[i] = 0.0;// branchLength;
//				returnValue = 1;
//
//				System.out.println("branch of length " + branchLength
//						+ " is above transition time " + i);
////				System.out.println("\t gets weight of " + weights[i]);
//				System.out.println("\t returnValue = " + returnValue);
//				System.out.println();
//
//			} else if (parentHeight < transitionTime) {
//				// TODO gets length 0
//
////				weights[i] = 0.0;
//				returnValue = 1;
//
//				System.out.println("branch of length " + branchLength
//						+ " is beneath transition time " + i);
////				System.out.println("\t gets weight of " + weights[i]);
//				System.out.println("\t returnValue = " + returnValue);
//				System.out.println();
//
//			} else {
//
//				System.err.println("bad juju");
//
//			}
//
//		}// END: transitionTimes loop	

package dr.app.beagle.evomodel.sitemodel;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import beagle.Beagle;
import dr.app.beagle.evomodel.substmodel.EigenDecomposition;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.treelikelihood.BufferIndexHelper;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

/**
 * @author Filip Bielejec
 * @author Marc A. Suchard
 * @version $Id$
 */
@SuppressWarnings("serial")
public class EpochBranchSubstitutionModel extends AbstractModel implements
		BranchSubstitutionModel, Citable {

	private final List<SubstitutionModel> substModelList;
	private final List<FrequencyModel> frequencyModelList;
	private final Parameter epochTimes;
	private int firstBuffer;
	private Map<Integer, double[]> convolutionMatricesMap = new HashMap<Integer, double[]>();

	public EpochBranchSubstitutionModel(List<SubstitutionModel> substModelList,
			List<FrequencyModel> frequencyModelList, Parameter epochTimes) {
		super("EpochBranchSubstitutionModel");

		// if (substModelList.size() != 2) {
		// throw new IllegalArgumentException(
		// "EpochBranchSubstitutionModel requires two SubstitutionModels");
		// }

		if (frequencyModelList.size() != 1) {
			throw new IllegalArgumentException(
					"EpochBranchSubstitutionModel requires one FrequencyModel");
		}

		this.substModelList = substModelList;
		this.frequencyModelList = frequencyModelList;
		this.epochTimes = epochTimes;

		for (SubstitutionModel model : substModelList) {
			addModel(model);
		}
		for (FrequencyModel model : frequencyModelList) {
			addModel(model);
		}

		addVariable(epochTimes);
	}

	public int getBranchIndex(final Tree tree, final NodeRef node,
			int bufferIndex) {
		// TODO return substModelList.size() if branch will require convolution
		// 3 cases: 0 - first transition, i-th to i+1-th transition, last
		// transition - infinity

		int nModels = substModelList.size();
		int lastTransitionTime = nModels - 2;

		double[] weights = new double[nModels + 0];
		double[] transitionTimes = epochTimes.getParameterValues();
		double parentHeight = tree.getNodeHeight(tree.getParent(node));
		double nodeHeight = tree.getNodeHeight(node);
		double branchLength = tree.getBranchLength(node);

		int returnValue = 0;

		///////////////////////
		//---TODO: TRIAL 2---//
		///////////////////////
		
//		//first case: 0th transition time
		if (parentHeight <= transitionTimes[0]) {

			weights[0] = branchLength;
			returnValue = 0;

		} else {

			if (nodeHeight < transitionTimes[0] && transitionTimes[0] <= parentHeight) {

				weights[0] = transitionTimes[0] - nodeHeight;
				returnValue = nModels;
				
			} else if (nodeHeight > transitionTimes[0]) {

				 weights[0] = 0;
				
			}// END: 0-th model overlap check

			//TODO: tu jest spierdolone
			//second case: i to i+1 transition times
			for (int i = 1; i <= lastTransitionTime; i++) {

				if (parentHeight <= transitionTimes[i]) {
				
					break;
					
				} else { 
					
					if (parentHeight >= transitionTimes[i] && transitionTimes[i] > nodeHeight) {
				
				    double startTime = Math.max(nodeHeight, transitionTimes[i - 1]);
	                double endTime = Math.min(parentHeight, transitionTimes[i]);
	                weights[i] = (endTime - startTime);
					
				} else if (nodeHeight > transitionTimes[i]) {
					
					//
					
				}
					
				}//END: bail out check
		
			}// END: i loop
			
			//third case: last transition time
			if (parentHeight >= transitionTimes[lastTransitionTime] && transitionTimes[lastTransitionTime] > nodeHeight) {
			
						weights[lastTransitionTime + 1] = parentHeight - transitionTimes[lastTransitionTime];
						returnValue = nModels;
			
					} else if (parentHeight <= transitionTimes[0]) {
			
						 weights[lastTransitionTime + 1] = 0;
			
					} else if (nodeHeight > transitionTimes[0]) {
			
						weights[lastTransitionTime + 1] = branchLength;
						returnValue = nModels-1;
						
					}// END: last transition time check
			
		}// END: if branch below first transition time bail out
		

		///////////////////////
		//---TODO: TRIAL 1---//
		///////////////////////
		
//		if (nodeHeight < transitionTimes[0] && transitionTimes[0] <= parentHeight) {
//
//			weights[0] = transitionTimes[0] - nodeHeight;
//			returnValue = nModels; // Branch overlaps 0th model
//
//		} else if (parentHeight <= transitionTimes[0]) {
//
//			weights[0] = branchLength;
//			returnValue = 0; // Whole branch in 0th model
//			// TODO Check this condition first and bail out of function if true
//
//		} else if (nodeHeight > transitionTimes[0]) {
//
//			// weights[0] = 0;
//
//		}// END: cutoff check
//
//		for (int i = 1; i <= lastTransitionTime; i++) {
//			if (nodeHeight <= transitionTimes[i]) {
//
//				weights[i] = transitionTimes[i] - transitionTimes[i - 1];
//
//			}// END: are we there yet check
//			
//			
//			
//		}// END: i loop
//
//		if (parentHeight >= transitionTimes[lastTransitionTime] && transitionTimes[lastTransitionTime] > nodeHeight) {
//
//			weights[lastTransitionTime + 1] = parentHeight - transitionTimes[lastTransitionTime];
////			returnValue = nModels;
//
//		} else if (parentHeight <= transitionTimes[0]) {
//
//			// weights[lastTransitionTime + 1] = 0;
//
//		} else if (nodeHeight > transitionTimes[0]) {
//
//			weights[lastTransitionTime + 1] = branchLength;
////			returnValue = nModels;
//			
//		}// END: last transition time check

		
		

		
		
		System.out.println("branch length: " + branchLength);
		System.out.println("return value: " + returnValue);
		printArray(weights);

		convolutionMatricesMap.put(bufferIndex, weights);

		return returnValue;
	}// END: getBranchIndex

	/*
	 * Return number of transition matrices buffers
	 */
	public int getExtraBufferCount() {
		return (2);
	}// END: getBufferCount

	@Override
	public void setFirstBuffer(int firstBufferCount) {
		firstBuffer = firstBufferCount;
	}// END: setFirstBuffer

	public EigenDecomposition getEigenDecomposition(int branchIndex,
			int categoryIndex) {
		return substModelList.get(branchIndex).getEigenDecomposition();
	}

	public SubstitutionModel getSubstitutionModel(int branchIndex,
			int categoryIndex) {
		return substModelList.get(branchIndex);
	}

	public double[] getStateFrequencies(int categoryIndex) {
		return frequencyModelList.get(categoryIndex).getFrequencies();
	}

	public boolean canReturnComplexDiagonalization() {
		for (SubstitutionModel model : substModelList) {
			if (model.canReturnComplexDiagonalization()) {
				return true;
			}
		}
		return false;
	}

	public int getEigenCount() {
		return substModelList.size() + 1; // Use an extra eigenIndex to identify
		// branches that need convolution
	}

	protected void handleModelChangedEvent(Model model, Object object, int index) {
		fireModelChanged();
	}

	@SuppressWarnings("unchecked")
	protected void handleVariableChangedEvent(Variable variable, int index,
			Parameter.ChangeType type) {
	}

	protected void storeState() {
	}

	protected void restoreState() {
	}

	protected void acceptState() {
	}

	public void updateTransitionMatrices(Beagle beagle, int eigenIndex,
			BufferIndexHelper bufferHelper, final int[] probabilityIndices,
			final int[] firstDerivativeIndices,
			final int[] secondDervativeIndices, final double[] edgeLengths,
			int count) {

		if (eigenIndex < substModelList.size()) {
			// Branch falls in a single category
			beagle.updateTransitionMatrices(bufferHelper
					.getOffsetIndex(eigenIndex), probabilityIndices,
					firstDerivativeIndices, secondDervativeIndices,
					edgeLengths, count);
		} else {
			// Branch requires convolution of two or more matrices

			// for(int i = 0; i < substModelList.size(); i++) {
			// int[] firstIndices = Arrays.copyOfRange(probabilityIndices,
			// eigenIndex, eigenIndex);
			// int[] secondIndices = Arrays.copyOfRange(probabilityIndices,
			// eigenIndex, eigenIndex);
			// int[] resultIndices = Arrays.copyOfRange(probabilityIndices,
			// eigenIndex, eigenIndex);
			// }//END: subst model loop

			// System.err.println("eigenIndex = " + eigenIndex);
			//			
			// double[] firstBranchLengths = { epochTimes.getParameterValue(0),
			// 73.7468 - epochTimes.getParameterValue(0) };
			// double[] secondBranchLengths = { epochTimes.getParameterValue(0)
			// - 10.0083, 45.2487
			// - epochTimes.getParameterValue(0) + 10.00832 };

			// double[] firstBranchLengths = { epochTimes.getParameterValue(0),
			// 73.7468 - epochTimes.getParameterValue(0) };
			// double[] secondBranchLengths = { epochTimes.getParameterValue(0)
			// - 10.0083, 45.2487
			// - epochTimes.getParameterValue(0) + 10.00832 };

			for (int i = 0; i < count; i++) {

				final int resultIndex = probabilityIndices[i];
				final int firstIndex = firstBuffer;
				final int secondIndex = firstBuffer + 1;

				// System.out.println("firstIndex " + firstIndex);
				// System.out.println("secondIndex " + secondIndex);

				// int[] firstProbIndices = { firstIndex };
				// double[] firstBranchLength = { 72 };
				// beagle.updateTransitionMatrices(bufferHelper.getOffsetIndex(1),
				// firstProbIndices, null, null, firstBranchLength, 1);
				// int[] secondProbIndices = { secondIndex };
				// double[] secondBranchLength = { 22 };
				// beagle.updateTransitionMatrices(bufferHelper.getOffsetIndex(0),
				// secondProbIndices, null, null, secondBranchLength, 1);
				// int[] resultProbIndices = { resultIndex };
				//				
				// beagle.convolveTransitionMatrices(firstProbIndices,
				// secondProbIndices, resultProbIndices, 1);
				// System.err.println("Writing to buffer " +
				// resultProbIndices[0]);

				double[] length = convolutionMatricesMap
						.get(probabilityIndices[i]);

				// System.err.println("Attempting pI = " +
				// probabilityIndices[i]);
				// System.err.println("Times: " + length[0] + " " + length[1]);

				int[] firstProbIndices = { firstIndex };
				beagle.updateTransitionMatrices(bufferHelper.getOffsetIndex(1),
						firstProbIndices, null, null,
						new double[] { length[1] }, 1);
				int[] secondProbIndices = { secondIndex };
				beagle.updateTransitionMatrices(bufferHelper.getOffsetIndex(0),
						secondProbIndices, null, null,
						new double[] { length[0] }, 1);
				int[] resultProbIndices = { resultIndex };

				beagle.convolveTransitionMatrices(firstProbIndices,
						secondProbIndices, resultProbIndices, 1);

				// System.err.println("Writing to buffer " +
				// resultProbIndices[0]);

				// System.out.println(Arrays.copyOfRange(firstBranchLengths, i,
				// i)[0]);

				// double[] tmp = new double[16*4];
				// beagle.getTransitionMatrix(firstIndex, tmp);
				// printMatrix(tmp);

				// System.out.println("first " + getSubstitutionModel(0,
				// 0).getVariable(0).getValue(0));
				// beagle.getTransitionMatrix(firstIndex, tmp);
				// printMatrix(tmp);
				//				
				// System.out.println("second " + getSubstitutionModel(1,
				// 0).getVariable(0).getValue(0));
				// beagle.getTransitionMatrix(secondIndex, tmp);
				// printMatrix(tmp);
				//				
				// System.out.println("third");
				// beagle.getTransitionMatrix(resultIndex, tmp);
				// printMatrix(tmp);

				// System.out.println(convolutionMatricesMap.get(firstIndex)[0]);
				// System.out.println(convolutionMatricesMap.get(firstIndex)[1]);
				// System.out.println(convolutionMatricesMap.get(secondIndex)[0]);
				// System.out.println(convolutionMatricesMap.get(secondIndex)[1]);

			}// END: count loop

		}// END: eigenIndex check
	}// END: updateTransitionMatrices

	/**
	 * @return a list of citations associated with this object
	 */
	public List<Citation> getCitations() {
		List<Citation> citations = new ArrayList<Citation>();
		citations.add(new Citation(new Author[] { new Author("F", "Bielejec"),
				new Author("P", "Lemey"), new Author("MA", "Suchard") },
				Citation.Status.IN_PREPARATION));
		return citations;
	}

	@Override
	public void setEigenDecomposition(Beagle beagle, int eigenIndex,
			BufferIndexHelper bufferHelper, int dummy) {

		if (eigenIndex < substModelList.size()) {

			EigenDecomposition ed = getEigenDecomposition(eigenIndex, dummy);

			beagle.setEigenDecomposition(bufferHelper
					.getOffsetIndex(eigenIndex), ed.getEigenVectors(), ed
					.getInverseEigenVectors(), ed.getEigenValues());

		}
	}

	public static void printMatrix(double[] matrix) {
		int stateCount = matrix.length;
		for (int row = 0; row < stateCount; row++) {
			System.out.print("| ");
			for (int col = 0; col < stateCount; col++)
				System.out.print(String.format("%.6f", matrix[col + row
						* stateCount])
						+ " ");
			System.out.print("|\n");
		}
		System.out.print("\n");

	}// END: printMatrix

	public static void printArray(double[] array) {
		for (int col = 0; col < array.length; col++) {
			System.out.println(String.format("%.6f", array[col]));
		}
		System.out.print("\n");
	}// END: printArray

}