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

package dr.app.beagle.evomodel.sitemodel;

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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;

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
	}// END: constructor

	/**
	 * @return number of extra transition matrices buffers to allocate
	 */
	public int getExtraBufferCount() {
		//TODO: count * 4
		return 100*4;
	}// END: getBufferCount

	public void setFirstBuffer(int firstBufferCount) {
		firstBuffer = firstBufferCount;
	}// END: setFirstBuffer

	public EigenDecomposition getEigenDecomposition(int branchIndex,
			int categoryIndex) {
		return substModelList.get(branchIndex).getEigenDecomposition();
	}// END: getEigenDecomposition

	public SubstitutionModel getSubstitutionModel(int branchIndex,
			int categoryIndex) {
		return substModelList.get(branchIndex);
	}// END: getSubstitutionModel

	public double[] getStateFrequencies(int categoryIndex) {
		return frequencyModelList.get(categoryIndex).getFrequencies();
	}// END: getStateFrequencies

	public int getEigenCount() {
		// Use an extra eigenIndex to identify branches that need convolution
		return substModelList.size() + 1;
	}// END: getEigenCount

	public void setEigenDecomposition(Beagle beagle, int eigenIndex,
			BufferIndexHelper bufferHelper, int dummy) {

		if (eigenIndex < substModelList.size()) {

			EigenDecomposition ed = getEigenDecomposition(eigenIndex, dummy);

			beagle.setEigenDecomposition(bufferHelper.getOffsetIndex(eigenIndex), 
					ed.getEigenVectors(), 
					ed.getInverseEigenVectors(), 
					ed.getEigenValues()
					);

		}// END: nModels check
	}// END: setEigenDecomposition

	public boolean canReturnComplexDiagonalization() {
		for (SubstitutionModel model : substModelList) {
			if (model.canReturnComplexDiagonalization()) {
				return true;
			}
		}
		return false;
	}// END: canReturnComplexDiagonalization

	protected void handleModelChangedEvent(Model model, Object object, int index) {
		fireModelChanged();
	}// END: handleModelChangedEvent

	@SuppressWarnings("unchecked")
	protected void handleVariableChangedEvent(Variable variable, int index,
			Parameter.ChangeType type) {
	}// END: handleVariableChangedEvent

	protected void storeState() {
	}// END: storeState

	protected void restoreState() {
	}// END: restoreState

	protected void acceptState() {
	}// END: acceptState

	public int getBranchIndex(final Tree tree, final NodeRef node,
			int bufferIndex) {

		int nModels = substModelList.size();
		int lastTransitionTime = nModels - 2;

		double[] weights = new double[nModels + 0];
		double[] transitionTimes = epochTimes.getParameterValues();
		double parentHeight = tree.getNodeHeight(tree.getParent(node));
		double nodeHeight = tree.getNodeHeight(node);
		double branchLength = tree.getBranchLength(node);

		int returnValue = 0;

		// TODO: simplify this logic, it's a mess
		if (parentHeight <= transitionTimes[0]) {

			weights[0] = branchLength;
			returnValue = 0;

		} else {

			// first case: 0th transition time
			if (nodeHeight < transitionTimes[0] && transitionTimes[0] <= parentHeight) {

				weights[0] = transitionTimes[0] - nodeHeight;
				returnValue = nModels;

			} else {

				weights[0] = 0;

			}// END: 0-th model check

			// second case: i to i+1 transition times
			for (int i = 1; i <= lastTransitionTime; i++) {

				if (nodeHeight < transitionTimes[i]) {

					if (parentHeight <= transitionTimes[i] && transitionTimes[i - 1] < nodeHeight) {

						weights[i] = branchLength;
						returnValue = i;

					} else {

						double startTime = Math.max(nodeHeight,
								transitionTimes[i - 1]);
						double endTime = Math.min(parentHeight,
								transitionTimes[i]);

						if (endTime < startTime) {

							weights[i] = 0;

						} else {

							weights[i] = (endTime - startTime);
							returnValue = nModels;

						}// END: negative weights check

					}// END: full branch in middle epoch check

				} else {

					weights[i] = 0;

				}// END: i-th model check

			}// END: i loop

			// third case: last transition time
			if (parentHeight >= transitionTimes[lastTransitionTime] && transitionTimes[lastTransitionTime] > nodeHeight) {

				weights[lastTransitionTime + 1] = parentHeight - transitionTimes[lastTransitionTime];
				returnValue = nModels;

			} else if (nodeHeight > transitionTimes[lastTransitionTime]) {

				weights[lastTransitionTime + 1] = branchLength;
				returnValue = nModels - 1;

			} else {

				weights[lastTransitionTime + 1] = 0;

			}// END: last transition time check

		}// END: if branch below first transition time bail out

//		 System.out.println("bufferIndex: " + bufferIndex);
//		 System.out.println("branch length: " + branchLength);
//		 System.out.println("return value: " + returnValue);
//		 printArray(weights, weights.length);

		convolutionMatricesMap.put(bufferIndex, weights);

		return returnValue;
	}// END: getBranchIndex

	public void updateTransitionMatrices(Beagle beagle, 
			int eigenIndex,
			BufferIndexHelper bufferHelper, 
			final int[] probabilityIndices,
			final int[] firstDerivativeIndices,
			final int[] secondDervativeIndices, 
			final double[] edgeLengths,
			int count // number of branches to update in paralell
	) {

		if (eigenIndex < substModelList.size()) {

//			System.out.println("\neigenIndex: " + bufferHelper.getOffsetIndex(eigenIndex) + "\n");
//			System.out.println("Populating indices: ");
//			printArray(probabilityIndices, probabilityIndices.length);
//			System.out.println("for edge lengths: ");
//			printArray(edgeLengths, edgeLengths.length);

			// Branches fall in a single category
			beagle.updateTransitionMatrices(bufferHelper.getOffsetIndex(eigenIndex), 
					probabilityIndices,
					firstDerivativeIndices, 
					secondDervativeIndices,
					edgeLengths, 
					count);

//			for(int k =0;k<probabilityIndices.length;k++){
//				
//				double tmp[] = new double[4 * 4 * 4];
//				beagle.getTransitionMatrix(probabilityIndices[k], // matrixIndex
//						tmp // outMatrix
//						);
//				
//				System.out.println(probabilityIndices[k]);
//				printMatrix(tmp, 4, 4);
//				}
			
			
		} else {

//			System.err.println(count);
//			Scanner sc = new Scanner(System.in);
//			System.out.println("Press Enter to continue");
//			sc.nextLine();
			
			// Branches require convolution of two or more matrices
			int[] firstBuffers = new int[count];
			int[] secondBuffers = new int[count];
			int[] firstExtraBuffers = new int[count];
			int[] secondExtraBuffers = new int[count];
			
			int[] resultBranchBuffers = new int[count];

			int[] probabilityBuffers = new int[count];
			int[] firstConvolutionBuffers = new int[count];
			int[] secondConvolutionBuffers = new int[count];
			int[] resultConvolutionBuffers = new int[count];
			
			for (int i = 0; i < count; i++) {

				firstBuffers[i] = firstBuffer + i;
				secondBuffers[i] = (firstBuffer + count) + i;
				firstExtraBuffers[i] = (firstBuffer + 2*count) + i;
				secondExtraBuffers[i] = (firstBuffer + 3*count) + i;
				
				resultBranchBuffers[i] = probabilityIndices[i];

			}// END: count loop

			for (int i = 0; i < substModelList.size(); i++) {

				int eigenBuffer = bufferHelper.getOffsetIndex(i);
				double[] weights = new double[count];
				
				System.out.println("\ni: " + i + ", eigenBuffer: " + eigenBuffer + "\n");

				for (int j = 0; j < count; j++) {

					int index = probabilityIndices[j];
					weights[j] = convolutionMatricesMap.get(index)[i];

				}// END: count loop


				if ((i == 1) && (i == (substModelList.size() - 1))) {
					
					probabilityBuffers = secondBuffers;
					
					firstConvolutionBuffers = firstBuffers;
					secondConvolutionBuffers = probabilityBuffers;
					resultConvolutionBuffers = resultBranchBuffers;

				} else if ((i == 1) && (i != (substModelList.size() - 1))) {

					probabilityBuffers = secondBuffers;
					
					firstConvolutionBuffers = firstBuffers;
					secondConvolutionBuffers = probabilityBuffers;
					resultConvolutionBuffers = firstExtraBuffers;
					
				} else if ((i != 1) && (i == (substModelList.size() - 1))) {
					
					// even
					if (i % 2 == 0) {

						probabilityBuffers = firstBuffers;

						firstConvolutionBuffers = firstExtraBuffers;
						secondConvolutionBuffers = probabilityBuffers;
						resultConvolutionBuffers = resultBranchBuffers;
						
						// odd
					} else {
						
						probabilityBuffers = secondBuffers;
						
						firstConvolutionBuffers = secondExtraBuffers;
						secondConvolutionBuffers = probabilityBuffers;
						resultConvolutionBuffers = resultBranchBuffers;
					}
					

				} else {

					// even
					if (i % 2 == 0) {

						probabilityBuffers = firstBuffers;

						firstConvolutionBuffers = firstExtraBuffers;
						secondConvolutionBuffers = probabilityBuffers;
						resultConvolutionBuffers = secondExtraBuffers;
						
						// odd
					} else {
						
						probabilityBuffers = secondBuffers;
						
						firstConvolutionBuffers = secondExtraBuffers;
						secondConvolutionBuffers = probabilityBuffers;
						resultConvolutionBuffers = firstExtraBuffers;

					}// END: even-odd check

				}// END: first-last buffer check

				System.out.println("Populating buffers: ");
				printArray(probabilityBuffers, probabilityBuffers.length);
				System.out.println("for weights: ");
				printArray(weights, weights.length);
				
				beagle.updateTransitionMatrices(eigenBuffer, // eigenIndex
						probabilityBuffers, // probabilityIndices
						null, // firstDerivativeIndices
						null, // secondDerivativeIndices
						weights, // edgeLengths
						count // count
						);
				
				if (i != 0) {

					System.out.println("convolving buffers: ");
					printArray(firstConvolutionBuffers, firstConvolutionBuffers.length);
					System.out.println("with buffers: ");
					printArray(secondConvolutionBuffers, secondConvolutionBuffers.length);
					System.out.println("into buffers: ");
					printArray(resultConvolutionBuffers, resultConvolutionBuffers.length);
					
					beagle.convolveTransitionMatrices(firstConvolutionBuffers, // A
							secondConvolutionBuffers, // B
							resultConvolutionBuffers, // C
							count // count
							);
					
				}// END: 0 eigen index check

			}// END: eigen indices loop

		}// END: eigenIndex check
	}// END: updateTransitionMatrices

	/**
	 * @return a list of citations associated with this object
	 */
	public List<Citation> getCitations() {
		List<Citation> citations = new ArrayList<Citation>();
		citations.add(new Citation(new Author[] { new Author("F", "Bielejec"),
				new Author("P", "Lemey"), new Author("G", "Baele"),
				new Author("MA", "Suchard") }, Citation.Status.IN_PREPARATION));
		return citations;
	}// END: getCitations

	// /////////////
	// ---DEBUG---//
	// /////////////

	public static void printMatrix(double[] matrix, int nrow, int ncol) {
		for (int row = 0; row < nrow; row++) {
//			System.out.print("| ");
			for (int col = 0; col < nrow; col++)
				System.out.print(String
						.format("%.20f", matrix[col + row * nrow])
						+ ", ");
			System.out.print("\n");
//			System.out.print("|\n");
		}
		System.out.print("\n");
	}// END: printMatrix

	public static void printMatrix(int[] matrix, int nrow, int ncol) {
		for (int row = 0; row < nrow; row++) {
			System.out.print("| ");
			for (int col = 0; col < nrow; col++)
				System.out.print(matrix[col + row * nrow] + " ");
			System.out.print("|\n");
		}
		System.out.print("\n");
	}// END: printMatrix

	public static void printArray(double[] array, int ncol) {
		for (int col = 0; col < ncol; col++) {
			System.out.println(String.format("%.4f", array[col]));
		}
		System.out.print("\n");
	}// END: printArray

	public static void printArray(int[] array, int nrow) {
		for (int col = 0; col < nrow; col++) {
			System.out.println(array[col]);
		}
		System.out.print("\n");
	}// END: printArray

	// //////////////////
	// ---END: DEBUG---//
	// //////////////////

}// END: class