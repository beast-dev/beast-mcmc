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

		int nModels = substModelList.size();
		int lastTransitionTime = nModels - 2;

		double[] weights = new double[nModels + 0];
		double[] transitionTimes = epochTimes.getParameterValues();
		double parentHeight = tree.getNodeHeight(tree.getParent(node));
		double nodeHeight = tree.getNodeHeight(node);
		double branchLength = tree.getBranchLength(node);

		int returnValue = 0;

		// /////////////////////
		// ---TODO: TRIAL 2---//
		// /////////////////////
		// TODO: simplify this logic, it's a mess

		// first case: 0th transition time
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

			// second case: i to i+1 transition times
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

				}// END: bail out check

			}// END: i loop

			// third case: last transition time
			if (parentHeight >= transitionTimes[lastTransitionTime] && transitionTimes[lastTransitionTime] > nodeHeight) {

				weights[lastTransitionTime + 1] = parentHeight - transitionTimes[lastTransitionTime];
				returnValue = nModels;

			} else if (parentHeight <= transitionTimes[lastTransitionTime]) {

				weights[lastTransitionTime + 1] = 0;

			} else if (nodeHeight > transitionTimes[lastTransitionTime]) {

				weights[lastTransitionTime + 1] = branchLength;
				returnValue = nModels - 1;

			}// END: last transition time check

		}// END: if branch below first transition time bail out

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
		// Use an extra eigenIndex to identify
		// branches that need convolution
		return substModelList.size() + 1; 
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
			beagle.updateTransitionMatrices(bufferHelper.getOffsetIndex(eigenIndex),
					                        probabilityIndices,
					                        firstDerivativeIndices, 
					                        secondDervativeIndices,
					                        edgeLengths, 
					                        count);
		} else {

			// Branch requires convolution of two or more matrices
			for (int i = 0; i < count; i++) {

				final int resultIndex = probabilityIndices[i];
				final int firstIndex = firstBuffer;
				final int secondIndex = firstBuffer + 1;

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
				new Author("P", "Lemey"), new Author("G", "Baele"),
				new Author("MA", "Suchard") }, Citation.Status.IN_PREPARATION));
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
		int stateCount = 4;// matrix.length;
		for (int row = 0; row < stateCount; row++) {
			System.out.print("| ");
			for (int col = 0; col < stateCount; col++)
				System.out.print(String.format("%.20f", matrix[col + row
						* stateCount])
						+ " ");
			System.out.print("|\n");
		}
		System.out.print("\n");

	}// END: printMatrix

	public static void printArray(double[] array) {
		for (int col = 0; col < array.length; col++) {
			System.out.println(String.format("%.20f", array[col]));
		}
		System.out.print("\n");
	}// END: printArray

}// END: class