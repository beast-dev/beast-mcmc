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
//                            returnValue = i;
							returnValue = nModels;

                        }//END: negative weights check

                    }//END: full branch in middle epoch check

                } else {

                    weights[i] = 0;

                }// END: i-th model check

            }//END: i loop

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

//        System.out.println("branch length: " + branchLength);
//        System.out.println("return value: " + returnValue);
//        printArray(weights, weights.length);

        convolutionMatricesMap.put(bufferIndex, weights);

        return returnValue;
    }// END: getBranchIndex

    /*
      * Return number of transition matrices buffers
      */
    public int getExtraBufferCount() {
        return (2);
    }// END: getBufferCount

    //	@Override
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
                                         int count //number of branches to update in paralell
                                         ) {
	
        if (eigenIndex < substModelList.size()) {

            // Branches fall in a single category
            beagle.updateTransitionMatrices(bufferHelper.getOffsetIndex(eigenIndex),
                    probabilityIndices,
                    firstDerivativeIndices,
                    secondDervativeIndices,
                    edgeLengths,
                    count);
        } else {
        	
			// Branches require convolution of two or more matrices
			for (int i = 0; i < count; i++) {

				boolean wasZero = false;

				int index = probabilityIndices[i];
				double[] weights = convolutionMatricesMap.get(index);
				int[] resultProbIndices = { index };
				int[] firstProbIndices = { firstBuffer };
				int[] secondProbIndices = { firstBuffer + 1 };

				System.out.println("branch: " + index);
				System.out.println("weights:");
				printArray(weights, weights.length);

				for (int j = 0; j < weights.length - 1; j++) {

					if (weights[j] != 0.0 && !wasZero) {

						System.out.println("populating matrix index " + firstProbIndices[0] + " for length " + weights[j + 1]);

						beagle.updateTransitionMatrices(bufferHelper
								.getOffsetIndex(j + 1), // eigenIndex
								firstProbIndices, // probabilityIndices
								null, // firstDerivativeIndices
								null, // secondDerivativeIndices
								new double[] { weights[j + 1] }, // edgeLengths
								1 // count
								);

						if (j == 0) {

							System.out.println("populating matrix index " + secondProbIndices[0] + " for length " + weights[j]);

							beagle.updateTransitionMatrices(bufferHelper
									.getOffsetIndex(j), // eigenIndex
									secondProbIndices, // probabilityIndices
									null, // firstDerivativeIndices
									null, // secondDerivativeIndices
									new double[] { weights[j] }, // edgeLengths
									1 // count
									);

						} else {

							System.out.println("copying matrix index " + resultProbIndices[0] + " into matrix index " + secondProbIndices[0]);

							// TODO: check if this works as expected
							// memcpy resultIndex into secondIndex
							double tmp[] = new double[4 * 4 * 1];
							beagle.getTransitionMatrix(resultProbIndices[0], // matrixIndex
									tmp // outMatrix
									);

							beagle.setTransitionMatrix(secondProbIndices[0], // matrixIndex
									tmp, // inMatrix
									1.0); // paddedValue
						}

						System.out.println("convolving matrix index " + firstProbIndices[0] + " with " + secondProbIndices[0] + " into " + resultProbIndices[0]);

						beagle.convolveTransitionMatrices(firstProbIndices, // A
								secondProbIndices, // B
								resultProbIndices, // C
								1 // count
								);

					} else if (weights[j + 1] != 0 && !wasZero) {

						System.out.println("zero weight detected");
						wasZero = true;

						System.out.println("populating matrix index " + firstProbIndices[0] + " for length " + weights[j + 2]);

						beagle.updateTransitionMatrices(bufferHelper
								.getOffsetIndex(j + 1), // eigenIndex
								firstProbIndices, // probabilityIndices
								null, // firstDerivativeIndices
								null, // secondDerivativeIndices
								new double[] { weights[j + 2] }, // edgeLengths
								1 // count
								);

						System.out.println("populating matrix index " + secondProbIndices[0] + " for length " + weights[j + 1]);

						beagle.updateTransitionMatrices(bufferHelper
								.getOffsetIndex(j), // eigenIndex
								secondProbIndices, // probabilityIndices
								null, // firstDerivativeIndices
								null, // secondDerivativeIndices
								new double[] { weights[j + 1] }, // edgeLengths
								1 // count
								);

						System.out.println("convolving matrix index " + firstProbIndices[0] + " with " + secondProbIndices[0] + " into " + resultProbIndices[0]);

						beagle.convolveTransitionMatrices(firstProbIndices, // A
								secondProbIndices, // B
								resultProbIndices, // C
								1 // count
								);

						j++;
					} else if (wasZero) {

						// break from loop when next weights are zero
						System.out.println("Done here. Bailing out");
						break;

					}// END: weight check
				}// END: weight loop

				System.out.println();
			}// END: branches loop
        	
        	
        	
//             // Branch requires convolution of two or more matrices
//            for (int i = 0; i < count; i++) {
//
//                final int resultIndex = probabilityIndices[i];
//                final int firstIndex = firstBuffer;
//                final int secondIndex = firstBuffer + 1;
//
//                double[] weights = convolutionMatricesMap.get(probabilityIndices[i]);
//
//                 
//                int[] firstProbIndices = { firstIndex };
//                beagle.updateTransitionMatrices(bufferHelper.getOffsetIndex(1), // eigenIndex
//                                                firstProbIndices, // probabilityIndices
//                                                null, // firstDerivativeIndices
//                                                null, // secondDerivativeIndices
//                                                new double[]{ weights[1] }, // edgeLengths
//                                                1 // count
//                                                );
//                
//                int[] secondProbIndices = { secondIndex };
//                beagle.updateTransitionMatrices(bufferHelper.getOffsetIndex(0), // eigenIndex
//                                                secondProbIndices, // probabilityIndices
//                                                null, // firstDerivativeIndices
//                                                null, // secondDerivativeIndices
//                                                new double[]{ weights[0] }, // edgeLengths
//                                                1 // count
//                                                );
//                
//                int[] resultProbIndices = { resultIndex };
//                beagle.convolveTransitionMatrices(firstProbIndices, // A
//                                                  secondProbIndices, // B
//                                                  resultProbIndices, // C
//                                                  1 // count
//                                                  );
//			      
//            }// END: count loop

        }// END: eigenIndex check
    }// END: updateTransitionMatrices

    /**
     * @return a list of citations associated with this object
     */
    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(new Citation(new Author[]{new Author("F", "Bielejec"),
                new Author("P", "Lemey"), new Author("G", "Baele"),
                new Author("MA", "Suchard")}, Citation.Status.IN_PREPARATION));
        return citations;
    }

    //	@Override
    public void setEigenDecomposition(Beagle beagle, int eigenIndex,
                                      BufferIndexHelper bufferHelper, int dummy) {

        if (eigenIndex < substModelList.size()) {

            EigenDecomposition ed = getEigenDecomposition(eigenIndex, dummy);

            beagle.setEigenDecomposition(bufferHelper
                    .getOffsetIndex(eigenIndex), ed.getEigenVectors(), ed
                    .getInverseEigenVectors(), ed.getEigenValues());

        }
    }

    public static void printMatrix(double[] matrix, int nrow, int ncol) {
        for (int row = 0; row < nrow; row++) {
            System.out.print("| ");
            for (int col = 0; col < nrow; col++)
                System.out.print(String.format("%.20f", matrix[col + row
                        * nrow])
                        + " ");
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
    
    public static void printArray(int[] array) {
        for (int col = 0; col < array.length; col++) {
            System.out.println(array[col]);
        }
        System.out.print("\n");
    }// END: printArray

}// END: class