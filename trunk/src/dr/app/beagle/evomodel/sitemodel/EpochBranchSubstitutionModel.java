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
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.*;

/**
 * @author Filip Bielejec
 * @author Marc A. Suchard
 * @version $Id$
 */
@SuppressWarnings("serial")
public class EpochBranchSubstitutionModel extends AbstractModel implements
        BranchSubstitutionModel, Citable {

    public static final boolean TRY_EPOCH = true;

    private final List<SubstitutionModel> substModelList;
    private final List<FrequencyModel> frequencyModelList;
    private final Parameter epochTimes;
    private int firstBuffer;
    private Map<Integer, double[]> convolutionMatricesMap = new HashMap<Integer, double[]>();
    private int requestedBuffers;

    public EpochBranchSubstitutionModel(List<SubstitutionModel> substModelList,
                                        List<FrequencyModel> frequencyModelList, Parameter epochTimes) {

        super("EpochBranchSubstitutionModel");

        if (frequencyModelList.size() != 1) {
            throw new IllegalArgumentException(
                    "EpochBranchSubstitutionModel requires one FrequencyModel");
        }

        this.substModelList = substModelList;
        this.frequencyModelList = frequencyModelList;
        this.epochTimes = epochTimes;
        this.requestedBuffers = 0;

        for (SubstitutionModel model : substModelList) {
            addModel(model);
        }

        for (FrequencyModel model : frequencyModelList) {
            addModel(model);
        }

        addVariable(epochTimes);
    }// END: Constructor

    /**
     * @return number of extra transition matrices buffers to allocate
     */
    public int getExtraBufferCount(TreeModel treeModel) {

//        // loop over the tree to determine the count
//        double[] transitionTimes = epochTimes.getParameterValues();
//        int rootId = treeModel.getRoot().getNumber();
//        int count = 0;
//        for (NodeRef node : treeModel.getNodes()) {
//
//            if (node.getNumber() != rootId) {
//
//                double nodeHeight = treeModel.getNodeHeight(node);
//                double parentHeight = treeModel.getNodeHeight(treeModel
//                        .getParent(node));
//
//                for (int i = 0; i < transitionTimes.length; i++) {
//
//                    if (nodeHeight <= transitionTimes[i] && transitionTimes[i] < parentHeight) {
//                        count++;
//                        break;
//                    }// END: transition time check check
//
//                }// END: transition times loop
//            }// END: root check
//        }// END: nodes loop

        requestedBuffers = 20;
        
//        System.out.println("fixed count = " + count);
//        System.out.println("extra buffers = " + requestedBuffers);
        
        return requestedBuffers;
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

    /**
     * Calculate weights that branch spends in each substitution model
     *
     * @param tree
     * @param node
     * @return nModels if branch needs convolution, subst model index if not
     */
    public int getBranchIndex(final Tree tree, final NodeRef node,
                              int bufferIndex) {

        int nModels = substModelList.size();
        int lastTransitionTime = nModels - 2;

        double[] weights = new double[nModels];
        double[] transitionTimes = epochTimes.getParameterValues();
        double parentHeight = tree.getNodeHeight(tree.getParent(node));
        double nodeHeight = tree.getNodeHeight(node);
        double branchLength = tree.getBranchLength(node);

        int returnValue = 0;

        if (parentHeight <= transitionTimes[0]) {

            weights[0] = branchLength;
            returnValue = 0;

        } else {

            // first case: 0-th transition time
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

        convolutionMatricesMap.put(bufferIndex, weights);

//        System.out.println("bufferIndex: " + bufferIndex);
//        printArray(weights, weights.length);
        
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
        	
            // Branches fall in a single category
            beagle.updateTransitionMatrices(bufferHelper.getOffsetIndex(eigenIndex),
                    probabilityIndices,
                    firstDerivativeIndices,
                    secondDervativeIndices,
                    edgeLengths,
                    count);

            // ////////////////////////////////////////////////////

//			for (int k = 0; k < probabilityIndices.length; k++) {
//
//				double tmp[] = new double[4 * 4 * 4];
//				beagle.getTransitionMatrix(probabilityIndices[k], // matrixIndex
//						tmp // outMatrix
//						);
//
//				System.out.println(probabilityIndices[k]);
//				printMatrix(tmp, 4, 4);
//			}

            // ////////////////////////////////////////////////////

        } else {

            // Branches require convolution of two or more matrices
        	int stepSize = requestedBuffers/4 ;
        	
          System.out.println("stepSize: " + stepSize);
          System.out.println("count from tree = " + count);
//          System.out.println("convolutionMatricesMap.size() = " + convolutionMatricesMap.size());
//        	System.out.println("probabilityIndices ");
//        	printArray(probabilityIndices, probabilityIndices.length);
          
        	int step = 0;
        	while(step < count) {

				System.out.println("step: " + step);
        		
            int[] firstBuffers = new int[stepSize];
            int[] secondBuffers = new int[stepSize];
            int[] firstExtraBuffers = new int[stepSize];
            int[] secondExtraBuffers = new int[stepSize];

            int[] resultBranchBuffers = new int[stepSize];

            int[] probabilityBuffers = new int[stepSize];
            int[] firstConvolutionBuffers = new int[stepSize];
            int[] secondConvolutionBuffers = new int[stepSize];
            int[] resultConvolutionBuffers = new int[stepSize];
            
            for (int i = 0; i < stepSize; i++) {  

                firstBuffers[i] = firstBuffer + i;
                secondBuffers[i] = (firstBuffer + stepSize) + i;
                firstExtraBuffers[i] = (firstBuffer + 2 * stepSize) + i;
                secondExtraBuffers[i] = (firstBuffer + 3 * stepSize) + i;

					if (i < count) {
						resultBranchBuffers[i] = probabilityIndices[i + step];
					}
                
            }// END: stepSize loop

        	////////////////////////////////////////////////////////////
            
//        	System.out.println("resultBranchBuffers ");
//        	printArray(resultBranchBuffers, resultBranchBuffers.length);
        	
        	////////////////////////////////////////////////////////////

				for (int i = 0; i < substModelList.size(); i++) {

					int eigenBuffer = bufferHelper.getOffsetIndex(i);
					double[] weights = new double[stepSize];

					for (int j = 0; j < stepSize; j++) {

						if ((step + j) < count) {

							int index = probabilityIndices[j + step];
//							System.out.println("step + j: " + (step + j) + " index: " + index);
							weights[j] = convolutionMatricesMap.get(index)[i];

						}// END: index padding check

					}// END: stepSize loop

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

                checkBuffers(probabilityBuffers);
                int operationsCount = Math.min(stepSize, (count - step));
                
//              System.out.println("eigenBuffer: " + eigenBuffer);                
                System.out.println("Populating buffers: ");
                printArray(probabilityBuffers, operationsCount);
                System.out.println("for weights: ");
                printArray(weights, operationsCount);
                
                beagle.updateTransitionMatrices(eigenBuffer, // eigenIndex
                        probabilityBuffers, // probabilityIndices
                        null, // firstDerivativeIndices
                        null, // secondDerivativeIndices
                        weights, // edgeLengths
                        operationsCount // count
                );

                if (i != 0) {

					System.out.println("convolving buffers: ");
					printArray(firstConvolutionBuffers, operationsCount);
					System.out.println("with buffers: ");
					printArray(secondConvolutionBuffers, operationsCount);
					System.out.println("into buffers: ");
					printArray(resultConvolutionBuffers, operationsCount);    
					
                    beagle.convolveTransitionMatrices(firstConvolutionBuffers, // A
                            secondConvolutionBuffers, // B
                            resultConvolutionBuffers, // C
                            operationsCount // count
                    );

					}// END: 0-th eigen index check

				}// END: eigen indices loop        	
        	
         	step += stepSize;
        	}// END: step loop

//        	System.exit(-1);
        	System.out.println("=========================================================");
		}// END: eigenIndex check

        // ////////////////////////////////////////////////////

//		for (int k = 0; k < probabilityIndices.length; k++) {
//
//			double tmp[] = new double[4 * 4 * 4];
//			beagle.getTransitionMatrix(probabilityIndices[k], // matrixIndex
//					tmp // outMatrix
//					);
//
//			System.out.println(probabilityIndices[k]);
//			printMatrix(tmp, 4, 4);
//		}

        // ////////////////////////////////////////////////////


    }// END: updateTransitionMatrices

    private void checkBuffers(int[] probabilityBuffers) {
        for (int buffer : probabilityBuffers) {
            if (buffer >= firstBuffer + requestedBuffers) {
                System.err.println("Programming error: requesting use of BEAGLE transition matrix buffer not allocated.");
                System.err.println("Allocated: 0 to " + (firstBuffer + requestedBuffers - 1));
                System.err.println("Requested = " + buffer);
                System.err.println("Please complain to Button-Boy");
                
//                System.exit(-1);
                
            }
        }
    }//END: checkBuffers
    
	private int countFilledBuffers(int[] buffers) {
		int nonZero = 0;
		for (int i = 0; i < buffers.length; i++) {
			if (buffers[i] != 0) {
				nonZero++;
			}
		}

		return nonZero;
	}// END: countFilledBuffers
	
	private int countTrueFilledBuffers(int[] buffers) {
		int zero = 0;
		for (int i = buffers.length - 1; i >= 0; i--) {

			if (buffers[i] == 0) {
				zero++;
			} else {
				break;
			}

		}

		return buffers.length - zero;
	}

	private int countFilledWeights(double[] weights) {
		int nonZero = 0;
		for (int i = 0; i < weights.length; i++) {
			if (weights[i] != 0.0) {
				nonZero++;
			}
		}

		return nonZero;
	}// END: countFilledWeights
    
    /**
     * @return a list of citations associated with this object
     */
    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(new Citation(new Author[]{new Author("F", "Bielejec"),
                new Author("P", "Lemey"), new Author("G", "Baele"),
                new Author("MA", "Suchard")}, Citation.Status.IN_PREPARATION));
        return citations;
    }// END: getCitations

    // /////////////
    // ---DEBUG---//
    // /////////////

    public static void printMatrix(double[][] matrix, int nrow, int ncol) {
        for (int row = 0; row < nrow; row++) {
            for (int col = 0; col < nrow; col++)
                System.out.print(String.format(Locale.US, "%.20f", matrix[col + row * nrow]) + " ");
            System.out.print("\n");
        }
        System.out.print("\n");
    }// END: printMatrix

    public static void printMatrix(double[] matrix, int nrow, int ncol) {
        for (int row = 0; row < nrow; row++) {
            System.out.print("| ");
            for (int col = 0; col < nrow; col++)
                System.out.print(String.format(Locale.US, "%.20f", matrix[col + row * nrow]) + " ");
//			System.out.print("\n");
            System.out.print("|\n");
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
            System.out.println(String.format(Locale.US, "%.20f", array[col]));
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