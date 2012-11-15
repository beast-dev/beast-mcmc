/*
 * EpochBranchSubstitutionModel.java
 *
 * Copyright (C) 2002-2012 Alexei Drummond, Andrew Rambaut & Marc A. Suchard
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
import dr.evomodel.branchratemodel.BranchRateModel;
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

@Deprecated // Switching to BranchModel
public class EpochBranchSubstitutionModel extends AbstractModel implements
        BranchSubstitutionModel, Citable {

    // /////////////
    // ---DEBUG---//
    // /////////////
	
	private static final boolean DEBUG_EPOCH = false;
	private static Integer stateCount = null;
	private static Integer categoryCount = null;
	
    // //////////////////
    // ---END: DEBUG---//
    // //////////////////
	
    public static final boolean TRY_EPOCH = true;
    public static final String EPOCH_BRANCH_SUBSTITUTION_MODEL = "EpochBranchSubstitutionModel";
    
    private final List<SubstitutionModel> substModelList;
    private final List<FrequencyModel> frequencyModelList;
    private final BranchRateModel branchRateModel;
    private final Parameter epochTimes;
    private int firstBuffer;
    private Map<Integer, double[]> convolutionMatricesMap = new HashMap<Integer, double[]>();
    private int requestedBuffers;

    public EpochBranchSubstitutionModel(List<SubstitutionModel> substModelList,
                                        List<FrequencyModel> frequencyModelList, 
                                        BranchRateModel branchRateModel, 
                                        Parameter epochTimes) {

        super(EPOCH_BRANCH_SUBSTITUTION_MODEL);

        if (frequencyModelList.size() != 1) {
            throw new IllegalArgumentException(
                    "EpochBranchSubstitutionModel requires one FrequencyModel");
        }

        this.substModelList = substModelList;
        this.frequencyModelList = frequencyModelList;
        this.epochTimes = epochTimes;
        this.requestedBuffers = 0;
        this.branchRateModel = branchRateModel;
        
        for (SubstitutionModel model : substModelList) {
            addModel(model);
        }

        for (FrequencyModel model : frequencyModelList) {
            addModel(model);
        }

        
        if(DEBUG_EPOCH) {
        	
        	stateCount = frequencyModelList.get(0).getDataType().getStateCount();
        	categoryCount = 4;
        	
        }//END: DEBUG_EPOCH
        
        addVariable(epochTimes);
    }// END: Constructor

    /**
     * @return number of extra transition matrices buffers to allocate
     */
    public int getExtraBufferCount(TreeModel treeModel) {
    	
        requestedBuffers = 100;
        System.out.println("Allocating " + requestedBuffers + " extra buffers.");
        
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

    @SuppressWarnings("rawtypes")
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

                        double startTime = Math.max(nodeHeight, transitionTimes[i - 1]);
                        double endTime = Math.min(parentHeight, transitionTimes[i]);

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

		if (branchRateModel != null) {
			weights = scaleArray(weights, branchRateModel.getBranchRate(tree, node));
		}

        convolutionMatricesMap.put(bufferIndex, weights);

		if (DEBUG_EPOCH) {

			System.out.println("bufferIndex: " + bufferIndex);
			System.out.println("weights: ");
			printArray(weights, weights.length);

		}// END: DEBUG_EPOCH
    		
        return returnValue;
    }// END: getBranchIndex

    public void updateTransitionMatrices(Beagle beagle,
                                         int eigenIndex,
                                         BufferIndexHelper bufferHelper,
                                         final int[] probabilityIndices,
                                         final int[] firstDerivativeIndices,
                                         final int[] secondDervativeIndices,
                                         final double[] edgeLengths,
                                         int count // number of branches to update in parallel
    ) {

        if (eigenIndex < substModelList.size()) {

			if (DEBUG_EPOCH) {

				System.out.println("Branch falls in a single category");
				System.out.println("eigenIndex: " + eigenIndex);
				System.out.println("Populating buffers: ");
				printArray(probabilityIndices, count);
				System.out.println("for weights: ");
				printArray(edgeLengths, count);

			}//END: DEBUG_EPOCH
        
            // Branches fall in a single category
            beagle.updateTransitionMatrices(bufferHelper.getOffsetIndex(eigenIndex),
                    probabilityIndices,
                    firstDerivativeIndices,
                    secondDervativeIndices,
                    edgeLengths,
                    count);

			if (DEBUG_EPOCH) {

				System.out.println("Transition probabilities from model: ");
				
				for (int k = 0; k < probabilityIndices.length; k++) {

					double tmp[] = new double[categoryCount * stateCount * stateCount];
					beagle.getTransitionMatrix(probabilityIndices[k], // matrixIndex
							tmp // outMatrix
					);

					System.out.println(probabilityIndices[k]);
					printMatrix(tmp, stateCount, stateCount);
				}

			}// END: DEBUG_EPOCH
            	
        } else {

            // Branches require convolution of two or more matrices
        	int stepSize = requestedBuffers/4 ;

			if (DEBUG_EPOCH) {

				System.out.println("Branch requires convolution");
				System.out.println("stepSize: " + stepSize);
//				System.out.println("count from tree = " + count);
//				System.out.println("convolutionMatricesMap.size() = " + convolutionMatricesMap.size());
				System.out.println("probabilityIndices: ");
				printArray(probabilityIndices, probabilityIndices.length);

			}//END: DEBUG_EPOCH
        		
        	int step = 0;
        	while(step < count) {

				if (DEBUG_EPOCH) {

					System.out.println("step: " + step);

				}//END: DEBUG_EPOCH
        		
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

				if (DEBUG_EPOCH) {

					System.out.println("resultBranchBuffers ");
					printArray(resultBranchBuffers, resultBranchBuffers.length);

				}//END: DEBUG_EPOCH

				for (int i = 0; i < substModelList.size(); i++) {

					int eigenBuffer = bufferHelper.getOffsetIndex(i);
					double[] weights = new double[stepSize];

					for (int j = 0; j < stepSize; j++) {

						if ((step + j) < count) {

							int index = probabilityIndices[j + step];
							
							if (DEBUG_EPOCH) {
							
							System.out.println("step + j: " + (step + j) + " index: " + index);

							}
							
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
  
					if (DEBUG_EPOCH) {

						System.out.println("eigenBuffer: " + eigenBuffer);
						System.out.println("Populating buffers: ");
						printArray(probabilityBuffers, operationsCount);
						System.out.println("for weights: ");
						printArray(weights, operationsCount);

					}//END: DEBUG_EPOCH
            		
                beagle.updateTransitionMatrices(eigenBuffer, // eigenIndex
                        probabilityBuffers, // probabilityIndices
                        null, // firstDerivativeIndices
                        null, // secondDerivativeIndices
                        weights, // edgeLengths
                        operationsCount // count
                );

                if (i != 0) {

                	if (DEBUG_EPOCH) {
                	
					System.out.println("convolving buffers: ");
					printArray(firstConvolutionBuffers, operationsCount);
					System.out.println("with buffers: ");
					printArray(secondConvolutionBuffers, operationsCount);
					System.out.println("into buffers: ");
					printArray(resultConvolutionBuffers, operationsCount);    

                	}//END: DEBUG_EPOCH
                	
                    beagle.convolveTransitionMatrices(firstConvolutionBuffers, // A
                            secondConvolutionBuffers, // B
                            resultConvolutionBuffers, // C
                            operationsCount // count
                    );

					}// END: 0-th eigen index check

				}// END: eigen indices loop        	
        	
         	step += stepSize;
        	}// END: step loop

		}// END: eigenIndex check

    	if (DEBUG_EPOCH) {
        
		System.out.println("Transition probabilities from model:");    		
		
		for (int k = 0; k < probabilityIndices.length; k++) {

			double tmp[] = new double[categoryCount * stateCount * stateCount];
			beagle.getTransitionMatrix(probabilityIndices[k], // matrixIndex
					tmp // outMatrix
					);

			System.out.println(probabilityIndices[k]);
			printMatrix(tmp, stateCount, stateCount);
		}

    	}//END: DEBUG_EPOCH

    }// END: updateTransitionMatrices

    private void checkBuffers(int[] probabilityBuffers) {
        for (int buffer : probabilityBuffers) {
            if (buffer >= firstBuffer + requestedBuffers) {
                System.err.println("Programming error: requesting use of BEAGLE transition matrix buffer not allocated.");
                System.err.println("Allocated: 0 to " + (firstBuffer + requestedBuffers - 1));
                System.err.println("Requested = " + buffer);
                System.err.println("Please complain to Button-Boy");
            }
        }
    }//END: checkBuffers
    
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
    
	public static void printArray(double[] array) {
		for (int i = 0; i < array.length; i++) {
			System.out.println(String.format(Locale.US, "%.10f", array[i]));
		}
		System.out.print("\n");
	}// END: printArray
    
	public static void printArray(int[] array) {
		for (int i = 0; i < array.length; i++) {
			System.out.println(array[i]);
		}
	}// END: printArray

    public static void printArray(double[] array, int nrow) {
        for (int row = 0; row < nrow; row++) {
            System.out.println(String.format(Locale.US, "%.10f", array[row]));
        }
        System.out.print("\n");
    }// END: printArray

    public static void printArray(int[] array, int nrow) {
        for (int row = 0; row < nrow; row++) {
            System.out.println(array[row]);
        }
        System.out.print("\n");
    }// END: printArray
	
	public static void print2DArray(double[][] array) {
		for (int row = 0; row < array.length; row++) {
			System.out.print("| ");
			for (int col = 0; col < array[row].length; col++) {
				System.out.print(String.format(Locale.US, "%.10f", array[row][col]) + " ");
			}
			System.out.print("|\n");
		}
		System.out.print("\n");
	}// END: print2DArray

	public static void print2DArray(int[][] array) {
		for (int row = 0; row < array.length; row++) {
			for (int col = 0; col < array[row].length; col++) {
				System.out.print(array[row][col] + " ");
			}
			System.out.print("\n");
		}
	}// END: print2DArray
    
    public static void printMatrix(double[][] matrix, int nrow, int ncol) {
        for (int row = 0; row < nrow; row++) {
            for (int col = 0; col < nrow; col++)
                System.out.print(String.format(Locale.US, "%.10f", matrix[col + row * nrow]) + " ");
            System.out.print("\n");
        }
        System.out.print("\n");
    }// END: printMatrix

    public static void printMatrix(double[] matrix, int nrow, int ncol) {
        for (int row = 0; row < nrow; row++) {
            System.out.print("| ");
            for (int col = 0; col < nrow; col++)
                System.out.print(String.format(Locale.US, "%.10f", matrix[col + row * nrow]) + " ");
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

	public double[] scaleArray(double[] array, double scalar) {

		for (int i = 0; i < array.length; i++) {
			array[i] = array[i] * scalar;
		}
		return array;

	}// END: scaleArray
    
    // //////////////////
    // ---END: DEBUG---//
    // //////////////////

    public int getExtraBufferCount_old(TreeModel treeModel) {

        // loop over the tree to determine the count
        double[] transitionTimes = epochTimes.getParameterValues();
        int rootId = treeModel.getRoot().getNumber();
		int count = 0;
		
        for (NodeRef node : treeModel.getNodes()) {

            if (node.getNumber() != rootId) {

                double nodeHeight = treeModel.getNodeHeight(node);
                double parentHeight = treeModel.getNodeHeight(treeModel
                        .getParent(node));

                for (int i = 0; i < transitionTimes.length; i++) {

                    if (nodeHeight <= transitionTimes[i] && transitionTimes[i] < parentHeight) {
                        count++;
                        break;
                    }// END: transition time check check

                }// END: transition times loop
            }// END: root check
        }// END: nodes loop

        requestedBuffers = count * 4;
        
        System.out.println("Allocating " + requestedBuffers + " extra buffers.");
        
        return requestedBuffers;
    }// END: getBufferCount_old
	
}// END: class
