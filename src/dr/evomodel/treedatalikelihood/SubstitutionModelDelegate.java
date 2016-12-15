/*
 * SubstitutionModelDelegate.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood;

import beagle.Beagle;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evolution.tree.Tree;
import dr.util.Timer;

import java.io.Serializable;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Filip Bielejec
 * @author Marc A. Suchard
 * @version $Id$
 */
public final class SubstitutionModelDelegate implements EvolutionaryProcessDelegate, Serializable {
    private static final boolean DEBUG = false;
    private static final boolean RUN_IN_SERIES = false;
    public static final boolean MEASURE_RUN_TIME = false;

    public double updateTime;
    public double convolveTime;

    private static final int BUFFER_POOL_SIZE_DEFAULT = 100;

    private final Tree tree;
    private final List<SubstitutionModel> substitutionModelList;
    private final BranchModel branchModel;

    private final int eigenCount;
    private final int nodeCount;

    private final int extraBufferCount;
    private final int reserveBufferIndex;

    private final BufferIndexHelper eigenBufferHelper;
    private final BufferIndexHelper matrixBufferHelper;

    private Deque<Integer> availableBuffers = new ArrayDeque<Integer>();

    /**
     * A class which handles substitution models including epoch models where multiple
     * substitution models on a branch are convolved.
     * @param tree
     * @param branchModel Describes which substitution models use on each branch
     */
    public SubstitutionModelDelegate(Tree tree, BranchModel branchModel) {
        this(tree, branchModel, 0, BUFFER_POOL_SIZE_DEFAULT);
    }

    /**
     * A class which handles substitution models including epoch models where multiple
     * substitution models on a branch are convolved.
     * @param tree
     * @param branchModel Describes which substitution models use on each branch
     * @param partitionNumber which data partition is this (used to offset eigen and matrix buffer numbers)
     */
    public SubstitutionModelDelegate(Tree tree, BranchModel branchModel, int partitionNumber) {
        this(tree, branchModel, partitionNumber, BUFFER_POOL_SIZE_DEFAULT);
    }

    public SubstitutionModelDelegate(Tree tree, BranchModel branchModel, int partitionNumber, int bufferPoolSize) {

        if (MEASURE_RUN_TIME) {
            updateTime = 0;
            convolveTime = 0;
        }

        this.tree = tree;

        this.substitutionModelList = branchModel.getSubstitutionModels();

        this.branchModel = branchModel;

        eigenCount = substitutionModelList.size();
        nodeCount = tree.getNodeCount();

        // two eigen buffers for each decomposition for store and restore.
        eigenBufferHelper = new BufferIndexHelper(eigenCount, 0, partitionNumber);

        // two matrices for each node less the root
        matrixBufferHelper = new BufferIndexHelper(nodeCount, 0, partitionNumber);

        this.extraBufferCount = branchModel.requiresMatrixConvolution() ?
                (bufferPoolSize > 0 ? bufferPoolSize : BUFFER_POOL_SIZE_DEFAULT) : 0;

        if (branchModel.requiresMatrixConvolution() && this.extraBufferCount < eigenCount) {
            throw new RuntimeException("SubstitutionModelDelegate requires at least " + eigenCount + " extra buffers to convolve matrices");
        }

        for (int i = 0; i < extraBufferCount; i++) {
            pushAvailableBuffer(i + matrixBufferHelper.getBufferCount());
        }

        // one extra created as a reserve
        // which is used to free up buffers when the avail stack is empty.
        reserveBufferIndex = matrixBufferHelper.getBufferCount() + extraBufferCount;

        if (DEBUG) {
            System.out.println("Creating reserve buffer with index: "
                    + reserveBufferIndex);
        }

    }// END: Constructor

    @Override
    public boolean canReturnComplexDiagonalization() {
        for (SubstitutionModel model : substitutionModelList) {
            if (model.canReturnComplexDiagonalization()) {
                return true;
            }
        }
        return false;
    }

    @Override
    public int getEigenBufferCount() {
        return eigenBufferHelper.getBufferCount();
    }

    @Override
    public int getMatrixBufferCount() {
        // plus one for the reserve buffer
        return matrixBufferHelper.getBufferCount() + extraBufferCount + 1;
    }

    @Override
    public int getSubstitutionModelCount() {
        return substitutionModelList.size();
    }

    @Override
    public SubstitutionModel getSubstitutionModel(int index) {
        return substitutionModelList.get(index);
    }

    @Override
    public int getEigenIndex(int bufferIndex) {
        return eigenBufferHelper.getOffsetIndex(bufferIndex);
    }

    @Override
    public int getMatrixIndex(int branchIndex) {
        return matrixBufferHelper.getOffsetIndex(branchIndex);
    }

    @Override
    public double[] getRootStateFrequencies() {
        return branchModel.getRootFrequencyModel().getFrequencies();
    }// END: getStateFrequencies


    @Override
    public void updateSubstitutionModels(Beagle beagle, boolean flipBuffers) {
        for (int i = 0; i < eigenCount; i++) {
            if (flipBuffers) {
                eigenBufferHelper.flipOffset(i);
            }

            EigenDecomposition ed = substitutionModelList.get(i).getEigenDecomposition();

            beagle.setEigenDecomposition(
                    eigenBufferHelper.getOffsetIndex(i),
                    ed.getEigenVectors(),
                    ed.getInverseEigenVectors(),
                    ed.getEigenValues());
        }
    }

    @Override
    public void updateTransitionMatrices(Beagle beagle, int[] branchIndices, double[] edgeLength, int updateCount, boolean flipBuffers) {

        int[][] probabilityIndices = new int[eigenCount][updateCount];
        double[][] edgeLengths = new double[eigenCount][updateCount];

        int[] counts = new int[eigenCount];

        List<Deque<Integer>> convolutionList = new ArrayList<Deque<Integer>>();

        for (int i = 0; i < updateCount; i++) {

            BranchModel.Mapping mapping = branchModel.getBranchModelMapping(tree.getNode(branchIndices[i]));
            int[] order = mapping.getOrder();
            double[] weights = mapping.getWeights();

            if (order.length == 1) {
                int k = order[0];
                if (flipBuffers) {
                    matrixBufferHelper.flipOffset(branchIndices[i]);
                }
                probabilityIndices[k][counts[k]] = matrixBufferHelper.getOffsetIndex(branchIndices[i]);
                edgeLengths[k][counts[k]] = edgeLength[i];
                counts[k]++;
            } else {
                double sum = 0.0;
                for (double w : weights) {
                    sum += w;
                }

                if (getAvailableBufferCount() < order.length) {
                    // too few buffers available, process what we have and continue...
                    if (flipBuffers) { throw new UnsupportedOperationException("flipping not implemented for Epoch models"); }
                    computeTransitionMatrices(beagle, probabilityIndices, edgeLengths, counts);
                    convolveMatrices(beagle, convolutionList);

                    // reset the counts
                    for (int k = 0; k < eigenCount; k++) {
                        counts[k] = 0;
                    }
                }

                Deque<Integer> bufferIndices = new ArrayDeque<Integer>();
                for (int j = 0; j < order.length; j++) {

                    int buffer = popAvailableBuffer();

                    if (buffer < 0) {
                        // no buffers available
                        throw new RuntimeException("Ran out of buffers for transition matrices - computing current list.");
                    }

                    int k = order[j];
                    probabilityIndices[k][counts[k]] = buffer;
                    edgeLengths[k][counts[k]] = weights[j] * edgeLength[i] / sum;
                    counts[k]++;

                    bufferIndices.add(buffer);
                }
                bufferIndices.add(matrixBufferHelper.getOffsetIndex(branchIndices[i]));

                convolutionList.add(bufferIndices);
            }// END: if convolution needed

        }// END: i loop

        computeTransitionMatrices(beagle, probabilityIndices, edgeLengths, counts);
        convolveMatrices(beagle, convolutionList);

    }// END: updateTransitionMatrices

    @Override
    public void flipTransitionMatrices(int[] branchIndices, int updateCount) {
        for (int i = 0; i < updateCount; i++) {
            matrixBufferHelper.flipOffset(branchIndices[i]);
        }
    }

    private void computeTransitionMatrices(Beagle beagle, int[][] probabilityIndices, double[][] edgeLengths, int[] counts) {

        Timer timer;
        if (MEASURE_RUN_TIME) {
            timer = new Timer();
            timer.start();
        }

        if (DEBUG) {
            System.out.print("Computing matrices:");
        }

        for (int i = 0; i < eigenCount; i++) {
            if (DEBUG) {
                for (int j = 0; j < counts[i]; j++) {
//                    System.out.print(" " + probabilityIndices[i][j]);
                    System.out.print(" " + probabilityIndices[i][j] + " (" + edgeLengths[i][j] + ")");
                }
            }
            if (counts[i] > 0) {
                beagle.updateTransitionMatrices(eigenBufferHelper.getOffsetIndex(i),
                        probabilityIndices[i],
                        null, // firstDerivativeIndices
                        null, // secondDerivativeIndices
                        edgeLengths[i],
                        counts[i]);
            }
        }

        if (DEBUG) {
            System.out.println();
        }

        if (MEASURE_RUN_TIME) {
            timer.stop();
            double timeInSeconds = timer.toSeconds();
            updateTime += timeInSeconds;
        }

    }//END: computeTransitionMatrices

    private void convolveMatrices(Beagle beagle, List<Deque<Integer>> convolutionList) {

        Timer timer;
        if (MEASURE_RUN_TIME) {
            timer = new Timer();
            timer.start();
        }

        while (convolutionList.size() > 0) {
            int[] firstConvolutionBuffers = new int[nodeCount];
            int[] secondConvolutionBuffers = new int[nodeCount];
            int[] resultConvolutionBuffers = new int[nodeCount];
            int operationsCount = 0;

            List<Deque<Integer>> empty = new ArrayList<Deque<Integer>>();

            for (Deque<Integer> convolve : convolutionList) {

                if (convolve.size() > 3) {
                    firstConvolutionBuffers[operationsCount] = convolve.pop();
                    secondConvolutionBuffers[operationsCount] = convolve.pop();

                    int buffer;
                    boolean done;

                    do {
                        done = true;

                        buffer = popAvailableBuffer();

                        if (buffer < 0) {
                            // no buffers available
//                        throw new RuntimeException("All out of buffers");

                            // we have run out of buffers, process what we have and continue...
                            if (DEBUG) {
                                System.out.println("Ran out of buffers for convolving - computing current list.");
                                System.out.print("Convolving " + operationsCount + " matrices:");
                                for (int i = 0; i < operationsCount; i++) {
                                    System.out.print(" " + firstConvolutionBuffers[i] + "*" + secondConvolutionBuffers[i] + "->" + resultConvolutionBuffers[i]);
                                }
                                System.out.println();
                            }

                            if (operationsCount > 0) {

                                convolveAndRelease(beagle, firstConvolutionBuffers, secondConvolutionBuffers, resultConvolutionBuffers, operationsCount);

                                // copy the uncompleted operation back down to the beginning of the operations list
                                firstConvolutionBuffers[0] = firstConvolutionBuffers[operationsCount];
                                secondConvolutionBuffers[0] = secondConvolutionBuffers[operationsCount];

                                // reset the operation count
                                operationsCount = 0;
                                done = false;

                                // there should be enough spare buffers to get a resultConvolutionBuffer for this operation now
                            } else {
                                // only one partially setup operation so there would be none to free up
                                // in this case we will use the reserve buffer
                                resultConvolutionBuffers[operationsCount] = getReserveBuffer();
                                convolveAndRelease(beagle, firstConvolutionBuffers, secondConvolutionBuffers, resultConvolutionBuffers, 1);
                                convolve.push(getReserveBuffer());
                                done = true; // break out of the do loop
                            }
                        }
                    } while (!done);

                    if (buffer >= 0) {
                        // if the buffer is still negative then the loop above will have used the reserve buffer
                        // to complete the convolution.
                        resultConvolutionBuffers[operationsCount] = buffer;
                        convolve.push(buffer);
                        operationsCount++;
                    }

                } else if (convolve.size() == 3) {
                    firstConvolutionBuffers[operationsCount] = convolve.pop();
                    secondConvolutionBuffers[operationsCount] = convolve.pop();
                    resultConvolutionBuffers[operationsCount] = convolve.pop();
                    operationsCount++;
                } else {
                    throw new RuntimeException("Unexpected convolve list size");
                }

                if (convolve.size() == 0) {
                    empty.add(convolve);
                }
            }

            if (DEBUG) {
                System.out.print("Convolving " + operationsCount+ " matrices:");
                for (int i = 0; i < operationsCount; i++) {
                    System.out.print(" " + firstConvolutionBuffers[i] + "*" + secondConvolutionBuffers[i] + "->" + resultConvolutionBuffers[i]);
                }
                System.out.println();
            }

            convolveAndRelease(beagle, firstConvolutionBuffers, secondConvolutionBuffers, resultConvolutionBuffers, operationsCount);

            convolutionList.removeAll(empty);
        }

        if (MEASURE_RUN_TIME) {
            timer.stop();
            double timeInSeconds = timer.toSeconds();
            convolveTime += timeInSeconds;
        }

    }// END: convolveTransitionMatrices

    private void convolveAndRelease(Beagle beagle, int[] firstConvolutionBuffers, int[] secondConvolutionBuffers, int[] resultConvolutionBuffers, int operationsCount) {

        if (RUN_IN_SERIES) {
            if (operationsCount > 1) {
                throw new RuntimeException("Unable to convolve matrices in series");
            }
        }

        beagle.convolveTransitionMatrices(firstConvolutionBuffers, // A
                secondConvolutionBuffers, // B
                resultConvolutionBuffers, // C
                operationsCount // count
        );

        for (int i = 0; i < operationsCount; i++) {
            if (firstConvolutionBuffers[i] >= matrixBufferHelper.getBufferCount() && firstConvolutionBuffers[i] != reserveBufferIndex) {
                pushAvailableBuffer(firstConvolutionBuffers[i]);
            }
            if (secondConvolutionBuffers[i] >= matrixBufferHelper.getBufferCount() && secondConvolutionBuffers[i] != reserveBufferIndex) {
                pushAvailableBuffer(secondConvolutionBuffers[i]);
            }
        }

    }//END: convolveAndRelease

    private int getAvailableBufferCount() {
        if (RUN_IN_SERIES) {
            return 0;
        } else {
            return availableBuffers.size();
        }
    }

    private int popAvailableBuffer() {
        if (availableBuffers.isEmpty()) {
            return -1;
        }
        return availableBuffers.pop();
    }

    /**
     * the reserve buffer is one extra buffer used to free up some spare buffers
     * @return
     */
    private int getReserveBuffer() {
        return reserveBufferIndex;
    }

    private void pushAvailableBuffer(int index) {
        availableBuffers.push(index);
    }

    @Override
    public void storeState() {
        eigenBufferHelper.storeState();
        matrixBufferHelper.storeState();
    }

    @Override
    public void restoreState() {
        eigenBufferHelper.restoreState();
        matrixBufferHelper.restoreState();
    }

}// END: class
