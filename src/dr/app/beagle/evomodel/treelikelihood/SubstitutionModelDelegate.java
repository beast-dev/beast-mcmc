/*
 * SubstitutionModelDelegate.java
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

package dr.app.beagle.evomodel.treelikelihood;

import beagle.Beagle;
import dr.app.beagle.evomodel.branchmodel.BranchModel;
import dr.app.beagle.evomodel.substmodel.EigenDecomposition;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.evolution.tree.Tree;

import java.util.*;

/**
 * @author Andrew Rambaut
 * @author Filip Bielejec
 * @author Marc A. Suchard
 * @version $Id$
 */
public final class SubstitutionModelDelegate {

    private static final boolean DEBUG = false;

    private static final int BUFFER_POOL_SIZE = 100;

    private final Tree tree;
    private final List<SubstitutionModel> substitutionModelList;
    private final BranchModel branchModel;

    private final int eigenCount;
    private final int nodeCount;

    private final int extraBufferCount;

    private final BufferIndexHelper eigenBufferHelper;
    private BufferIndexHelper matrixBufferHelper;

    private Deque<Integer> availableBuffers = new ArrayDeque<Integer>();

    public SubstitutionModelDelegate(Tree tree, BranchModel branchModel) {

        this.tree = tree;

        this.substitutionModelList = branchModel.getSubstitutionModels();

        this.branchModel = branchModel;

        eigenCount = substitutionModelList.size();
        nodeCount = tree.getNodeCount();

        // two eigen buffers for each decomposition for store and restore.
        eigenBufferHelper = new BufferIndexHelper(eigenCount, 0);

        // two matrices for each node less the root
        matrixBufferHelper = new BufferIndexHelper(nodeCount, 0);

        this.extraBufferCount = branchModel.requiresMatrixConvolution() ? BUFFER_POOL_SIZE : 0;

        for (int i = 0; i < extraBufferCount; i++) {
            pushAvailableBuffer(i + matrixBufferHelper.getBufferCount());
        }

    }// END: Constructor

    public boolean canReturnComplexDiagonalization() {
        return substitutionModelList.get(0).getEigenDecomposition().canReturnComplexDiagonalization();
    }

    public int getEigenBufferCount() {
        return eigenBufferHelper.getBufferCount();
    }

    public int getMatrixBufferCount() {
        return matrixBufferHelper.getBufferCount() + extraBufferCount;
    }

    public int getSubstitutionModelCount() {
        return substitutionModelList.size();
    }

    public SubstitutionModel getSubstitutionModel(int index) {
        return substitutionModelList.get(index);
    }

    public void updateSubstitutionModels(Beagle beagle) {
        for (int i = 0; i < eigenCount; i++) {
            eigenBufferHelper.flipOffset(i);

            EigenDecomposition ed = substitutionModelList.get(i).getEigenDecomposition();

            beagle.setEigenDecomposition(
                    eigenBufferHelper.getOffsetIndex(i),
                    ed.getEigenVectors(),
                    ed.getInverseEigenVectors(),
                    ed.getEigenValues());
        }
    }

    public void updateTransitionMatrices(Beagle beagle, int[] branchIndices, double[] edgeLength, int updateCount) {

        int[][] probabilityIndices = new int[eigenCount][updateCount];
        double[][] edgeLengths = new double[eigenCount][updateCount];

        int[] counts = new int[eigenCount];

        List<Deque<Integer>> convolutionList = new ArrayList<Deque<Integer>>();

        for (int i = 0; i < updateCount; i++) {

            BranchModel.Mapping mapping = branchModel.getBranchModelMapping(tree.getNode(branchIndices[i]));
            int[] order = mapping.getOrder();
            double[] weights = mapping.getWeights();

            if (order.length == 1) {
                probabilityIndices[order[0]][counts[order[0]]] = matrixBufferHelper.getOffsetIndex(branchIndices[i]);
                edgeLengths[order[0]][counts[order[0]]] = edgeLength[i];
                counts[order[0]] ++;
            } else {
                double sum = 0.0;
                for (double w : weights) {
                    sum += w;
                }

                Deque<Integer> bufferIndices = new ArrayDeque<Integer>();
                for (int j = 0; j < order.length; j++) {

                    int buffer;
                    boolean done;

                    do {
                        done = true;

                        buffer = popAvailableBuffer();

                        if (buffer < 0) {
                            // no buffers available
                            if (DEBUG) {
                                System.out.println("Ran out of buffers for transition matrices - computing current list.");
                            }
                            // we have run out of buffers, process what we have and continue...
                            computeTransitionMatrices(beagle, probabilityIndices, edgeLengths, counts);
                            convolveMatrices(beagle, convolutionList);

                            // reset the counts
                            for (int k = 0; k < eigenCount; k ++) {
                                counts[k] = 0;
                            }

                            done = false;
                        }
                    } while (!done);

                    probabilityIndices[order[j]][counts[order[j]]] = buffer;
                    edgeLengths[order[j]][counts[order[j]]] = weights[j] * edgeLength[i] / sum;
                    counts[order[j]]++;

                    bufferIndices.add(buffer);
                }
                bufferIndices.add(matrixBufferHelper.getOffsetIndex(branchIndices[i]));

                convolutionList.add(bufferIndices);
            }

        }

        computeTransitionMatrices(beagle, probabilityIndices, edgeLengths, counts);
        convolveMatrices(beagle, convolutionList);
    }

    private void computeTransitionMatrices(Beagle beagle, int[][] probabilityIndices, double[][] edgeLengths, int[] counts) {
        if (DEBUG) {
            System.out.print("Computing matrices:");
        }

        for (int i = 0; i < eigenCount; i++) {
            if (DEBUG) {
                for (int j = 0; j < counts[i]; j++) {
                    System.out.print(" " + probabilityIndices[i][j]);
                }
            }

            beagle.updateTransitionMatrices(eigenBufferHelper.getOffsetIndex(i),
                    probabilityIndices[i],
                    null, // firstDerivativeIndices
                    null, // secondDerivativeIndices
                    edgeLengths[i],
                    counts[i]);
        }

        if (DEBUG) {
            System.out.println();
        }
    }

    private void convolveMatrices(Beagle beagle, List<Deque<Integer>> convolutionList) {
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
                                System.out.print("Convolving matrices:");
                                for (int i = 0; i < operationsCount; i++) {
                                    System.out.print(" " + firstConvolutionBuffers[i] + "*" + secondConvolutionBuffers[i] + "->" + resultConvolutionBuffers[i]);
                                }
                                System.out.println();
                            }

                            beagle.convolveTransitionMatrices(firstConvolutionBuffers, // A
                                    secondConvolutionBuffers, // B
                                    resultConvolutionBuffers, // C
                                    operationsCount // count
                            );

                            operationsCount = 0;
                            done = false;
                        }
                    } while(!done);

                    resultConvolutionBuffers[operationsCount] = buffer;
                    convolve.push(buffer);
                    operationsCount ++;

                } else  if (convolve.size() == 3) {
                    firstConvolutionBuffers[operationsCount] = convolve.pop();
                    secondConvolutionBuffers[operationsCount] = convolve.pop();
                    resultConvolutionBuffers[operationsCount] = convolve.pop();
                    operationsCount ++;
                } else {
                    throw new RuntimeException("Unexpected convolve list size");
                }

                if (convolve.size() == 0) {
                    empty.add(convolve);
                }
            }

            if (DEBUG) {
                System.out.print("Convolving matrices:");
                for (int i = 0; i < operationsCount; i++) {
                    System.out.print(" " + firstConvolutionBuffers[i] + "*" + secondConvolutionBuffers[i] + "->" + resultConvolutionBuffers[i]);
                }
                System.out.println();
            }

            beagle.convolveTransitionMatrices(firstConvolutionBuffers, // A
                    secondConvolutionBuffers, // B
                    resultConvolutionBuffers, // C
                    operationsCount // count
            );

            for (int i = 0; i < operationsCount; i++) {
                if (firstConvolutionBuffers[i] >= matrixBufferHelper.getBufferCount()) {
                    pushAvailableBuffer(firstConvolutionBuffers[i]);
                }
                if (secondConvolutionBuffers[i] >= matrixBufferHelper.getBufferCount()) {
                    pushAvailableBuffer(secondConvolutionBuffers[i]);
                }
            }

            convolutionList.removeAll(empty);
        }
    }

    private int popAvailableBuffer() {
        if (availableBuffers.isEmpty()) {
            return -1;
        }
        return availableBuffers.pop();
    }

    private void pushAvailableBuffer(int index) {
        availableBuffers.push(index);
    }

    public double[] getRootStateFrequencies() {
        return substitutionModelList.get(0).getFrequencyModel().getFrequencies();
    }// END: getStateFrequencies

    public void flipMatrixBuffer(int branchIndex) {
        matrixBufferHelper.flipOffset(branchIndex);
    }

    public int getMatrixIndex(int branchIndex) {
        return matrixBufferHelper.getOffsetIndex(branchIndex);
    }

    public void storeState() {
        eigenBufferHelper.storeState();
        matrixBufferHelper.storeState();
    }

    public void restoreState() {
        eigenBufferHelper.restoreState();
        matrixBufferHelper.restoreState();
    }

}// END: class
