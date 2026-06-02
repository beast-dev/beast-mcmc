/*
 * ProcessOnCoalescentIntervalDelegate.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.coalescent.basta;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static beagle.basta.BeagleBasta.BASTA_OPERATION_SIZE;

/**
 * ProcessOnCoalescentIntervalDelegate - interface for a plugin delegate for the likelihood based on coalescent intervals
 *
 * @author Marc A. Suchard
 * @author Guy Baele
 */
public interface ProcessOnCoalescentIntervalDelegate {

    class BranchIntervalOperationList {

        private static final boolean CACHE_FRIENDLY = true;

        private final List<BranchIntervalOperation> operations;
        private final List<Integer> starts;

        private final int[] nodeMap;
        private final int unmapped;
        private int used;

        public BranchIntervalOperationList(int numUnmapped, int maxCapacity) {
            this.operations = new ArrayList<>();
            this.starts = new ArrayList<>();
            this.unmapped = numUnmapped;

            if (CACHE_FRIENDLY) {
                this.nodeMap = new int[maxCapacity];
                Arrays.fill(nodeMap, -1);
                this.used = unmapped;
            } else {
                this.nodeMap = null;
            }
        }

        public void addOperation(BranchIntervalOperation operation) {
            operations.add(operation);
        }

        public void addStart() {
            starts.add(operations.size());
        }

        public List<BranchIntervalOperation> getOperations() {
            return operations;
        }

        public int map(int index) {
            if (!CACHE_FRIENDLY || index < unmapped) {
                return index;
            } else {
                if (nodeMap[index] == -1) {
                    nodeMap[index] = used;
                    ++used;
                }
                return nodeMap[index];
            }
        }

        public List<Integer> getStarts() {
            return starts;
        }

        public void clear() {
            operations.clear();
            starts.clear();

            Arrays.fill(nodeMap, -1);
            used = unmapped;
        }

        public void vectorize(int[] operations,
                              int[] intervals,
                              double[] lengths) {
            int k = 0;
            for (BranchIntervalOperation op : this.operations) {

                operations[k    ] = map(op.outputBuffer);
                operations[k + 1] = map(op.inputBuffer1);
                operations[k + 2] = op.inputMatrix1;
                operations[k + 3] = map(op.inputBuffer2);
                operations[k + 4] = op.inputMatrix2;
                operations[k + 5] = map(op.accBuffer1);
                operations[k + 6] = map(op.accBuffer2);
                operations[k + 7] = op.intervalNumber;

                k += BASTA_OPERATION_SIZE;
            }

            int i = 0;
            for (int end = starts.size() - 1; i < end; ++i) {
                int start = starts.get(i);
                intervals[i] = start;
                lengths[i] = this.operations.get(start).intervalLength;
            }
            intervals[i] = starts.get(i);
        }
    }

    final class BranchIntervalOperation {
        BranchIntervalOperation(int outputBuffer,
                                int inputBuffer1,
                                int inputBuffer2,
                                int inputMatrix1,
                                int inputMatrix2,
                                int accBuffer1,
                                int accBuffer2,
                                double intervalLength,
                                int executionOrder,
                                int intervalNumber) {
            this.outputBuffer = outputBuffer;
            this.inputBuffer1 = inputBuffer1;
            this.inputBuffer2 = inputBuffer2;
            this.inputMatrix1 = inputMatrix1;
            this.inputMatrix2 = inputMatrix2;
            this.accBuffer1 = accBuffer1;
            this.accBuffer2 = accBuffer2;
            this.intervalLength = intervalLength;
            this.executionOrder = executionOrder;
            this.intervalNumber = intervalNumber;
        }

        public String toString() {
            return intervalNumber + ":" + outputBuffer + " <- " +
                    inputBuffer1 + " (" + inputMatrix1 + ") + " +
                    inputBuffer2 + " (" + inputMatrix2 +  ") (" + intervalLength + ") ["+
                    accBuffer1 + " + " + accBuffer2 + "] @ " + executionOrder;
        }

        public final int outputBuffer;
        public final int inputBuffer1;
        public final int inputBuffer2;
        public final int inputMatrix1;
        public final int inputMatrix2;
        public final int accBuffer1;
        public final int accBuffer2;
        public final double intervalLength;
        public final int executionOrder;
        public final int intervalNumber;
    }

    final class TransitionMatrixOperation {
        TransitionMatrixOperation(int outputBuffer, int decompositionBuffer, double time) {
            this.outputBuffer = outputBuffer;
            this.decompositionBuffer = decompositionBuffer;
            this.time = time;
        }

        public String toString() { return outputBuffer + " <- " + decompositionBuffer + " (" + time + ")"; }

        public final int outputBuffer;
        public final int decompositionBuffer;
        public final double time;
    }
}
