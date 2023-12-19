/*
 * ProcessOnTreeDelegate.java
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

package dr.evomodel.coalescent.basta;

/**
 * ProcessOnCoalescentIntervalDelegate - interface for a plugin delegate for the likelihood based on coalescent intervals
 *
 * @author Marc A. Suchard
 * @author Guy Baele
 * @version $Id$
 */
public interface ProcessOnCoalescentIntervalDelegate {

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

    final class OtherOperation {
        OtherOperation(int nodeNumber, int leftChild, int rightChild) {
            this.nodeNumber = nodeNumber;
            this.leftChild = leftChild;
            this.rightChild = rightChild;
        }

        public int getNodeNumber() {
            return nodeNumber;
        }

        public int getLeftChild() {
            return leftChild;
        }

        public int getRightChild() {
            return rightChild;
        }

        public String toString() {
            return nodeNumber + "(" + leftChild + "," + rightChild + ")";
        }

        private final int nodeNumber;
        private final int leftChild;
        private final int rightChild;
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
