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

import dr.evomodel.treedatalikelihood.TreeTraversal;

import java.util.List;

/**
 * ProcessOnCoalescentIntervalDelegate - interface for a plugin delegate for the likelihood based on coalescent intervals
 *
 * @author Marc A. Suchard
 * @author Guy Baele
 * @version $Id$
 */
public interface ProcessOnCoalescentIntervalDelegate {

    final class BranchIntervalOperation {
        BranchIntervalOperation(int intervalNumber, double intervalLength, int intervalOrder) {
            this.intervalNumber = intervalNumber;
            this.intervalLength = intervalLength;
            this.intervalOrder = intervalOrder;
        }

        public int getIntervalNumber() {
            return intervalNumber;
        }

        public double getIntervalLength() {
            return intervalLength;
        }

        public String toString() {
            return intervalNumber + ":" + intervalLength;
        }

        private final int intervalNumber;
        private final double intervalLength;
        private final int intervalOrder;
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
}
