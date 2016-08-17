/*
 * DataLikelihoodDelegate.java
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

import dr.inference.model.Model;

import java.util.List;

/**
 * DataLikelihoodDelegate - interface for a plugin delegate for the data likelihood.
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public interface DataLikelihoodDelegate extends Model {

    double calculateLikelihood(List<BranchOperation> branchOperations, List<NodeOperation> nodeOperations, int rootNodeNumber) throws LikelihoodUnderflowException;

    void makeDirty();

    void storeState();

    void restoreState();

    final class BranchOperation {
        public BranchOperation(int branchNumber, double branchLength) {
            this.branchNumber = branchNumber;
            this.branchLength = branchLength;
        }

        public int getBranchNumber() {
            return branchNumber;
        }

        public double getBranchLength() {
            return branchLength;
        }

        private final int branchNumber;
        private final double branchLength;
    }

    final class NodeOperation {
        public NodeOperation(int nodeNumber, int leftChild, int rightChild) {
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

        private final int nodeNumber;
        private final int leftChild;
        private final int rightChild;
    }

    class LikelihoodUnderflowException extends Exception { }
}
