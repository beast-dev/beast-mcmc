/*
 * AbstractTreeOperator.java
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

package dr.evomodel.operators;

import dr.evomodel.tree.TreeModel;
import dr.evolution.tree.*;
import dr.inference.model.Bounds;
import dr.inference.operators.SimpleMCMCOperator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Andrew Rambaut
 */
public abstract class AbstractTreeOperator extends SimpleMCMCOperator {

    private long transitions = 0;

    /**
     * @return the number of transitions since last call to reset().
     */
    public long getTransitions() {
        return transitions;
    }

    /**
     * Set the number of transitions since last call to reset(). This is used
     * to restore the state of the operator
     *
     * @param transitions number of transition
     */
    public void setTransitions(long transitions) {
        this.transitions = transitions;
    }

    public double getTransistionProbability() {
        final long accepted = getAcceptCount();
        final long rejected = getRejectCount();
        final long transition = getTransitions();
        return (double) transition / (double) (accepted + rejected);
    }

    /* exchange sub-trees whose root are i and j */
    protected void exchangeNodes(TreeModel tree, NodeRef i, NodeRef j,
                                 NodeRef iP, NodeRef jP) {

        tree.beginTreeEdit();
        tree.removeChild(iP, i);
        tree.removeChild(jP, j);
        tree.addChild(jP, i);
        tree.addChild(iP, j);

        tree.endTreeEdit();
    }

    public void reset() {
        super.reset();
        transitions = 0;
    }

    /**
     * @param tree   the tree
     * @param parent the parent
     * @param child  the child that you want the sister of
     * @return the other child of the given parent.
     */
    protected NodeRef getOtherChild(Tree tree, NodeRef parent, NodeRef child) {
        if( tree.getChild(parent, 0) == child ) {
            return tree.getChild(parent, 1);
        } else {
            return tree.getChild(parent, 0);
        }
    }

    public void addTreeUpdateCounts(int operationCount, int matrixCount) {
        operationCounts.add(operationCount);
        matrixCounts.add(matrixCount);
    }


    public void writeOperationCounts() {
        int sumOpCount = 0;
        int minOpCount = Integer.MAX_VALUE;
        int maxOpCount = Integer.MIN_VALUE;
        int sumMatCount = 0;
        int minMatCount = Integer.MAX_VALUE;
        int maxMatCount = Integer.MIN_VALUE;
        for (int i = 1; i < operationCounts.size(); i++) {
            int opCount = operationCounts.get(i);
            int matCount = matrixCounts.get(i);
            sumOpCount += opCount;
            if (opCount > maxOpCount) maxOpCount = opCount;
            if (opCount < minOpCount) minOpCount = opCount;
            sumMatCount += matCount;
            if (matCount > maxMatCount) maxMatCount = matCount;
            if (matCount < minMatCount) minMatCount = matCount;
        }

        System.out.println("Operator: " + getOperatorName());
        System.out.println("Operation count: " + operationCounts.size());
        System.out.println("Op count: " + ((double)sumOpCount) / (operationCounts.size() - 1));
        System.out.println("Op min,max: " + minOpCount + ", " + maxOpCount);
        System.out.println("Mat count: " + ((double)sumMatCount) / (matrixCounts.size() - 1));
        System.out.println("Mat min,max: " + minMatCount + ", " + maxMatCount);
        System.out.println();
    }

    List<Integer> operationCounts = new ArrayList<>();
    List<Integer> matrixCounts = new ArrayList<>();
}
