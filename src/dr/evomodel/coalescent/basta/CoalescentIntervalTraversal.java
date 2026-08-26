/*
 * CoalescentIntervalTraversal.java
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

import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.coalescent.StructuredCoalescentSchedule;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.TreeTraversal;

import java.util.ArrayList;
import java.util.List;

import static dr.evomodel.coalescent.basta.ProcessOnCoalescentIntervalDelegate.*;

/**
 * @author Marc A Suchard
 * @author Guy Baele
 */
public class CoalescentIntervalTraversal extends TreeTraversal {

    private final BigFastTreeIntervals treeIntervals;
    private final boolean checkForZeroLengthIntervals;
    private final BastaBranchIntervalOperationBuilder operationBuilder;

    private int currentLikelihoodInterval;

    protected CoalescentIntervalTraversal(final Tree tree,
                                          final BigFastTreeIntervals treeIntervals,
                                          final BranchRateModel branchRateModel,
                                          final int numberSubIntervals) {
        this(tree, treeIntervals, branchRateModel, numberSubIntervals, true);
    }

    protected CoalescentIntervalTraversal(final Tree tree,
                                          final BigFastTreeIntervals treeIntervals,
                                          final BranchRateModel branchRateModel,
                                          final int numberSubIntervals,
                                          final boolean checkForZeroLengthIntervals) {
        super(tree, branchRateModel, TraversalType.REVERSE_LEVEL_ORDER);

        assert tree instanceof TreeModel;

        this.treeIntervals = treeIntervals;
        this.checkForZeroLengthIntervals = checkForZeroLengthIntervals;
        this.operationBuilder = new BastaBranchIntervalOperationBuilder(
                tree, branchRateModel, numberSubIntervals, checkForZeroLengthIntervals,
                branchIntervalOperations, matrixOperations, intervalStarts);
    }

    @Override
    public final void dispatchTreeTraversalCollectBranchAndNodeOperations() {
        matrixOperations.clear();
        branchIntervalOperations.clear();
        intervalStarts.clear();

        if (traversalType == TraversalType.REVERSE_LEVEL_ORDER) {
            traverseReverseCoalescentLevelOrder();
        } else {
            assert false : "Unknown traversal type";
        }
    }

    public List<BranchIntervalOperation> getBranchIntervalOperations() {
        return branchIntervalOperations;
    }

    public List<TransitionMatrixOperation> getMatrixOperations() {
        return matrixOperations;
    }

    public int getCoalescentIntervalCount() {
        return currentLikelihoodInterval + 1;
    }

    public List<Integer> getIntervalStarts() {
        return intervalStarts;
    }

    private void traverseReverseCoalescentLevelOrder() {
        final StructuredCoalescentSchedule schedule =
                StructuredCoalescentSchedule.fromTreeIntervals(treeModel, treeIntervals, true,
                        checkForZeroLengthIntervals);
        operationBuilder.build(schedule);
        currentLikelihoodInterval = operationBuilder.getCoalescentIntervalCount() - 1;
    }

    private final List<BranchIntervalOperation> branchIntervalOperations = new ArrayList<>();
    private final List<TransitionMatrixOperation> matrixOperations = new ArrayList<>();
    private final List<Integer> intervalStarts = new ArrayList<>();

}
