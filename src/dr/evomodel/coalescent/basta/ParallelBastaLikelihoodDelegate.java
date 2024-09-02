/*
 * ParallelBastaLikelihoodDelegate.java
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

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

/**
 * @author Marc A. Suchard
 */
public class ParallelBastaLikelihoodDelegate extends GenericBastaLikelihoodDelegate {

    private static final int MIN_BRANCH_TASKS = 50;
    private static final int MIN_MATRIX_TASKS = 1000000;

    private final int threadCount;
    private final ExecutorService pool;
    private final double[][] temp;

    public ParallelBastaLikelihoodDelegate(String name,
                                           Tree tree,
                                           int stateCount,
                                           int threadCount,
                                           boolean transpose) {
        super(name, tree, stateCount, transpose);

        if (threadCount > 1) {
            pool = Executors.newFixedThreadPool(threadCount);
        } else if (threadCount < -1) {
            pool = Executors.newCachedThreadPool();
        } else {
            throw new IllegalArgumentException("Illegal threadCount value");
        }

        this.threadCount = Math.abs(threadCount);
        this.temp = new double[this.threadCount][stateCount * stateCount];
    }

    @Override
    protected void computeInnerBranchIntervalOperations(List<BranchIntervalOperation> branchIntervalOperations,
                                                      int start, int end) {

        int totalTasks = end - start;
//        System.err.println(totalTasks);

        if (totalTasks <= MIN_BRANCH_TASKS) {
            super.computeInnerBranchIntervalOperations(branchIntervalOperations, start, end);
        } else {
            forkJoin((s, e, t) -> super.computeInnerBranchIntervalOperations(branchIntervalOperations, s, e),
                    start, end);
        }
    }

    @Override
    protected void computeInnerTransitionProbabilityOperationsGrad(List<TransitionMatrixOperation> matrixOperations, int start, int end) {
        int totalTasks = end - start;

        if (totalTasks <= MIN_MATRIX_TASKS) {
            super.computeInnerTransitionProbabilityOperationsGrad(matrixOperations, start, end);
        } else {
            forkJoin((s, e, t) -> super.computeInnerTransitionProbabilityOperationsGrad(matrixOperations, s, e),
                    start, end);
        }
    }

    @Override
    protected void computeInnerBranchIntervalOperationsGrad(List<BranchIntervalOperation> branchIntervalOperations, List<TransitionMatrixOperation> matrixOperations, int start, int end){

        int totalTasks = end - start;
//        System.err.println(totalTasks);

        if (totalTasks <= MIN_BRANCH_TASKS) {
            super.computeInnerBranchIntervalOperationsGrad(branchIntervalOperations, matrixOperations, start, end);
        } else {
            forkJoin((s, e, t) -> super.computeInnerBranchIntervalOperationsGrad(branchIntervalOperations, matrixOperations, s, e),
                    start, end);
        }
    }

    @Override
    protected void computeInnerTransitionProbabilityOperations(List<TransitionMatrixOperation> matrixOperations,
                                                               int start, int end, double[] temp) {
        int totalTasks = end - start;

        if (totalTasks <= MIN_MATRIX_TASKS) {
            super.computeInnerTransitionProbabilityOperations(matrixOperations, start, end, temp);
        } else {
            forkJoin((s, e, t) -> super.computeInnerTransitionProbabilityOperations(matrixOperations, s, e, this.temp[t]),
                    start, end);
        }
    }

    private interface RangeCallable {
        void execute(int start, int end, int task);
    }

    private void forkJoin(RangeCallable f, int start, int end) {
        final int totalTasks = end - start;

        int numTasksPerThread = totalTasks / threadCount;
        if (totalTasks % threadCount != 0) {
            ++numTasksPerThread;
        }

        List<Callable<Object>> tasks = new ArrayList<>(threadCount);
        for (int i = 0, startTask = start; startTask < end; ++i, startTask += numTasksPerThread) {
            final int thisStart = startTask;
            final int thisEnd = Math.min(end, startTask + numTasksPerThread);
            final int thisTask = i;

            tasks.add(() -> {
                f.execute(thisStart, thisEnd, thisTask);
                return null;
            });

        }
        try {
            pool.invokeAll(tasks);
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        }
    }
}
