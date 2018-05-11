/*
 * PrefetchSubtreeLeapOperator.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.operators;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.PrefetchableLikelihood;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.PrefetchableOperator;

import java.util.logging.Logger;

/**
 * PrefetchSubtreeLeapOperator
 *
 * @author Andrew Rambaut
 */

public class PrefetchSubtreeLeapOperator extends SubtreeLeapOperator implements PrefetchableOperator {

    // for debugging purposes switch this on to calculate the moves in series.
    private static final boolean NO_PARALLEL_PREFETCH = false;

    // provide debugging strings to stderr
    private static final boolean PREFETCH_DEBUG = false;

    private final PrefetchableLikelihood prefetchableLikelihood;
    private final int prefetchCount;
    private int currentPrefetch;
    private final Instance[] instances;
    private final double[] logHastingsRatios;

    /**
     * Constructor
     *
     * @param tree   the tree
     * @param weight the weight
     * @param size   scaling on a unit Gaussian to draw the patristic distance from
     * @param accP
     * @param mode   coercion mode
     */
    public PrefetchSubtreeLeapOperator(PrefetchableLikelihood prefetchableLikelihood,
                                       TreeModel tree, double weight, double size, double accP, CoercionMode mode) {
        super(tree, weight, size, accP, mode);

        final Logger logger = Logger.getLogger("dr.evomodel");
        if (!NO_PARALLEL_PREFETCH) {
            logger.info("\nUsing Parallel-Prefetch SubtreeLeapOperator");
        } else {
            logger.info("\nUsing Non-Parallel-Prefetch SubtreeLeapOperator");
        }

        this.prefetchableLikelihood = prefetchableLikelihood;
        this.prefetchCount = prefetchableLikelihood.getPrefetchCount();
        instances = new Instance[prefetchCount];
        logHastingsRatios = new double[prefetchCount];

        currentPrefetch = -1;
    }

    @Override
    public boolean prefetchingDone() {
        return currentPrefetch == prefetchCount - 1;
    }

    @Override
    public double doOperation() {

        if (currentPrefetch < 0) {
            // Run the operator N times, cache the treeDataLikelihood evaluations.
            // store the hastings ratios. TreeDataLikelihood will also store likelihoods.
            currentPrefetch = 0;
            if (PREFETCH_DEBUG) {
                System.out.println("Prefetch: Drawing " + prefetchCount + " subtree leaps");
            }

            for (int i = 0; i < prefetchCount; i++) {
                instances[i] = drawOperation();
            }

            if (!NO_PARALLEL_PREFETCH) {

                getTreeModel().pushState();

                for (int i = 0; i < prefetchCount; i++) {

//                    System.out.println("PSTL, before push:     " +getTreeModel().getNewick());
                    // temporarily store state of tree

                    prefetchableLikelihood.startPrefetchOperation(i);
                    prefetchableLikelihood.setIgnoreTreeEvents(false);
                    logHastingsRatios[i] = applyInstance(instances[i]);
                    prefetchableLikelihood.finishPrefetchOperation(i);

//                    System.out.println("PSTL, after operation: " +getTreeModel().getNewick());

                    // restore temporary state
                    getTreeModel().popState();
//                    System.out.println("PSTL, after pop:       " +getTreeModel().getNewick());
//                    System.out.println();
                }

                prefetchableLikelihood.prefetchLogLikelihoods();
            }
        }

        if (PREFETCH_DEBUG) {
            System.out.println("Prefetch: performing operation " + (currentPrefetch + 1));
        }

        if (NO_PARALLEL_PREFETCH) {
            // A debugging option where there is no parallel processing of the operator
            // moves but they are simply done in sequence. The N operator instances
            // are drawn so that the random number sequence is conserved (allowing for
            // comparison with the parallel approach) and these are then just applied
            // in sequence as doOperation is called.
//            System.out.println("PSTL, before apply:     " +getTreeModel().getNewick());
            return applyInstance(instances[currentPrefetch]);
//            System.out.println("PSTL, after apply:     " +getTreeModel().getNewick());
        } else {
            prefetchableLikelihood.setPrefetchLikelihood(currentPrefetch);

            // Need to get the tree back to the appropriate one so the tree priors calculate
            // on the correct one.
            prefetchableLikelihood.setIgnoreTreeEvents(true);
            applyInstance(instances[currentPrefetch]);
        }

        return logHastingsRatios[currentPrefetch];
    }

    @Override
    public void accept(double deviation) {
        super.accept(deviation);

        // clear the cached likelihoods, restore treemodel, TDL and apply successful operation instance
        if (!NO_PARALLEL_PREFETCH) {
            prefetchableLikelihood.acceptPrefetch(currentPrefetch);

            // This is not necessary any more - the operation would have been re-applied in the last call
            // of doOperation().
            // the original tree would have been popped so re-apply the accepted operation but
            // tell the treeDataLikelihood to ignore the tree changed events.
//            prefetchableLikelihood.setIgnoreTreeEvents(true);
//            applyInstance(instances[currentPrefetch]);

            //            System.out.println("PSTL, after accept:    " + getTreeModel().getNewick());
        } else {
            prefetchableLikelihood.releaseBufferIndices(0);
        }

        if (PREFETCH_DEBUG) {
            System.out.println("Prefetch: accepted");
        }

        // done with this prefetch set of operations
        currentPrefetch = -1;
    }

    @Override
    public void reject() {
        super.reject();

        // move on to the next prefetched
        currentPrefetch += 1;
        if (currentPrefetch == prefetchCount) {
            // depleted all the prefetched likelihoods so flag for a full
            // restoreState
            if (!NO_PARALLEL_PREFETCH) {
                prefetchableLikelihood.rejectAllPrefetches();
            }

            if (PREFETCH_DEBUG) {
                System.out.println("Prefetch: rejected all prefetches");
            }

            currentPrefetch = -1;
        } else if (PREFETCH_DEBUG) {
            System.out.println("Prefetch: rejected");
        }

    }

    public String getOperatorName() {
        return "Prefetch" + super.getOperatorName();
    }

}
