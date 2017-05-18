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

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.PrefetchTreeDataLikelihood;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.Prefetchable;

/**
 * PrefetchSubtreeLeapOperator
 *
 * @author Andrew Rambaut
 */

public class PrefetchSubtreeLeapOperator extends SubtreeLeapOperator implements Prefetchable {

    private static final boolean NO_PARALLEL_PREFETCH = false;
    private static final boolean PREFETCH_DEBUG = false;

    private final PrefetchTreeDataLikelihood treeDataLikelihood;
    private final int prefetchCount;
    private int currentPrefetch;
    private final double[] hastingsRatios;
    private final Instance[] instances;

    /**
     * Constructor
     *
     * @param tree   the tree
     * @param weight the weight
     * @param size   scaling on a unit Gaussian to draw the patristic distance from
     * @param accP
     * @param mode   coercion mode
     */
    public PrefetchSubtreeLeapOperator(PrefetchTreeDataLikelihood treeDataLikelihood,
                                       TreeModel tree, double weight, double size, double accP, CoercionMode mode) {
        super(tree, weight, size, accP, mode);

        this.treeDataLikelihood = treeDataLikelihood;
        this.prefetchCount = treeDataLikelihood.getPrefetchCount();
        hastingsRatios = new double[prefetchCount];
        instances = new Instance[prefetchCount];

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
                System.err.println("Prefetch: Drawing " + prefetchCount + " subtree leaps");
            }

            if (NO_PARALLEL_PREFETCH) {
                // A debugging option where there is no parallel processing of the operator
                // moves but they are simply done in sequence. The N operator instances
                // are drawn so that the random number sequence is conserved (allowing for
                // comparison with the parallel approach) and these are then just applied
                // in sequence as doOperation is called.

                for (int i = 0; i < prefetchCount; i++) {
                    instances[i] = drawOperation();
                    hastingsRatios[i] =  instances[i].logHastingsRatio;
                }
            } else {

                for (int i = 0; i < prefetchCount; i++) {
                    treeDataLikelihood.setCurrentPrefetch(i);

                    // temporarily store state of tree
                    getTreeModel().pushState();

                    hastingsRatios[i] = super.doOperation();

                    treeDataLikelihood.collectOperations();

                    // restore temporary state
                    getTreeModel().popState();

                    instances[i] = getLastInstance();
                }

                treeDataLikelihood.prefetchLogLikelihoods();
            }
        }

        if (PREFETCH_DEBUG) {
            System.err.println("Prefetch: performing operation " + (currentPrefetch + 1));
        }

        if (NO_PARALLEL_PREFETCH) {
            applyInstance(instances[currentPrefetch]);
        } else {
            treeDataLikelihood.setCurrentPrefetch(currentPrefetch);
        }

        return hastingsRatios[currentPrefetch];
    }

    @Override
    public void accept(double deviation) {
        super.accept(deviation);

        // clear the cached likelihoods, restore treemodel, TDL and apply successful operation instance
        if (!NO_PARALLEL_PREFETCH) {

            // would be better to adopt the partial likelihoods of the successful move.
            // getTreeModel().quietlyRestoreModelState();
            // treeDataLikelihood.quietlyRestoreModelState();

            // todo - rather than rerunning the winning operation and then needing
            // todo   update all the touched partial likelihoods, it would be better
            // todo   shuffle the buffers in the deligate
            treeDataLikelihood.setCurrentPrefetch(-1);
            applyInstance(instances[currentPrefetch]);

            // update the log likelihood (partials etc).
            treeDataLikelihood.getLogLikelihood();
        }

        if (PREFETCH_DEBUG) {
            System.err.println("Prefetch: accepted");
        }

        currentPrefetch = -1;
    }

    @Override
    public void reject() {
        super.reject();

        // move on to the next prefetched
        currentPrefetch += 1;
        if (currentPrefetch == prefetchCount) {
            currentPrefetch = -1;
        }

        if (PREFETCH_DEBUG) {
            System.err.println("Prefetch: rejected");
        }

    }

    public String getOperatorName() {
        return "Prefetch" + super.getOperatorName();
    }

}
