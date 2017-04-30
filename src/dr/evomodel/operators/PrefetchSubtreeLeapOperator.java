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
        return currentPrefetch == prefetchCount;
    }

    @Override
    public double doOperation() {

        if (currentPrefetch < 0) {
            // Run the operator N times, cache the treeDataLikelihood evaluations.
            // store the hastings ratios. TreeDataLikelihood will also store likelihoods.
            currentPrefetch = 0;

            for (int i = 0; i < prefetchCount; i++) {
                treeDataLikelihood.setCurrentPrefetch(i);

                hastingsRatios[i] = super.doOperation();

                treeDataLikelihood.collectOperations();

                getTreeModel().restoreModelState();

                instances[i] = getLastInstance();
            }

            treeDataLikelihood.prefetchLogLikelihoods();
        }

        treeDataLikelihood.setCurrentPrefetch(currentPrefetch);

        return hastingsRatios[currentPrefetch];
    }

    @Override
    public void accept(double deviation) {
        super.accept(deviation);
        // clear the cached likelihoods, restore treemodel, TDL and apply successful operation instance
        getTreeModel().restoreModelState();
        treeDataLikelihood.restoreModelState();

        applyInstance(instances[currentPrefetch]);
        
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

    }

    public String getOperatorName() {
        return "Prefetch" + super.getOperatorName();
    }

}
