/*
 * CanonicalBMTreeLikelihoodIntegrator.java
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

package dr.evomodel.treedatalikelihood.continuous.canonical.adapter;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalRootPrior;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTreeMessagePasser;
import dr.evomodel.treedatalikelihood.continuous.canonical.SequentialCanonicalOUMessagePasser;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipObservation;

/**
 * Likelihood-only canonical tree backend for homogeneous Brownian motion.
 */
public final class CanonicalBMTreeLikelihoodIntegrator {

    private final CanonicalTreeMessagePasser passer;
    private final CanonicalBranchTransitionProvider transitionProvider;
    private final CanonicalRootPrior rootPrior;
    private final Tree tipTree;
    private final ContinuousTraitPartialsProvider tipDataModel;
    private final CanonicalTipObservation[] ownedTipObservations;
    private boolean tipObservationsDirty;
    private boolean storedTipObservationsDirty;

    public CanonicalBMTreeLikelihoodIntegrator(final Tree tree,
                                               final MultivariateDiffusionModel diffusionModel,
                                               final ContinuousTraitPartialsProvider dataModel,
                                               final ConjugateRootTraitPrior conjugatePrior,
                                               final BranchRateModel rateModel,
                                               final ContinuousRateTransformation rateTransformation) {
        final int dim = diffusionModel.getDimension();
        this.transitionProvider = new HomogeneousCanonicalBMBranchTransitionProvider(
                tree,
                diffusionModel,
                rateModel,
                rateTransformation);
        this.rootPrior = new CanonicalConjugateRootPriorAdapter(conjugatePrior, dim);
        this.passer = SequentialCanonicalOUMessagePasser.likelihoodOnly(tree, dim);
        this.tipTree = tree;
        this.tipDataModel = dataModel;
        this.ownedTipObservations = allocateCanonicalTipObservations(tree.getExternalNodeCount(), dim);
        this.tipObservationsDirty = true;
        this.storedTipObservationsDirty = true;
        syncOwnedTipObservations();
    }

    public double calculateLogLikelihood() {
        syncOwnedTipObservations();
        return passer.computePostOrderLogLikelihood(transitionProvider, rootPrior);
    }

    public void markTipObservationsDirty() {
        tipObservationsDirty = true;
    }

    public void markModelDirty() {
        if (transitionProvider instanceof HomogeneousCanonicalBMBranchTransitionProvider) {
            ((HomogeneousCanonicalBMBranchTransitionProvider) transitionProvider).makeDirty();
        }
    }

    public void storeState() {
        passer.storeState();
        transitionProvider.storeState();
        storedTipObservationsDirty = tipObservationsDirty;
    }

    public void restoreState() {
        passer.restoreState();
        transitionProvider.restoreState();
        tipObservationsDirty = storedTipObservationsDirty;
    }

    public void acceptState() {
        passer.acceptState();
        transitionProvider.acceptState();
        storedTipObservationsDirty = tipObservationsDirty;
    }

    private void syncOwnedTipObservations() {
        if (!tipObservationsDirty) {
            return;
        }
        CanonicalTipObservationAdapter.fillTipObservations(
                tipTree,
                tipDataModel,
                passer.getDimension(),
                ownedTipObservations);
        for (int tipIdx = 0; tipIdx < ownedTipObservations.length; tipIdx++) {
            passer.setTipObservation(tipIdx, ownedTipObservations[tipIdx]);
        }
        tipObservationsDirty = false;
    }

    private static CanonicalTipObservation[] allocateCanonicalTipObservations(final int tipCount, final int dim) {
        final CanonicalTipObservation[] observations = new CanonicalTipObservation[tipCount];
        for (int tipIdx = 0; tipIdx < tipCount; tipIdx++) {
            observations[tipIdx] = new CanonicalTipObservation(dim);
        }
        return observations;
    }
}
