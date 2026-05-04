/*
 * CanonicalTreeStateStore.java
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

package dr.evomodel.treedatalikelihood.continuous.canonical.traversal;

import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipObservation;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipObservationModel;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.IdentityCanonicalTipObservationModel;

/**
 * Owns canonical tree message buffers and their MCMC snapshots.
 *
 * <p>The message passer computes with these buffers, but this class owns the
 * storage layout and state lifecycle so traversal and gradient code can stay
 * focused on algorithms rather than buffer bookkeeping.</p>
 */
public final class CanonicalTreeStateStore {

    public CanonicalTipObservation[] tipObservations;
    private CanonicalTipObservation[] storedTipObservations;
    public CanonicalTipObservationModel[] tipObservationModels;
    private CanonicalTipObservationModel[] storedTipObservationModels;

    public CanonicalGaussianState[] postOrder;
    private CanonicalGaussianState[] storedPostOrder;

    public CanonicalGaussianState[] branchAboveParent;
    private CanonicalGaussianState[] storedBranchAboveParent;

    public boolean hasPostOrderState;
    private boolean storedHasPostOrderState;

    public boolean hasPreOrderState;
    private boolean storedHasPreOrderState;

    public double lastRootDiffusionScale;
    private double storedLastRootDiffusionScale;

    public boolean hasFixedRootValue;
    private boolean storedHasFixedRootValue;

    public double[] fixedRootValue;
    private double[] storedFixedRootValue;

    private final int nodeCount;
    private final int dimension;

    public CanonicalTreeStateStore(final int nodeCount, final int dimension) {
        this.nodeCount = nodeCount;
        this.dimension = dimension;
        this.tipObservations = allocateTipObservations(nodeCount, dimension);
        this.storedTipObservations = allocateTipObservations(nodeCount, dimension);
        this.tipObservationModels = allocateTipObservationModels(nodeCount, dimension);
        this.storedTipObservationModels = allocateTipObservationModels(nodeCount, dimension);
        this.postOrder = allocateStates(nodeCount, dimension);
        this.storedPostOrder = allocateStates(nodeCount, dimension);
        this.branchAboveParent = allocateStates(nodeCount, dimension);
        this.storedBranchAboveParent = allocateStates(nodeCount, dimension);
        this.fixedRootValue = new double[dimension];
        this.storedFixedRootValue = new double[dimension];
        this.hasPostOrderState = false;
        this.storedHasPostOrderState = false;
        this.hasPreOrderState = false;
        this.storedHasPreOrderState = false;
        this.lastRootDiffusionScale = 0.0;
        this.storedLastRootDiffusionScale = 0.0;
        this.hasFixedRootValue = false;
        this.storedHasFixedRootValue = false;
    }

    public void storeState() {
        storedHasPostOrderState = hasPostOrderState;
        storedHasPreOrderState = hasPreOrderState;
        storedLastRootDiffusionScale = lastRootDiffusionScale;
        storedHasFixedRootValue = hasFixedRootValue;
        System.arraycopy(fixedRootValue, 0, storedFixedRootValue, 0, dimension);
        for (int i = 0; i < nodeCount; i++) {
            storedTipObservations[i].copyFrom(tipObservations[i]);
            storedTipObservationModels[i] = tipObservationModels[i].copy();
            copyState(postOrder[i], storedPostOrder[i]);
            copyState(branchAboveParent[i], storedBranchAboveParent[i]);
        }
    }

    public void restoreState() {
        final CanonicalGaussianState[] tmpPost = postOrder;
        postOrder = storedPostOrder;
        storedPostOrder = tmpPost;

        final CanonicalGaussianState[] tmpAbove = branchAboveParent;
        branchAboveParent = storedBranchAboveParent;
        storedBranchAboveParent = tmpAbove;

        final boolean tmpHasPost = hasPostOrderState;
        hasPostOrderState = storedHasPostOrderState;
        storedHasPostOrderState = tmpHasPost;

        final boolean tmpHasPre = hasPreOrderState;
        hasPreOrderState = storedHasPreOrderState;
        storedHasPreOrderState = tmpHasPre;

        final double tmpRootScale = lastRootDiffusionScale;
        lastRootDiffusionScale = storedLastRootDiffusionScale;
        storedLastRootDiffusionScale = tmpRootScale;

        final boolean tmpHasFixedRoot = hasFixedRootValue;
        hasFixedRootValue = storedHasFixedRootValue;
        storedHasFixedRootValue = tmpHasFixedRoot;

        final double[] tmpRootValue = fixedRootValue;
        fixedRootValue = storedFixedRootValue;
        storedFixedRootValue = tmpRootValue;

        final CanonicalTipObservation[] tmpTips = tipObservations;
        tipObservations = storedTipObservations;
        storedTipObservations = tmpTips;

        final CanonicalTipObservationModel[] tmpTipModels = tipObservationModels;
        tipObservationModels = storedTipObservationModels;
        storedTipObservationModels = tmpTipModels;
    }

    public void setTipObservation(final int tipIndex, final CanonicalTipObservation observation) {
        tipObservations[tipIndex].copyFrom(observation);
        tipObservationModels[tipIndex] =
                IdentityCanonicalTipObservationModel.fromObservation(tipObservations[tipIndex]);
    }

    public void setTipObservationModel(final int tipIndex, final CanonicalTipObservationModel observationModel) {
        if (observationModel == null) {
            throw new IllegalArgumentException("Tip observation model must not be null");
        }
        if (observationModel.getLatentDimension() != dimension) {
            throw new IllegalArgumentException(
                    "Tip observation model latent dimension mismatch: "
                            + observationModel.getLatentDimension() + " vs " + dimension);
        }
        tipObservationModels[tipIndex] = observationModel.copy();
        if (observationModel instanceof IdentityCanonicalTipObservationModel) {
            ((IdentityCanonicalTipObservationModel) observationModel).copyTo(tipObservations[tipIndex]);
        } else {
            tipObservations[tipIndex].setMissing();
        }
    }

    private static CanonicalGaussianState[] allocateStates(final int count, final int dimension) {
        final CanonicalGaussianState[] out = new CanonicalGaussianState[count];
        for (int i = 0; i < count; i++) {
            out[i] = new CanonicalGaussianState(dimension);
        }
        return out;
    }

    private static CanonicalTipObservation[] allocateTipObservations(final int count, final int dimension) {
        final CanonicalTipObservation[] out = new CanonicalTipObservation[count];
        for (int i = 0; i < count; i++) {
            out[i] = new CanonicalTipObservation(dimension);
        }
        return out;
    }

    private static CanonicalTipObservationModel[] allocateTipObservationModels(final int count,
                                                                              final int dimension) {
        final CanonicalTipObservationModel[] out = new CanonicalTipObservationModel[count];
        for (int i = 0; i < count; i++) {
            out[i] = IdentityCanonicalTipObservationModel.missing(dimension);
        }
        return out;
    }

    private static void copyState(final CanonicalGaussianState source, final CanonicalGaussianState target) {
        CanonicalGaussianMessageOps.copyState(source, target);
    }
}
