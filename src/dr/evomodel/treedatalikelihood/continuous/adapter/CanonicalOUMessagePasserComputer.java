/*
 * CanonicalOUMessagePasserComputer.java
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

package dr.evomodel.treedatalikelihood.continuous.adapter;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalRootPrior;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalTipObservation;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalTreeMessagePasser;
import dr.evomodel.treedatalikelihood.continuous.integration.SequentialCanonicalOUMessagePasser;
import dr.inference.model.Parameter;

/**
 * Wiring class for the canonical OU tree passer path.
 *
 * <p>This mirrors {@link GaussianMessagePasserComputer}, but targets the new
 * canonical OU APIs. This class is intentionally canonical-first: callers supply
 * exact tip observations through {@link CanonicalTipObservation} rather than
 * through the moment-form tree-partials layer.
 */
public final class CanonicalOUMessagePasserComputer {

    private final CanonicalTreeMessagePasser passer;
    private final CanonicalBranchTransitionProvider transitionProvider;
    private final CanonicalRootPrior rootPrior;
    private final double[][] scratchTraitCovariance;
    private final int rootIndex;

    public CanonicalOUMessagePasserComputer(final Tree tree,
                                            final MultivariateElasticModel elasticModel,
                                            final MultivariateDiffusionModel diffusionModel,
                                            final CanonicalTipObservation[] tipObservations,
                                            final ConjugateRootTraitPrior conjugatePrior,
                                            final Parameter stationaryMean,
                                            final BranchRateModel rateModel) {
        this(tree, elasticModel, diffusionModel, tipObservations, conjugatePrior, stationaryMean, rateModel, null);
    }

    public CanonicalOUMessagePasserComputer(final Tree tree,
                                            final MultivariateElasticModel elasticModel,
                                            final MultivariateDiffusionModel diffusionModel,
                                            final CanonicalTipObservation[] tipObservations,
                                            final ConjugateRootTraitPrior conjugatePrior,
                                            final Parameter stationaryMean,
                                            final BranchRateModel rateModel,
                                            final ContinuousRateTransformation rateTransformation) {
        final int dim = diffusionModel.getDimension();
        this.transitionProvider = new HomogeneousCanonicalOUBranchTransitionProvider(
                tree, elasticModel, diffusionModel, stationaryMean, rateModel, rateTransformation);
        this.rootPrior = new CanonicalConjugateRootPriorAdapter(conjugatePrior, dim);
        this.passer = new SequentialCanonicalOUMessagePasser(tree, dim);
        this.scratchTraitCovariance = new double[dim][dim];
        this.rootIndex = tree.getRoot().getNumber();
        loadTips(tipObservations);
    }

    public CanonicalOUMessagePasserComputer(final CanonicalTreeMessagePasser passer,
                                            final CanonicalBranchTransitionProvider transitionProvider,
                                            final CanonicalRootPrior rootPrior) {
        this.passer = passer;
        this.transitionProvider = transitionProvider;
        this.rootPrior = rootPrior;
        this.scratchTraitCovariance = new double[passer.getDimension()][passer.getDimension()];
        this.rootIndex = -1;
    }

    public double computeLogLikelihood() {
        return passer.computePostOrderLogLikelihood(transitionProvider, rootPrior);
    }

    public void computeGradientA(final double[] gradA) {
        passer.computePostOrderLogLikelihood(transitionProvider, rootPrior);
        passer.computePreOrder(transitionProvider, rootPrior);
        passer.computeGradientA(transitionProvider, gradA);
    }

    public void computeGradientQ(final double[] gradQ) {
        passer.computePostOrderLogLikelihood(transitionProvider, rootPrior);
        passer.computePreOrder(transitionProvider, rootPrior);
        passer.computeGradientQ(transitionProvider, gradQ);
    }

    public void computeGradientMu(final double[] gradMu) {
        passer.computePostOrderLogLikelihood(transitionProvider, rootPrior);
        passer.computePreOrder(transitionProvider, rootPrior);
        passer.computeGradientMu(transitionProvider, gradMu);
    }

    public void computeGradientRootMean(final double[] gradRootMean) {
        for (int i = 0; i < gradRootMean.length; i++) {
            gradRootMean[i] = 0.0;
        }
        if (rootIndex < 0) {
            throw new IllegalStateException("Root index is unavailable for this canonical OU message passer computer.");
        }
        passer.computePostOrderLogLikelihood(transitionProvider, rootPrior);
        transitionProvider.fillTraitCovariance(scratchTraitCovariance);
        rootPrior.accumulateRootMeanGradient(
                passer.getPostOrderState(rootIndex),
                scratchTraitCovariance,
                gradRootMean);
    }

    public void computeGradientBranchLengths(final double[] gradT) {
        passer.computePostOrderLogLikelihood(transitionProvider, rootPrior);
        passer.computePreOrder(transitionProvider, rootPrior);
        passer.computeGradientBranchLengths(transitionProvider, gradT);
    }

    public void reloadTips(final CanonicalTipObservation[] tipObservations) {
        loadTips(tipObservations);
    }

    public CanonicalTreeMessagePasser getPasser() {
        return passer;
    }

    public void storeState() {
        passer.storeState();
        transitionProvider.storeState();
    }

    public void restoreState() {
        passer.restoreState();
        transitionProvider.restoreState();
    }

    public void acceptState() {
        passer.acceptState();
        transitionProvider.acceptState();
    }

    private void loadTips(final CanonicalTipObservation[] tipObservations) {
        if (tipObservations.length != passer.getTipCount()) {
            throw new IllegalArgumentException(
                    "Tip observation count mismatch: " + tipObservations.length + " vs " + passer.getTipCount());
        }
        for (int tipIdx = 0; tipIdx < tipObservations.length; tipIdx++) {
            if (tipObservations[tipIdx] == null) {
                throw new IllegalArgumentException("tipObservations[" + tipIdx + "] must not be null");
            }
            passer.setTipObservation(tipIdx, tipObservations[tipIdx]);
        }
    }
}
