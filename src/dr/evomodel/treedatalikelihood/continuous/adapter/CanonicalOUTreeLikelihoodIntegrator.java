/*
 * CanonicalOUTreeLikelihoodIntegrator.java
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
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalOUTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalRootPrior;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalTipObservation;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalTransitionCacheDiagnostics;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalTreeMessagePasser;
import dr.evomodel.treedatalikelihood.continuous.integration.SequentialCanonicalOUMessagePasser;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;
import dr.inference.model.Parameter;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockCanonicalParameterization;

/**
 * Integrator-style backend for the canonical OU tree passer path.
 * This class owns the canonical message buffers, reusable
 * tip observations, gradient caches, and MCMC state lifecycle for the tree path.
 */
public final class CanonicalOUTreeLikelihoodIntegrator implements CanonicalOUIntegrator {

    private static final int CACHE_DIAGNOSTIC_REPORT_INTERVAL =
            Integer.getInteger("beast.debug.canonicalTransitionCacheReportEvery", 100);

    private final CanonicalTreeMessagePasser passer;
    private final CanonicalBranchTransitionProvider transitionProvider;
    private final CanonicalRootPrior rootPrior;
    private final double[][] scratchTraitCovariance;
    private final int rootIndex;
    private final double[] cachedGradientA;
    private final double[] cachedGradientQ;
    private final double[] cachedGradientMu;
    private final double[] storedGradientA;
    private final double[] storedGradientQ;
    private final double[] storedGradientMu;
    private final Tree tipTree;
    private final ContinuousTraitPartialsProvider tipDataModel;
    private final CanonicalTipObservation[] ownedTipObservations;
    private boolean jointGradientCacheDirty;
    private boolean storedJointGradientCacheDirty;
    private boolean modelDirty;
    private boolean storedModelDirty;
    private boolean tipObservationsDirty;
    private boolean storedTipObservationsDirty;
    private int cacheDiagnosticReportCounter;

    /*
     * Standalone constructor for callers that already have canonical tip
     * observations and want explicit control over reloading them.
     */
    public CanonicalOUTreeLikelihoodIntegrator(final Tree tree,
                                               final MultivariateElasticModel elasticModel,
                                               final MultivariateDiffusionModel diffusionModel,
                                               final CanonicalTipObservation[] tipObservations,
                                               final ConjugateRootTraitPrior conjugatePrior,
                                               final Parameter stationaryMean,
                                               final BranchRateModel rateModel) {
        this(tree, elasticModel, diffusionModel, tipObservations, conjugatePrior, stationaryMean, rateModel, null);
    }

    /*
     * Delegate-owned constructor. The integrator owns reusable canonical tip
     * observations and refreshes them lazily from the data model.
     */
    public CanonicalOUTreeLikelihoodIntegrator(final Tree tree,
                                               final MultivariateElasticModel elasticModel,
                                               final MultivariateDiffusionModel diffusionModel,
                                               final ContinuousTraitPartialsProvider dataModel,
                                               final ConjugateRootTraitPrior conjugatePrior,
                                               final Parameter stationaryMean,
                                               final BranchRateModel rateModel,
                                               final ContinuousRateTransformation rateTransformation) {
        this(tree,
                elasticModel,
                diffusionModel,
                allocateCanonicalTipObservations(tree.getExternalNodeCount(), diffusionModel.getDimension()),
                conjugatePrior,
                stationaryMean,
                rateModel,
                rateTransformation,
                dataModel);
    }

    /*
     * Standalone constructor with an explicit rate transformation.
     */
    public CanonicalOUTreeLikelihoodIntegrator(final Tree tree,
                                               final MultivariateElasticModel elasticModel,
                                               final MultivariateDiffusionModel diffusionModel,
                                               final CanonicalTipObservation[] tipObservations,
                                               final ConjugateRootTraitPrior conjugatePrior,
                                               final Parameter stationaryMean,
                                               final BranchRateModel rateModel,
                                               final ContinuousRateTransformation rateTransformation) {
        this(tree, elasticModel, diffusionModel, tipObservations, conjugatePrior, stationaryMean, rateModel,
                rateTransformation, null);
    }

    private CanonicalOUTreeLikelihoodIntegrator(final Tree tree,
                                                final MultivariateElasticModel elasticModel,
                                                final MultivariateDiffusionModel diffusionModel,
                                                final CanonicalTipObservation[] tipObservations,
                                                final ConjugateRootTraitPrior conjugatePrior,
                                                final Parameter stationaryMean,
                                                final BranchRateModel rateModel,
                                                final ContinuousRateTransformation rateTransformation,
                                                final ContinuousTraitPartialsProvider dataModel) {
        final int dim = diffusionModel.getDimension();
        this.transitionProvider = new HomogeneousCanonicalOUBranchTransitionProvider(
                tree, elasticModel, diffusionModel, stationaryMean, rateModel, rateTransformation);
        this.rootPrior = new CanonicalConjugateRootPriorAdapter(conjugatePrior, dim);
        this.passer = new SequentialCanonicalOUMessagePasser(tree, dim);
        this.scratchTraitCovariance = new double[dim][dim];
        this.rootIndex = tree.getRoot().getNumber();
        final int selectionGradientDimension = inferSelectionGradientDimension(this.transitionProvider, this.passer);
        this.cachedGradientA = new double[selectionGradientDimension];
        this.cachedGradientQ = new double[dim * dim];
        this.cachedGradientMu = new double[dim];
        this.storedGradientA = new double[selectionGradientDimension];
        this.storedGradientQ = new double[dim * dim];
        this.storedGradientMu = new double[dim];
        this.tipTree = dataModel == null ? null : tree;
        this.tipDataModel = dataModel;
        this.ownedTipObservations = dataModel == null ? null : tipObservations;
        this.jointGradientCacheDirty = true;
        this.storedJointGradientCacheDirty = true;
        this.modelDirty = true;
        this.storedModelDirty = true;
        this.tipObservationsDirty = dataModel != null;
        this.storedTipObservationsDirty = this.tipObservationsDirty;
        this.cacheDiagnosticReportCounter = 0;
        if (dataModel == null) {
            loadTips(tipObservations);
        } else {
            syncOwnedTipObservations();
        }
    }

    /*
     * Unit-test constructor for injecting mock/counting passers and providers.
     */
    public CanonicalOUTreeLikelihoodIntegrator(final CanonicalTreeMessagePasser passer,
                                               final CanonicalBranchTransitionProvider transitionProvider,
                                               final CanonicalRootPrior rootPrior) {
        this.passer = passer;
        this.transitionProvider = transitionProvider;
        this.rootPrior = rootPrior;
        this.scratchTraitCovariance = new double[passer.getDimension()][passer.getDimension()];
        this.rootIndex = -1;
        final int selectionGradientDimension = inferSelectionGradientDimension(this.transitionProvider, this.passer);
        this.cachedGradientA = new double[selectionGradientDimension];
        this.cachedGradientQ = new double[passer.getDimension() * passer.getDimension()];
        this.cachedGradientMu = new double[passer.getDimension()];
        this.storedGradientA = new double[selectionGradientDimension];
        this.storedGradientQ = new double[passer.getDimension() * passer.getDimension()];
        this.storedGradientMu = new double[passer.getDimension()];
        this.tipTree = null;
        this.tipDataModel = null;
        this.ownedTipObservations = null;
        this.jointGradientCacheDirty = true;
        this.storedJointGradientCacheDirty = true;
        this.modelDirty = true;
        this.storedModelDirty = true;
        this.tipObservationsDirty = false;
        this.storedTipObservationsDirty = false;
        this.cacheDiagnosticReportCounter = 0;
    }

    @Override
    public double calculateLogLikelihood() {
        syncOwnedTipObservations();
        final double logLikelihood = passer.computePostOrderLogLikelihood(transitionProvider, rootPrior);
        modelDirty = false;
        reportTransitionCacheDiagnostics("logLikelihood");
        return logLikelihood;
    }

    @Override
    public void computeSelectionGradient(final double[] gradA) {
        syncOwnedTipObservations();
        ensureJointGradientCache();
        copyGradient(cachedGradientA, gradA, "selection");
    }

    @Override
    public void computeDiffusionGradient(final double[] gradQ) {
        syncOwnedTipObservations();
        ensureJointGradientCache();
        copyGradient(cachedGradientQ, gradQ, "diffusion");
    }

    @Override
    public void computeStationaryMeanGradient(final double[] gradMu) {
        syncOwnedTipObservations();
        ensureJointGradientCache();
        copyGradient(cachedGradientMu, gradMu, "stationary mean");
    }

    @Override
    public void computeGradientRootMean(final double[] gradRootMean) {
        syncOwnedTipObservations();
        for (int i = 0; i < gradRootMean.length; i++) {
            gradRootMean[i] = 0.0;
        }
        if (rootIndex < 0) {
            throw new IllegalStateException("Root index is unavailable for this canonical OU tree likelihood integrator.");
        }
        if (modelDirty || jointGradientCacheDirty) {
            passer.computePostOrderLogLikelihood(transitionProvider, rootPrior);
            modelDirty = false;
        }
        transitionProvider.fillTraitCovariance(scratchTraitCovariance);
        rootPrior.accumulateRootMeanGradient(
                passer.getPostOrderState(rootIndex),
                scratchTraitCovariance,
                gradRootMean);
    }

    @Override
    public void computeGradientBranchLengths(final double[] gradT) {
        syncOwnedTipObservations();
        passer.computePostOrderLogLikelihood(transitionProvider, rootPrior);
        passer.computePreOrder(transitionProvider, rootPrior);
        modelDirty = false;
        passer.computeGradientBranchLengths(transitionProvider, gradT);
        reportTransitionCacheDiagnostics("branchLengthGradient");
    }

    public void reloadTips(final CanonicalTipObservation[] tipObservations) {
        loadTips(tipObservations);
        markJointGradientCacheDirty();
    }

    @Override
    public void markTipObservationsDirty() {
        if (ownedTipObservations == null) {
            return;
        }
        tipObservationsDirty = true;
        markJointGradientCacheDirty();
    }

    @Override
    public void storeState() {
        passer.storeState();
        transitionProvider.storeState();
        storedJointGradientCacheDirty = jointGradientCacheDirty;
        storedModelDirty = modelDirty;
        storedTipObservationsDirty = tipObservationsDirty;
        if (!jointGradientCacheDirty) {
            System.arraycopy(cachedGradientA, 0, storedGradientA, 0, cachedGradientA.length);
            System.arraycopy(cachedGradientQ, 0, storedGradientQ, 0, cachedGradientQ.length);
            System.arraycopy(cachedGradientMu, 0, storedGradientMu, 0, cachedGradientMu.length);
        }
    }

    @Override
    public void restoreState() {
        passer.restoreState();
        transitionProvider.restoreState();
        jointGradientCacheDirty = storedJointGradientCacheDirty;
        modelDirty = storedModelDirty;
        tipObservationsDirty = storedTipObservationsDirty;
        if (!jointGradientCacheDirty) {
            System.arraycopy(storedGradientA, 0, cachedGradientA, 0, storedGradientA.length);
            System.arraycopy(storedGradientQ, 0, cachedGradientQ, 0, storedGradientQ.length);
            System.arraycopy(storedGradientMu, 0, cachedGradientMu, 0, storedGradientMu.length);
        }
    }

    @Override
    public void acceptState() {
        passer.acceptState();
        transitionProvider.acceptState();
        storedJointGradientCacheDirty = jointGradientCacheDirty;
        storedModelDirty = modelDirty;
        storedTipObservationsDirty = tipObservationsDirty;
        if (!jointGradientCacheDirty) {
            System.arraycopy(cachedGradientA, 0, storedGradientA, 0, cachedGradientA.length);
            System.arraycopy(cachedGradientQ, 0, storedGradientQ, 0, cachedGradientQ.length);
            System.arraycopy(cachedGradientMu, 0, storedGradientMu, 0, cachedGradientMu.length);
        }
    }

    @Override
    public void markModelDirty() {
        // Transition-provider cache invalidation is owned by the provider's
        // model listeners; this flag covers integrator-owned derived state.
        modelDirty = true;
        markJointGradientCacheDirty();
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

    private void syncOwnedTipObservations() {
        if (ownedTipObservations == null || !tipObservationsDirty) {
            return;
        }
        CanonicalTipObservationAdapter.fillTipObservations(
                tipTree, tipDataModel, passer.getDimension(), ownedTipObservations);
        loadTips(ownedTipObservations);
        tipObservationsDirty = false;
        markJointGradientCacheDirty();
    }

    private void ensureJointGradientCache() {
        if (!jointGradientCacheDirty) {
            return;
        }
        passer.computePostOrderLogLikelihood(transitionProvider, rootPrior);
        passer.computePreOrder(transitionProvider, rootPrior);
        passer.computeJointGradients(transitionProvider, cachedGradientA, cachedGradientQ, cachedGradientMu);
        modelDirty = false;
        jointGradientCacheDirty = false;
        reportTransitionCacheDiagnostics("jointGradient");
    }

    private void markJointGradientCacheDirty() {
        jointGradientCacheDirty = true;
    }

    private void reportTransitionCacheDiagnostics(final String label) {
        if (CACHE_DIAGNOSTIC_REPORT_INTERVAL <= 0) {
            return;
        }
        cacheDiagnosticReportCounter++;
        if (cacheDiagnosticReportCounter % CACHE_DIAGNOSTIC_REPORT_INTERVAL != 0) {
            return;
        }
        if (transitionProvider instanceof CanonicalTransitionCacheDiagnostics) {
            ((CanonicalTransitionCacheDiagnostics) transitionProvider)
                    .reportTransitionCacheDiagnostics(label);
        }
    }

    private static void copyGradient(final double[] source,
                                     final double[] destination,
                                     final String label) {
        if (destination.length != source.length) {
            throw new IllegalArgumentException(
                    "Requested " + label + " gradient length " + destination.length
                            + " does not match cached canonical length " + source.length + ".");
        }
        System.arraycopy(source, 0, destination, 0, source.length);
    }

    private static int inferSelectionGradientDimension(final CanonicalBranchTransitionProvider transitionProvider,
                                                       final CanonicalTreeMessagePasser passer) {
        if (!(transitionProvider instanceof CanonicalOUTransitionProvider)) {
            final int dim = passer.getDimension();
            return dim * dim;
        }

        final OUProcessModel processModel =
                ((CanonicalOUTransitionProvider) transitionProvider).getProcessModel();
        if (processModel.getSelectionMatrixParameterization()
                instanceof OrthogonalBlockCanonicalParameterization) {
            final OrthogonalBlockCanonicalParameterization parameterization =
                    (OrthogonalBlockCanonicalParameterization)
                            processModel.getSelectionMatrixParameterization();
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter =
                    (AbstractBlockDiagonalTwoByTwoMatrixParameter) parameterization.getMatrixParameter();
            if (!(blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider)) {
                throw new IllegalStateException(
                        "Orthogonal block native gradient requires an OrthogonalMatrixProvider rotation parameter.");
            }
            final OrthogonalMatrixProvider orthogonalRotation =
                    (OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter();
            return blockParameter.getBlockDiagonalNParameters()
                    + orthogonalRotation.getOrthogonalParameter().getDimension();
        }

        return processModel.getSelectionMatrixParameterization().getMatrixParameter().getRowDimension()
                * processModel.getSelectionMatrixParameterization().getMatrixParameter().getColumnDimension();
    }

    private static CanonicalTipObservation[] allocateCanonicalTipObservations(final int tipCount, final int dim) {
        final CanonicalTipObservation[] observations = new CanonicalTipObservation[tipCount];
        for (int tipIdx = 0; tipIdx < tipCount; tipIdx++) {
            observations[tipIdx] = new CanonicalTipObservation(dim);
        }
        return observations;
    }
}
