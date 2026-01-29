/*
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

package dr.evomodel.branchmodel;

import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.substmodel.*;
import dr.evomodel.tree.TreeModel;
import dr.inference.markovjumps.SericolaSeriesMarkovRewardFastModel;
import dr.inference.model.*;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * TransitionMatrixProviderBranchModel that uses SericolaSeriesMarkovRewardFastModel
 * to compute reward-conditioned transition matrices on each branch.
 *
 * Design goals:
 *  - BranchModel only:
 *      (i) packs branch rewards X and times,
 *      (ii) provides output buffers W[nodeNr],
 *      (iii) caches "knownTransitionMatrices" at the branch-model level.
 *
 * @author Filippo Monti
 */
public class RewardsAwareBranchModel extends AbstractModel
        implements TransitionMatrixProviderBranchModel, Citable, Reportable {

    public static final String REWARDS_AWARE_BRANCH_MODEL = "RewardsAwareBranchModel";

    private final Parameter rewardRates;
    private final SubstitutionModel underlyingSubstitutionModel;
    private final TreeModel tree;
    private final ArbitraryBranchRates branchRateModel;

    private final int nstates;
    private final int dim2;

    // Output transition matrices in ORIGINAL state order, per node index.
    // Root row is unused but kept for direct nodeNr indexing.
    private final double[][] W;

    // Packed inputs for Sericola: one entry per non-root branch (nodeCount - 1)
    private final double[] X;          // reward = (branchRate * branchLength)
    private final double[] times;      // branchLength
    private final double[][] Wpacked;  // references into W[nodeNr], one per branch

    // Sericola model (does all lazy updates and store/restore internally)
    private final SericolaSeriesMarkovRewardFastModel sericola;

    // Cache flag at this branch-model layer
    private boolean knownTransitionMatrices = false;

    // For compatibility with SubstitutionModelDelegate
    private List<SubstitutionModel> substitutionModels;

    // Optional debug/testing
    private boolean DUMMYTESTING = false;
    private boolean DEBUG = false;

    public RewardsAwareBranchModel(TreeModel tree,
                                   SubstitutionModel underlyingSubstitutionModel,
                                   Parameter rewardRates,
                                   ArbitraryBranchRates branchRateModel) {

        super(REWARDS_AWARE_BRANCH_MODEL);

        if (tree == null) throw new IllegalArgumentException("tree must be non-null");
        if (underlyingSubstitutionModel == null) {
            throw new IllegalArgumentException("RewardsAwareBranchModel must be provided with an underlying substitution model");
        }
        if (rewardRates == null) throw new IllegalArgumentException("rewardRates must be non-null");
        if (branchRateModel == null) throw new IllegalArgumentException("branchRateModel must be non-null");

        this.tree = tree;
        this.underlyingSubstitutionModel = underlyingSubstitutionModel;
        this.rewardRates = rewardRates;
        this.branchRateModel = branchRateModel;

        this.nstates = underlyingSubstitutionModel.getDataType().getStateCount();
        this.dim2 = nstates * nstates;

        final int nodeCount = tree.getNodeCount();
        final int branchCount = nodeCount - 1; // all non-root nodes

        this.W = new double[nodeCount][dim2];

        this.X = new double[branchCount];
        this.times = new double[branchCount];
        this.Wpacked = new double[branchCount][];

        // Create Sericola *as a BEAST model* so it listens to Q and rewardRates itself.
        // Choose epsilon to match your desired truncation accuracy.
        final double epsilon = 1e-10;
        this.sericola = new SericolaSeriesMarkovRewardFastModel(
                underlyingSubstitutionModel,
                rewardRates,
                nstates,
                epsilon
        );

        // Register dependencies:
        // We listen to tree + branchRateModel to know X/times changed,
        // and to sericola for reward/Q changes (already hooked internally).
        addModel(tree);
        addModel(branchRateModel);
        addModel(sericola);

        // NOTE: rewardRates & underlyingSubstitutionModel are already dependencies of sericola;
        // we do NOT add them again here to avoid duplicating event graphs.
    }

    // -------------------- Basic accessors --------------------

    public FrequencyModel getRootFrequencyModel() { return underlyingSubstitutionModel.getFrequencyModel(); }

    @Override
    public SubstitutionModel getRootSubstitutionModel() { return underlyingSubstitutionModel; }

    public TreeModel getTree() { return tree; }

    public Parameter getRewardRates() { return rewardRates; }

    public BranchRateModel getRateBranchModel() { return branchRateModel; }

    // -------------------- Main API --------------------

    @Override
    public double[] getTransitionMatrix(NodeRef branch) {
        return getTransitionMatrix(branch.getNumber());
    }

    public double[] getTransitionMatrix(int nodeNr) {
        if (DUMMYTESTING) {
            NodeRef node = tree.getNode(nodeNr);
            double t = tree.getBranchLength(node);
            getRootSubstitutionModel().getTransitionProbabilities(t, W[nodeNr]);
            return W[nodeNr];
        }
        computeTransitionMatrices();
        return W[nodeNr];
    }

    private void computeTransitionMatrices() {

        if (knownTransitionMatrices) return;
        if (DEBUG) System.out.println("RewardsAwareBranchModel: computeTransitionMatrices");

        final int nodeCount = tree.getNodeCount();

        int k = 0;
        for (int i = 0; i < nodeCount; i++) {
            final NodeRef node = tree.getNode(i);
            if (tree.isRoot(node)) continue;

            final int nodeNr = node.getNumber();
            final double t = tree.getBranchLength(node);
            final double rate = branchRateModel.getBranchRate(tree, node);

            times[k] = t;
            X[k] = rate * t;

            // Write results directly into W[nodeNr] (original order)
            Wpacked[k] = W[nodeNr];
            k++;
        }

        // Sericola handles all sorting/lazy caches and writes into ORIGINAL order by design
        sericola.computePdfInto(X, times, true, Wpacked);

        knownTransitionMatrices = true;
    }

    // -------------------- Branch model mapping --------------------

    @Override
    public Mapping getBranchModelMapping(NodeRef node) {
        final double[] weights = new double[]{1.0};
        final int[] order = new int[]{node.getNumber()};

        return new Mapping() {
            @Override
            public int[] getOrder() { return order; }

            @Override
            public double[] getWeights() { return weights; }
        };
    }

    @Override
    public boolean requiresMatrixConvolution() {
        return false;
    }

    // -------------------- Model events --------------------

    private boolean ignoreModelChangedEvent = false;

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (ignoreModelChangedEvent) return;

        // If any dependency changes, our W cache is invalid.
        // We do NOT manage sericola dirty flags here; sericola manages itself.
        knownTransitionMatrices = false;
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // We only added tree/branchRateModel/sericola as dependencies.
        // If you later add variables directly to this model, handle them here.
        knownTransitionMatrices = false;
        fireModelChanged();
    }

    // -------------------- MCMC store/restore --------------------

    private boolean storedKnownTransitionMatrices;

    @Override
    protected void storeState() {
        storedKnownTransitionMatrices = knownTransitionMatrices;
    }

    @Override
    protected void restoreState() {
        // Safest: force recomputation on demand.
        // (Even if storedKnownTransitionMatrices was true, the underlying state may have reverted.)
        knownTransitionMatrices = false;
    }

    @Override
    protected void acceptState() {
        // nothing
    }

    // -------------------- Citable / Reportable --------------------

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SUBSTITUTION_MODELS;
    }

    @Override
    public String getDescription() {
        return "Rewards Aware Branch model";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(
                new Citation(new Author[]{new Author("F", "Monti"),
                        new Author("MA", "Suchard")},
                        "Dependencies between CTMCs",
                        2025,
                        "TOBE",
                        1,
                        1,
                        1,
                        Citation.Status.IN_PRESS));
    }

    @Override
    public String getReport() {
        if (!DUMMYTESTING) computeTransitionMatrices();
        StringBuilder sb = new StringBuilder();
        sb.append("W matrix: ");
        for (int nodeNr = 0; nodeNr < tree.getNodeCount(); nodeNr++) {
            NodeRef node = tree.getNode(nodeNr);
            if (tree.isRoot(node)) continue;
            for (double val : W[nodeNr]) sb.append(val).append(" ");
        }
        return sb.toString();
    }

    // -------------------- Compatibility: SubstitutionModelDelegate --------------------

    @Override
    public List<SubstitutionModel> getSubstitutionModels() {
        if (substitutionModels == null) {
            buildSubstitutionModels();
        }
        return substitutionModels;
    }

    protected void buildSubstitutionModels() {

        substitutionModels = new ArrayList<>();

        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (tree.isRoot(node)) continue;

            ignoreModelChangedEvent = true;
            SubstitutionModel substitutionModel = new TransitionMatrixProvider(
                    "RewardsAwareSubstitutionModel",
                    underlyingSubstitutionModel.getDataType(),
                    underlyingSubstitutionModel.getFrequencyModel(),
                    W[node.getNumber()]
            );
            ignoreModelChangedEvent = false;

            substitutionModels.add(substitutionModel);
        }
    }

    class TransitionMatrixProvider extends ComplexSubstitutionModel {
        private final double[] transitionMatrixRef;

        public TransitionMatrixProvider(String name, DataType dataType,
                                        FrequencyModel freqModel, double[] transitionMatrixRef) {
            super(name, dataType, freqModel, null);
            this.transitionMatrixRef = transitionMatrixRef;
        }

        @Override
        public void getTransitionProbabilities(double distance, double[] matrix) {
            // Ensure current before exposing.
            if (!DUMMYTESTING) computeTransitionMatrices();
            System.arraycopy(transitionMatrixRef, 0, matrix, 0, transitionMatrixRef.length);
        }

        @Override protected void handleModelChangedEvent(Model model, Object object, int index) {}
        @Override protected void frequenciesChanged() {}
        @Override protected void ratesChanged() {}
        @Override protected void setupRelativeRates(double[] rates) {}
    }
}
