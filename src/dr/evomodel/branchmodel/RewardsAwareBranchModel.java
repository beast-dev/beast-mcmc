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
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * TransitionMatrixProviderBranchModel that uses SericolaSeriesMarkovRewardFastModelRho
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
    private final Parameter indicator;                  // 0/1, same indexing as rho
    private final IndexedParameter indexedParameter;    // same indexing as rho

    private final int nstates;
    private final int dim2;

    // Output transition matrices in ORIGINAL state order, per node index.
    // Root row is unused but kept for direct nodeNr indexing.
    private final double[][] W;
    private final double[][] Watomic;

    // Packed inputs for Sericola: one entry per non-root branch (nodeCount - 1)
    private final double[] X;          // total reward
    private final double[] times;      // branchLength
    private final double[][] Wpacked;  // references into W[nodeNr], one per branch

    private final SericolaSeriesMarkovRewardFastModel sericola;

    // Cache flag at this branch-model layer

    // For compatibility with SubstitutionModelDelegate
    private List<SubstitutionModel> substitutionModels;

    // Optional debug/testing
    private boolean DUMMYTESTING = false;
    private boolean DEBUG = false;

    private final int[] branchIndexToNodeNr;
    private final int[] nodeNrToBranchIndex;

    private final double[] atomicScale;
    private boolean atomicScalesDirty = true;
    private final int[] atomicNonZeroIndex;

    public RewardsAwareBranchModel(TreeModel tree,
                                   SubstitutionModel underlyingSubstitutionModel,
                                   Parameter rewardRates,
                                   Parameter indicator,
                                   IndexedParameter indexedParameter,
                                   ArbitraryBranchRates branchRateModel,
                                   Parameter extremeIndices,
                                   boolean conditional) {

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
        this.indicator = indicator;
        this.indexedParameter = indexedParameter;

        final int dim = branchRateModel.getRateParameter().getDimension();
        if (indicator.getDimension() != dim) {
            throw new IllegalArgumentException("indicator dimension must equal rho dimension (branchRateModel rate parameter).");
        }
        if (indexedParameter.getDimension() != dim) {
            throw new IllegalArgumentException("indexedParameter dimension must equal rho dimension (branchRateModel rate parameter).");
        }


        this.nstates = underlyingSubstitutionModel.getDataType().getStateCount();
        this.dim2 = nstates * nstates;

        final int nodeCount = tree.getNodeCount();
        final int branchCount = nodeCount - 1; // all non-root nodes

        this.W = new double[nodeCount][dim2];
        this.Watomic = new double[nodeCount][dim2];

        this.X = new double[branchCount];
        this.times = new double[branchCount];
        this.Wpacked = new double[branchCount][];

        this.atomicScale = new double[nodeCount];
        this.atomicNonZeroIndex = new int[nodeCount];
        Arrays.fill(this.atomicNonZeroIndex, -1);

        final double epsilon = 1e-10;
        this.sericola = new SericolaSeriesMarkovRewardFastModel( //only for cts values
                underlyingSubstitutionModel,
                rewardRates,
                nstates,
                epsilon,
                conditional
        );

        addModel(tree);
        addModel(branchRateModel);
        addModel(sericola);
        addVariable(indicator);
        addVariable(indexedParameter.getIndicesParameter());
        addVariable(indexedParameter.getValuesParameter());
        addVariable(extremeIndices);

        final int nNodes = tree.getNodeCount();
        final int nBranches = nNodes - 1;

        branchIndexToNodeNr = new int[nBranches];
        nodeNrToBranchIndex = new int[nNodes];
        Arrays.fill(nodeNrToBranchIndex, -1);

        int k = 0;
        for (int i = 0; i < nNodes; i++) {
            NodeRef node = tree.getNode(i);
            if (tree.isRoot(node)) continue;

            int nodeNr = node.getNumber();
            branchIndexToNodeNr[k] = nodeNr;
            nodeNrToBranchIndex[nodeNr] = k;
            k++;
        }
    }
    public int getNodeNumberForBranchIndex(int branchIndex) {
        return branchIndexToNodeNr[branchIndex];
    }
    private int getParameterIndexForNode(final int nodeNr) {
        final NodeRef node = tree.getNode(nodeNr);
        // This method name is from your commented code; if your version differs, swap accordingly.
        return branchRateModel.getParameterIndexFromNode(node);
    }

    public int getBranchIndexForNodeNumber(int nodeNr) {
        return nodeNrToBranchIndex[nodeNr];
    }
    public Parameter getIndicator() { return indicator; }
    public IndexedParameter getIndexedParameter() { return indexedParameter; }


    // -------------------- Basic accessors --------------------

    public FrequencyModel getRootFrequencyModel() { return underlyingSubstitutionModel.getFrequencyModel(); }

    @Override
    public SubstitutionModel getRootSubstitutionModel() { return underlyingSubstitutionModel; }

    public TreeModel getTree() { return tree; }

    public Parameter getRewardRates() { return rewardRates; }

    public BranchRateModel getRateBranchModel() { return branchRateModel; }

    public SericolaSeriesMarkovRewardFastModel getSericolaModel() { return sericola; }

    public double getUniformizationRate() {
        return sericola.getUniformizationRate();
    }

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
        final int p = getParameterIndexForNode(nodeNr);
        final boolean atomicOn = isOne(indicator.getParameterValue(p));
        if (atomicOn) {
            return getTransitionMatrixAtomic(nodeNr);
        } else {
            return getTransitionMatrixCts(nodeNr);
        }
    }

    public double[] getTransitionMatrixCts(int nodeNr) {
        computeCtsTransitionMatrices();
        if (DEBUG) {
            final double[] w = W[nodeNr];

            final double tol = 0.0; // or something like -1e-12 if you want numerical tolerance
            for (int i = 0; i < w.length; i++) {
                if (!(w[i] > tol)) {
                    throw new IllegalStateException(
                            "Transition matrix for node " + nodeNr +
                                    " contains non-positive entry at index " + i +
                                    ": value=" + w[i]
                    );
                }
            }
        }
        return W[nodeNr];
    }

    private void computeAtomicScales() {
        if (!atomicScalesDirty) return;

        final double lambda = sericola.getLambda();

        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (tree.isRoot(node)) continue;

            final int nodeNr = node.getNumber();
            final double t = tree.getBranchLength(node);
            atomicScale[nodeNr] = Math.exp(-lambda * t);
        }

        atomicScalesDirty = false;
    }

    public double[] getTransitionMatrixAtomic(int nodeNr) {
        computeAtomicScales();

        final int paramIndex = getParameterIndexForNode(nodeNr);
        final int atomState = indexedParameter.getIndexValue(paramIndex);
        final int newIndex = atomState * nstates + atomState;

        final double[] matrix = Watomic[nodeNr];
        final int oldIndex = atomicNonZeroIndex[nodeNr];

        if (oldIndex != -1 && oldIndex != newIndex) {
            matrix[oldIndex] = 0.0;
        }

        matrix[newIndex] = atomicScale[nodeNr];
        atomicNonZeroIndex[nodeNr] = newIndex;

        return matrix;
    }

    private static boolean isOne(final double x) {
        final long r = Math.round(x);
        return (Math.abs(x - r) <= 1e-9) && (r == 1L);
    }

    private void computeCtsTransitionMatrices() {

//        if (knownTransitionMatrices) return;
        if (!ctsMatricesDirty) return;
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
            X[k] = rate;

            // Write results directly into W[nodeNr] (original order)
            Wpacked[k] = W[nodeNr];
            k++;
        }

        if (false) {
            System.err.print("times = ");
            for (int i = 0; i < k; i++) {
                System.err.print(times[i] - 1);
                if (i < k - 1) System.err.print(" ");
            }
            System.err.println("]");
            throw new RuntimeException("finished printing");
        }

        // Sericola handles all sorting/lazy caches and writes into ORIGINAL order by design
//        sericola.computePdfIntoY(X, times, true, Wpacked);
        sericola.computePdfInto(X, times, true, Wpacked);
        ctsMatricesDirty = false;
    }
    public double[] getWPacked(int i) {
        return Wpacked[i];
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
        if (model == sericola || model == tree) {
            invalidateAtomicScales();
            invalidateCtsMatrices();
        }
        if (model == branchRateModel) {
            invalidateCtsMatrices();
        }

        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
    }

    // -------------------- MCMC store/restore --------------------

    private boolean storedKnownTransitionMatrices;

    private boolean ctsMatricesDirty = true;
    private boolean ctsMatricesDirtyDuringProposal = false;
    private boolean storedCtsMatricesDirty;

    private boolean atomicScalesDirtyDuringProposal = false;
    private boolean storedAtomicScalesDirty;

    private void invalidateCtsMatrices() {
        ctsMatricesDirty = true;
        ctsMatricesDirtyDuringProposal = true;
    }
    private void invalidateAtomicScales() {
        atomicScalesDirty = true;
        atomicScalesDirtyDuringProposal = true;
    }

    @Override
    protected void storeState() {
        storedCtsMatricesDirty = ctsMatricesDirty;
        ctsMatricesDirtyDuringProposal = false;

        storedAtomicScalesDirty = atomicScalesDirty;
        atomicScalesDirtyDuringProposal = false;
    }

    @Override
    protected void restoreState() {
        // If continuous inputs were touched during the rejected proposal,
        // the cached W arrays may contain proposal values, so force recomputation.
        ctsMatricesDirty = storedCtsMatricesDirty || ctsMatricesDirtyDuringProposal;
        ctsMatricesDirtyDuringProposal = false;

        atomicScalesDirty = storedAtomicScalesDirty || atomicScalesDirtyDuringProposal;
        atomicScalesDirtyDuringProposal = false;
    }

    @Override
    protected void acceptState() {
        ctsMatricesDirtyDuringProposal = false;
        atomicScalesDirtyDuringProposal = false;
    }

    public int getStateCount() {
        return nstates;
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
                        2026,
                        "TOBE",
                        1,
                        1,
                        1,
                        Citation.Status.IN_PRESS));
    }

    @Override
    public String getReport() {
        if (!DUMMYTESTING) computeCtsTransitionMatrices();
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

    // this method is here only for compatibility.
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
                    node.getNumber()
            );
            ignoreModelChangedEvent = false;

            substitutionModels.add(substitutionModel);
            return;  // returning early to avoid allocations
        }
    }

    class TransitionMatrixProvider extends ComplexSubstitutionModel {
        private final int nodeNr;

        public TransitionMatrixProvider(String name,
                                        DataType dataType,
                                        FrequencyModel freqModel,
                                        int nodeNr) {
            super(name, dataType, freqModel, null);
            this.nodeNr = nodeNr;
        }

        @Override
        public void getTransitionProbabilities(double distance, double[] matrix) {
            final double[] source;

            if (DUMMYTESTING) {
                NodeRef node = tree.getNode(nodeNr);
                double t = tree.getBranchLength(node);
                getRootSubstitutionModel().getTransitionProbabilities(t, matrix);
                return;
            } else {
                source = RewardsAwareBranchModel.this.getTransitionMatrix(nodeNr);
            }
            System.arraycopy(source, 0, matrix, 0, source.length);
        }

        @Override protected void handleModelChangedEvent(Model model, Object object, int index) {}
        @Override protected void frequenciesChanged() {}
        @Override protected void ratesChanged() {}
        @Override protected void setupRelativeRates(double[] rates) {}
    }
}
