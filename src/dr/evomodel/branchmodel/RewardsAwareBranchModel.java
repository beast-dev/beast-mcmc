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
import dr.evomodel.branchratemodel.RewardMixtureBranchRateModel;
import dr.evomodel.branchratemodel.RewardMixtureCategoryDecoder;
import dr.evomodel.branchratemodel.RewardMixtureCategoryDecoding;
import dr.evomodel.branchratemodel.RewardRates;
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
    public static final boolean DEFAULT_SERICOLA_SERIES_RESCALING = true;

    private final Parameter atomIndices;
    private final RewardRates rewardRates;
    private final SubstitutionModel underlyingSubstitutionModel;
    private final TreeModel tree;
    private final ArbitraryBranchRates branchRateModel;
    private final Parameter indicator;                  // 0/1, same indexing as rewardProportion
    private final RewardMixtureCategoryDecoding categoryDecoder;

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

    private final double[] infinitesimalMatrix;
    private final double[] stateNoJumpLogRate;
    private boolean noJumpLogRatesDirty = true;
    private final double[] atomicScale;
    private boolean atomicScalesDirty = true;
    private final int[] atomicNonZeroIndex;

    @Deprecated
    public RewardsAwareBranchModel(TreeModel tree,
                                   SubstitutionModel underlyingSubstitutionModel,
                                   Parameter rewardRatesValues,
                                   Parameter rewardRatesValuesInternal,
                                   Parameter rewardRatesMapping,
                                   Parameter indicator,
                                   ArbitraryBranchRates branchRateModel,  // TODO? use directly the RewardsAwareMixtureBranchRates instead of the more general ArbitraryBranchRates, to avoid redundant checks and mappings
                                   Parameter atomIndices,
                                   boolean conditional) {
        this(tree,
                underlyingSubstitutionModel,
                rewardRatesValues,
                rewardRatesValuesInternal,
                rewardRatesMapping,
                indicator,
                branchRateModel,
                atomIndices,
                conditional,
                DEFAULT_SERICOLA_SERIES_RESCALING);
    }

    @Deprecated
    public RewardsAwareBranchModel(TreeModel tree,
                                   SubstitutionModel underlyingSubstitutionModel,
                                   Parameter rewardRatesValues,
                                   Parameter rewardRatesValuesInternal,
                                   Parameter rewardRatesMapping,
                                   Parameter indicator,
                                   ArbitraryBranchRates branchRateModel,  // TODO? use directly the RewardsAwareMixtureBranchRates instead of the more general ArbitraryBranchRates, to avoid redundant checks and mappings
                                   Parameter atomIndices,
                                   boolean conditional,
                                   boolean sericolaSeriesRescaling) {
        this(tree,
                underlyingSubstitutionModel,
                new RewardRates(rewardRatesValues, null, rewardRatesValuesInternal, rewardRatesMapping),
                indicator,
                branchRateModel,
                atomIndices,
                conditional,
                sericolaSeriesRescaling);
    }

    public RewardsAwareBranchModel(TreeModel tree,
                                   SubstitutionModel underlyingSubstitutionModel,
                                   RewardRates rewardRates,
                                   Parameter indicator,
                                   ArbitraryBranchRates branchRateModel,  // TODO? use directly the RewardsAwareMixtureBranchRates instead of the more general ArbitraryBranchRates, to avoid redundant checks and mappings
                                   Parameter atomIndices,
                                   boolean conditional) {
        this(tree, underlyingSubstitutionModel, rewardRates, indicator, branchRateModel,
                atomIndices, null, null, null, conditional, DEFAULT_SERICOLA_SERIES_RESCALING);
    }

    public RewardsAwareBranchModel(TreeModel tree,
                                   SubstitutionModel underlyingSubstitutionModel,
                                   RewardRates rewardRates,
                                   Parameter indicator,
                                   ArbitraryBranchRates branchRateModel,  // TODO? use directly the RewardsAwareMixtureBranchRates instead of the more general ArbitraryBranchRates, to avoid redundant checks and mappings
                                   Parameter atomIndices,
                                   boolean conditional,
                                   boolean sericolaSeriesRescaling) {
        this(tree, underlyingSubstitutionModel, rewardRates, indicator, branchRateModel,
                atomIndices, null, null, null, conditional, sericolaSeriesRescaling);
    }

    public RewardsAwareBranchModel(TreeModel tree,
                                   SubstitutionModel underlyingSubstitutionModel,
                                   RewardRates rewardRates,
                                   Parameter categoryParameter,
                                   Parameter categoryCuts,
                                   ArbitraryBranchRates branchRateModel,
                                   boolean conditional) {
        this(tree, underlyingSubstitutionModel, rewardRates, null, branchRateModel,
                null, categoryParameter, categoryCuts, null, conditional, DEFAULT_SERICOLA_SERIES_RESCALING);
    }

    public RewardsAwareBranchModel(TreeModel tree,
                                   SubstitutionModel underlyingSubstitutionModel,
                                   RewardRates rewardRates,
                                   Parameter categoryParameter,
                                   Parameter categoryCuts,
                                   ArbitraryBranchRates branchRateModel,
                                   boolean conditional,
                                   boolean sericolaSeriesRescaling) {
        this(tree, underlyingSubstitutionModel, rewardRates, null, branchRateModel,
                null, categoryParameter, categoryCuts, null, conditional, sericolaSeriesRescaling);
    }

    /**
     * Categorical-state constructor accepting a pre-built decoder (e.g. a
     * PerBranchRewardMixtureCategoryDecoder shared with the branch-rate
     * model that already owns it) instead of building its own
     * RewardMixtureCategoryDecoder internally. Lets the dynamic
     * reward-category-ordering variant reuse one decoder instance -- and
     * one refreshEmbedding() cost -- across both the branch-rate model and
     * this branch model, rather than each maintaining a separate one over
     * the same categoryParameter/categoryCuts.
     */
    public RewardsAwareBranchModel(TreeModel tree,
                                   SubstitutionModel underlyingSubstitutionModel,
                                   RewardRates rewardRates,
                                   RewardMixtureCategoryDecoding externalCategoryDecoder,
                                   ArbitraryBranchRates branchRateModel,
                                   boolean conditional,
                                   boolean sericolaSeriesRescaling) {
        this(tree, underlyingSubstitutionModel, rewardRates, null, branchRateModel,
                null, externalCategoryDecoder.getCategoryParameter(), externalCategoryDecoder.getCutParameter(),
                externalCategoryDecoder, conditional, sericolaSeriesRescaling);
    }

    private RewardsAwareBranchModel(TreeModel tree,
                                    SubstitutionModel underlyingSubstitutionModel,
                                    RewardRates rewardRates,
                                    Parameter indicator,
                                    ArbitraryBranchRates branchRateModel,  // TODO? use directly the RewardsAwareMixtureBranchRates instead of the more general ArbitraryBranchRates, to avoid redundant checks and mappings
                                    Parameter atomIndices,
                                    Parameter categoryParameter,
                                    Parameter categoryCuts,
                                    RewardMixtureCategoryDecoding externalCategoryDecoder,
                                    boolean conditional,
                                    boolean sericolaSeriesRescaling) {

        super(REWARDS_AWARE_BRANCH_MODEL);
        if (tree == null) throw new IllegalArgumentException("tree must be non-null");
        if (underlyingSubstitutionModel == null) {
            throw new IllegalArgumentException("RewardsAwareBranchModel must be provided with an underlying substitution model");
        }
        if (rewardRates == null) throw new IllegalArgumentException("rewardRates must be non-null");
        if (branchRateModel == null) throw new IllegalArgumentException("branchRateModel must be non-null");
        final boolean useCategoricalState = categoryParameter != null || categoryCuts != null;
        if (useCategoricalState) {
            if (categoryParameter == null) throw new IllegalArgumentException("categoryParameter must be non-null");
            if (categoryCuts == null) throw new IllegalArgumentException("categoryCuts must be non-null");
            if (indicator != null || atomIndices != null) {
                throw new IllegalArgumentException(
                        "Provide either categorical mixture state or indicator/atomIndices, not both");
            }
        } else {
            if (indicator == null) throw new IllegalArgumentException("indicator must be non-null");
            if (atomIndices == null) throw new IllegalArgumentException("atomIndices must be non-null");
        }

        this.tree = tree;
        this.underlyingSubstitutionModel = underlyingSubstitutionModel;

        this.branchRateModel = branchRateModel;
        this.indicator = indicator;
        this.atomIndices = atomIndices;
        this.rewardRates = rewardRates;

        final int dim = branchRateModel.getRateParameter().getDimension();
        if (!useCategoricalState) {
            if (indicator.getDimension() != dim) {
                throw new IllegalArgumentException("indicator dimension must equal rewardProportion dimension (branchRateModel rate parameter).");
            }
            if (atomIndices.getDimension() != dim) {
                throw new IllegalArgumentException("atomIndices dimension must equal rewardProportion dimension (branchRateModel rate parameter).");
            }
        }

        this.nstates = underlyingSubstitutionModel.getDataType().getStateCount();
        this.dim2 = nstates * nstates;
        this.categoryDecoder = externalCategoryDecoder != null
                ? externalCategoryDecoder
                : (useCategoricalState
                        ? new RewardMixtureCategoryDecoder(categoryParameter, categoryCuts, nstates, dim)
                        : null);

        final int nodeCount = tree.getNodeCount();
        final int branchCount = nodeCount - 1; // all non-root nodes

        this.W = new double[nodeCount][dim2];
        this.Watomic = new double[nodeCount][dim2];

        this.X = new double[branchCount];
        this.times = new double[branchCount];
        this.Wpacked = new double[branchCount][];

        this.infinitesimalMatrix = new double[dim2];
        this.stateNoJumpLogRate = new double[nstates];
        this.atomicScale = new double[nodeCount];
        this.atomicNonZeroIndex = new int[nodeCount];
        Arrays.fill(this.atomicNonZeroIndex, -1);

        final double epsilon = 1e-10;
        this.sericola = new SericolaSeriesMarkovRewardFastModel( //only for cts values
                underlyingSubstitutionModel,
                rewardRates.getValues(),
                rewardRates.getVaryingValues(),
                rewardRates.getStateIndices(),
                nstates,
                epsilon,
                conditional,
                sericolaSeriesRescaling
        );

        addModel(tree);
        addModel(branchRateModel);
        addModel(sericola);
        if (categoryDecoder == null) {
            addVariable(indicator);
            addVariable(atomIndices);
        } else {
            addVariable(categoryDecoder.getCategoryParameter());
            addVariable(categoryDecoder.getCutParameter());
        }

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
    public int getParameterIndexForNode(final int nodeNr) {
        final NodeRef node = tree.getNode(nodeNr);
        return branchRateModel.getParameterIndexFromNode(node);
    }

    public int getBranchIndexForNodeNumber(int nodeNr) {
        return nodeNrToBranchIndex[nodeNr];
    }
    public Parameter getIndicator() { return indicator; }

    public RewardMixtureCategoryDecoding getCategoryDecoder() { return categoryDecoder; }

    // -------------------- Basic accessors --------------------

    public FrequencyModel getRootFrequencyModel() { return underlyingSubstitutionModel.getFrequencyModel(); }

    @Override
    public SubstitutionModel getRootSubstitutionModel() { return underlyingSubstitutionModel; }

    public TreeModel getTree() { return tree; }

    public RewardRates getRewardRates() { return rewardRates; }

    public BranchRateModel getRateBranchModel() { return branchRateModel; }

    public SericolaSeriesMarkovRewardFastModel getSericolaModel() { return sericola; }

    public boolean isSericolaSeriesRescalingEnabled() { return sericola.isSeriesRescalingEnabled(); }

    public Parameter getRewardRatesValues() { return rewardRates.getValues(); }

    public Parameter getRewardRatesInternal() { return rewardRates.getVaryingValues(); }

    public Parameter getRewardRatesMapping() { return rewardRates.getStateIndices(); }

    public double getContinuousRewardRawForBranch(int branchNodeNumber) {
        final NodeRef node = tree.getNode(branchNodeNumber);
        if (tree.isRoot(node)) {
            throw new IllegalArgumentException("Root node has no branch: " + branchNodeNumber);
        }
        if (branchRateModel instanceof RewardMixtureBranchRateModel) {
            return ((RewardMixtureBranchRateModel) branchRateModel).getContinuousRawReward(tree, node);
        }
        return branchRateModel.getUntransformedBranchRate(tree, node);
    }

    public double getRewardRateRawForState(int stateIndex) {
        if (stateIndex < 0 || stateIndex >= nstates) {
            throw new IllegalArgumentException("stateIndex out of range: " + stateIndex);
        }
        final int rewardRateIndex = (int) Math.round(rewardRates.getStateIndices().getParameterValue(stateIndex));
        if (rewardRateIndex < 0 || rewardRateIndex >= rewardRates.getValues().getDimension()) {
            throw new IllegalArgumentException(
                    "Reward-rate mapping for state " + stateIndex + " points outside rewardRates: " +
                            rewardRateIndex
            );
        }
        return rewardRates.getValues().getParameterValue(rewardRateIndex);
    }

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
        final boolean atomicOn = isAtomicForParameterIndex(p);
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

            final double tol = 0.0; // set higher ( e.g. -1e-12) if we want numerical tolerance
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

        ensureNoJumpLogRatesUpToDate();

        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (tree.isRoot(node)) continue;

            final int nodeNr = node.getNumber();
            final int paramIndex = getParameterIndexForNode(nodeNr);
            if (!isAtomicForParameterIndex(paramIndex)) {
                atomicScale[nodeNr] = 0.0;
                continue;
            }
            final int atomState = getAtomStateForParameterIndex(paramIndex);
            final double t = tree.getBranchLength(node);
            atomicScale[nodeNr] = Math.exp(stateNoJumpLogRate[atomState] * t);
        }

        atomicScalesDirty = false;
    }

    public double[] getTransitionMatrixAtomic(int nodeNr) {
        computeAtomicScales();

        final int paramIndex = getParameterIndexForNode(nodeNr);
        final int atomState = getAtomStateForParameterIndex(paramIndex);
        final int newIndex = atomState * nstates + atomState;

        final double[] matrix = Watomic[nodeNr];
        final int oldIndex = atomicNonZeroIndex[nodeNr];

        if (oldIndex != -1 && oldIndex != newIndex) {
            matrix[oldIndex] = 0.0;
        }

        matrix[newIndex] = atomicScale[nodeNr];
        atomicNonZeroIndex[nodeNr] = newIndex;

        if (DEBUG) {
            if (!(atomicScale[nodeNr] > 0.0)) {
                throw new IllegalStateException(
                        "Atomic transition scale for node " + nodeNr +
                                " (atomState=" + atomState + ", branch length=" +
                                tree.getBranchLength(tree.getNode(nodeNr)) + ")" +
                                " is non-positive: value=" + atomicScale[nodeNr]
                );
            }
        }

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
            final double rate = getContinuousBranchRate(tree, node);

            if (t < 0.0) {
                throw new IllegalArgumentException("Negative branch length for node " + nodeNr + ": " + t);
            }
            if (t == 0.0) {
                setZeroTimeContinuousTransition(nodeNr);
                continue;
            }

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

        if (k == branchIndexToNodeNr.length) {
            // Sericola handles all sorting/lazy caches and writes into ORIGINAL order by design.
            sericola.computePdfInto(X, times, true, Wpacked);
        } else if (k > 0) {
            sericola.computePdfInto(
                    Arrays.copyOf(X, k),
                    Arrays.copyOf(times, k),
                    true,
                    Arrays.copyOf(Wpacked, k));
        }
        ctsMatricesDirty = false;
    }

    private void setZeroTimeContinuousTransition(int nodeNr) {
        // P(X_t = j | X_0 = i, t = 0) = identity, not all-zero: an all-zero
        // matrix is not a valid stochastic transition matrix (rows don't sum
        // to 1) and produces a genuinely invalid (zero) likelihood contribution
        // wherever it's actually used -- see the run README / project log for
        // the ctmc_bm4d_timeseries scenario (c) diagnosis this fixed.
        Arrays.fill(W[nodeNr], 0.0);
        for (int i = 0; i < nstates; i++) {
            W[nodeNr][i * nstates + i] = 1.0;
        }
    }
    public double[] getWPacked(int i) {
        return Wpacked[i];
    }

    public boolean isAtomicBranch(int branchNodeNumber) {
        final int paramIndex = getParameterIndexForNode(branchNodeNumber);
        return isAtomicForParameterIndex(paramIndex);
    }

    public int getAtomicBranchState(int branchNodeNumber) {
        final int paramIndex = getParameterIndexForNode(branchNodeNumber);
        return getAtomStateForParameterIndex(paramIndex);
    }

    public double getAtomicBranchScale(int branchNodeNumber) {
        computeAtomicScales();
        return atomicScale[branchNodeNumber];
    }

    public double getAtomicBranchLogScaleForState(int branchNodeNumber, int stateIndex) {
        final NodeRef node = tree.getNode(branchNodeNumber);
        if (tree.isRoot(node)) {
            throw new IllegalArgumentException("Root node has no branch: " + branchNodeNumber);
        }
        return getAtomicLogScaleForState(stateIndex, tree.getBranchLength(node));
    }

    public double getAtomicLogScaleForState(int stateIndex, double branchLength) {
        if (stateIndex < 0 || stateIndex >= nstates) {
            throw new IllegalArgumentException("stateIndex out of range: " + stateIndex);
        }
        if (branchLength < 0.0) {
            throw new IllegalArgumentException("branchLength must be non-negative: " + branchLength);
        }
        ensureNoJumpLogRatesUpToDate();
        return stateNoJumpLogRate[stateIndex] * branchLength;
    }

    private void ensureNoJumpLogRatesUpToDate() {
        if (!noJumpLogRatesDirty) return;

        underlyingSubstitutionModel.getInfinitesimalMatrix(infinitesimalMatrix);
        for (int state = 0; state < nstates; state++) {
            final double qii = infinitesimalMatrix[state * nstates + state];
            if (Double.isNaN(qii) || Double.isInfinite(qii) || qii > 1.0e-12) {
                throw new IllegalStateException(
                        "Invalid infinitesimal diagonal for no-jump atom at state " + state + ": " + qii);
            }
            stateNoJumpLogRate[state] = qii > 0.0 ? 0.0 : qii;
        }

        noJumpLogRatesDirty = false;
    }

    private int getAtomStateForParameterIndex(int paramIndex) {
        if (categoryDecoder != null) {
            return categoryDecoder.getAtomicState(paramIndex);
        }
        final double raw = atomIndices.getParameterValue(paramIndex);
        final int state = (int) Math.round(raw);
        if (Math.abs(raw - state) > 1.0e-9 || state < 0 || state >= nstates) {
            throw new IllegalArgumentException(
                    "atomIndices must contain integer state indices in [0, " + (nstates - 1) +
                            "], found " + raw + " at parameter index " + paramIndex);
        }
        return state;
    }

    private boolean isAtomicForParameterIndex(final int paramIndex) {
        if (categoryDecoder != null) {
            return categoryDecoder.isAtomic(paramIndex);
        }
        return isOne(indicator.getParameterValue(paramIndex));
    }

    private double getContinuousBranchRate(final TreeModel tree, final NodeRef node) {
        if (branchRateModel instanceof RewardMixtureBranchRateModel) {
            final RewardMixtureBranchRateModel rewardMixtureRates =
                    (RewardMixtureBranchRateModel) branchRateModel;
            return rewardMixtureRates.getBranchRateForRawReward(
                    tree,
                    node,
                    rewardMixtureRates.getContinuousRawReward(tree, node));
        }
        return branchRateModel.getBranchRate(tree, node);
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
        if (model == sericola) {
            invalidateNoJumpLogRates();
            invalidateAtomicScales();
            invalidateCtsMatrices();
        }
        if (model == tree) {
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
        if (variable == atomIndices) {
            invalidateAtomicScales();
//            invalidateCtsMatrices(); // safe, even if a bit conservative
        } else if (categoryDecoder != null && variable == categoryDecoder.getCategoryParameter()) {
            invalidateAtomicScales();
        } else if (categoryDecoder != null && variable == categoryDecoder.getCutParameter()) {
            categoryDecoder.refreshEmbedding();
            invalidateAtomicScales();
        }
        fireModelChanged();
    }

    // -------------------- MCMC store/restore --------------------

    private boolean storedKnownTransitionMatrices;

    private boolean ctsMatricesDirty = true;
    private boolean ctsMatricesDirtyDuringProposal = false;
    private boolean storedCtsMatricesDirty;

    private boolean atomicScalesDirtyDuringProposal = false;
    private boolean storedAtomicScalesDirty;

    private boolean noJumpLogRatesDirtyDuringProposal = false;
    private boolean storedNoJumpLogRatesDirty;

    private void invalidateCtsMatrices() {
        ctsMatricesDirty = true;
        ctsMatricesDirtyDuringProposal = true;
    }
    private void invalidateAtomicScales() {
        atomicScalesDirty = true;
        atomicScalesDirtyDuringProposal = true;
    }
    private void invalidateNoJumpLogRates() {
        noJumpLogRatesDirty = true;
        noJumpLogRatesDirtyDuringProposal = true;
    }

    @Override
    protected void storeState() {
        storedCtsMatricesDirty = ctsMatricesDirty;
        ctsMatricesDirtyDuringProposal = false;

        storedAtomicScalesDirty = atomicScalesDirty;
        atomicScalesDirtyDuringProposal = false;

        storedNoJumpLogRatesDirty = noJumpLogRatesDirty;
        noJumpLogRatesDirtyDuringProposal = false;
    }

    @Override
    protected void restoreState() {
        if (categoryDecoder != null) {
            categoryDecoder.refreshEmbedding();
        }
        // If continuous inputs were touched during the rejected proposal,
        // the cached W arrays may contain proposal values, so force recomputation.
        ctsMatricesDirty = storedCtsMatricesDirty || ctsMatricesDirtyDuringProposal;
        ctsMatricesDirtyDuringProposal = false;

        atomicScalesDirty = storedAtomicScalesDirty || atomicScalesDirtyDuringProposal;
        atomicScalesDirtyDuringProposal = false;

        noJumpLogRatesDirty = storedNoJumpLogRatesDirty || noJumpLogRatesDirtyDuringProposal;
        noJumpLogRatesDirtyDuringProposal = false;
    }

    @Override
    protected void acceptState() {
        ctsMatricesDirtyDuringProposal = false;
        atomicScalesDirtyDuringProposal = false;
        noJumpLogRatesDirtyDuringProposal = false;
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
