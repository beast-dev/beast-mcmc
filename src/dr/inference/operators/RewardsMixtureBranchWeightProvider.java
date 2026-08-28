package dr.inference.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DiscreteDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;

import java.util.Arrays;

/**
 * Computes exact branch-local reward-mixture weights.
 *
 * For a branch, category 0 is the continuous reward candidate and categories
 * 1..K are atomic no-jump candidates for states 0..K-1.
 *
 * The primary reward-process likelihood handled through the Java discrete
 * delegate is currently intended for one driving process at a time. Dependent
 * CTMC contributions are added separately through RewardDependentEdgeEvidenceProvider
 * instances and should stay BEAGLE-backed/pattern-complete.
 */
public final class RewardsMixtureBranchWeightProvider {

    private final RewardsAwareBranchModel rewardsAwareBranchModel;
    private final TreeDataLikelihood treeDataLikelihood;
    private final DiscreteDataLikelihoodDelegate discreteDelegate;
    private final RewardDependentEdgeEvidenceProvider[] dependentEvidenceProviders;
    private final Tree tree;
    private final int nstates;
    private final int branchCount;

    private final int[] nodeNumberByParameterIndex;
    private final double[] prePartial;
    private final double[] postPartial;
    private final double[] preScales;
    private final double[] postScales;
    private final double[] logAtomicWeights;
    private final double[] continuousMatrixScratch;
    private final double[] singleCtsRewardScratch;
    private final double[] singleCtsTimeScratch;
    private final double[][] singleCtsMatrixScratch;
    private final double[][] cachedPrePartialByParameterIndex;
    private final double[][] cachedPostPartialByParameterIndex;
    private final double[] cachedPreScaleByParameterIndex;
    private final double[] cachedPostScaleByParameterIndex;
    private final double[] cachedSinglePrePartial;
    private final double[] cachedSinglePostPartial;
    private final RewardsMixtureBranchResamplingHelper.BranchWeights[] cachedBranchWeightsByParameterIndex;
    private final int[] cachedBranchWeightsEpoch;

    private boolean operationCacheActive = false;
    private boolean cacheAllBranchMessages = false;
    private boolean cacheSingleBranchMessages = false;
    private boolean likelihoodMessagesFresh = false;
    private int singleBranchCacheParameterIndex = -1;
    private double cachedSinglePreScale = Double.NaN;
    private double cachedSinglePostScale = Double.NaN;
    private int operationCacheEpoch = 1;

    public RewardsMixtureBranchWeightProvider(final RewardsAwareBranchModel rewardsAwareBranchModel,
                                              final TreeDataLikelihood treeDataLikelihood,
                                              final TreeDataLikelihood[] dependentTreeDataLikelihoods,
                                              final TreeDataLikelihood[] dependentContinuousTreeDataLikelihoods) {
        this(
                rewardsAwareBranchModel,
                treeDataLikelihood,
                requireDiscreteDelegate(treeDataLikelihood),
                createDependentEvidenceProviders(
                        rewardsAwareBranchModel,
                        dependentTreeDataLikelihoods,
                        dependentContinuousTreeDataLikelihoods));
    }

    public RewardsMixtureBranchWeightProvider(final RewardsAwareBranchModel rewardsAwareBranchModel,
                                              final TreeDataLikelihood treeDataLikelihood,
                                              final DiscreteDataLikelihoodDelegate discreteDelegate,
                                              final RewardDependentEdgeEvidenceProvider[] dependentEvidenceProviders) {
        if (rewardsAwareBranchModel == null) {
            throw new IllegalArgumentException("rewardsAwareBranchModel must be non-null");
        }
        if (treeDataLikelihood == null) {
            throw new IllegalArgumentException("treeDataLikelihood must be non-null");
        }
        if (discreteDelegate == null) {
            throw new IllegalArgumentException("discreteDelegate must be non-null");
        }
        if (treeDataLikelihood.getTree().getNodeCount() != rewardsAwareBranchModel.getTree().getNodeCount()) {
            throw new IllegalArgumentException(
                    "TreeDataLikelihood and RewardsAwareBranchModel must use trees with the same node count.");
        }

        this.rewardsAwareBranchModel = rewardsAwareBranchModel;
        this.treeDataLikelihood = treeDataLikelihood;
        this.discreteDelegate = discreteDelegate;
        this.dependentEvidenceProviders = dependentEvidenceProviders == null
                ? new RewardDependentEdgeEvidenceProvider[0]
                : Arrays.copyOf(dependentEvidenceProviders, dependentEvidenceProviders.length);
        this.tree = rewardsAwareBranchModel.getTree();
        this.nstates = rewardsAwareBranchModel.getStateCount();
        this.branchCount = tree.getNodeCount() - 1;

        this.nodeNumberByParameterIndex = new int[branchCount];
        initializeBranchMappings();

        this.prePartial = new double[nstates];
        this.postPartial = new double[nstates];
        this.preScales = new double[discreteDelegate.getPatternCount()];
        this.postScales = new double[discreteDelegate.getPatternCount()];
        this.logAtomicWeights = new double[nstates];
        this.continuousMatrixScratch = new double[nstates * nstates];
        this.singleCtsRewardScratch = new double[1];
        this.singleCtsTimeScratch = new double[1];
        this.singleCtsMatrixScratch = new double[1][];
        this.cachedPrePartialByParameterIndex = new double[branchCount][nstates];
        this.cachedPostPartialByParameterIndex = new double[branchCount][nstates];
        this.cachedPreScaleByParameterIndex = new double[branchCount];
        this.cachedPostScaleByParameterIndex = new double[branchCount];
        this.cachedSinglePrePartial = new double[nstates];
        this.cachedSinglePostPartial = new double[nstates];
        this.cachedBranchWeightsByParameterIndex =
                new RewardsMixtureBranchResamplingHelper.BranchWeights[branchCount];
        this.cachedBranchWeightsEpoch = new int[branchCount];
    }

    public Tree getTree() {
        return tree;
    }

    public int getStateCount() {
        return nstates;
    }

    public int getCategoryCount() {
        return nstates + 1;
    }

    public int getBranchCount() {
        return branchCount;
    }

    public int getParameterIndexForNode(final int branchNodeNumber) {
        return rewardsAwareBranchModel.getParameterIndexForNode(branchNodeNumber);
    }

    public int getNodeNumberForParameterIndex(final int parameterIndex) {
        if (parameterIndex < 0 || parameterIndex >= branchCount) {
            throw new IllegalArgumentException("Branch parameter index out of range: " + parameterIndex);
        }
        return nodeNumberByParameterIndex[parameterIndex];
    }

    public void refreshRewardCategoryEmbedding() {
        refreshRewardCategoryEmbedding(
                RewardMixturePerformanceStats.OperationCacheClearReason.DYNAMIC_EMBEDDING_REFRESH);
    }

    public void refreshRewardCategoryEmbedding(
            final RewardMixturePerformanceStats.OperationCacheClearReason reason) {
        rewardsAwareBranchModel.refreshCategoryDecoderEmbedding();
        clearOperationCache(reason);
    }

    public void beginOperationCache() {
        operationCacheActive = true;
        cacheAllBranchMessages = true;
        cacheSingleBranchMessages = false;
        singleBranchCacheParameterIndex = -1;
        RewardMixturePerformanceStats.recordOperationCacheBegin();
        clearOperationCache(RewardMixturePerformanceStats.OperationCacheClearReason.OPERATION_START);
    }

    public void beginSingleBranchOperationCache(final int parameterIndex) {
        getNodeNumberForParameterIndex(parameterIndex);
        operationCacheActive = true;
        cacheAllBranchMessages = false;
        cacheSingleBranchMessages = true;
        singleBranchCacheParameterIndex = parameterIndex;
        RewardMixturePerformanceStats.recordOperationCacheBegin();
        clearOperationCache(RewardMixturePerformanceStats.OperationCacheClearReason.OPERATION_START);
    }

    public void clearOperationCache() {
        clearOperationCache(RewardMixturePerformanceStats.OperationCacheClearReason.UNKNOWN);
    }

    public void clearOperationCache(final RewardMixturePerformanceStats.OperationCacheClearReason reason) {
        RewardMixturePerformanceStats.recordOperationCacheClear(reason);
        likelihoodMessagesFresh = false;
        advanceOperationCacheEpoch();
    }

    public void refreshLikelihoodMessages() {
        final long refreshStart = RewardMixturePerformanceStats.startTimer();
        try {
            long start = RewardMixturePerformanceStats.startTimer();
            discreteDelegate.updatePostOrdersFromTreeDataLikelihood(treeDataLikelihood);
            RewardMixturePerformanceStats.recordPostOrderRefresh(
                    RewardMixturePerformanceStats.elapsed(start));

            start = RewardMixturePerformanceStats.startTimer();
            discreteDelegate.ensurePreOrderComputed();
            RewardMixturePerformanceStats.recordPreOrderRefresh(
                    RewardMixturePerformanceStats.elapsed(start));

            for (RewardDependentEdgeEvidenceProvider provider : dependentEvidenceProviders) {
                start = RewardMixturePerformanceStats.startTimer();
                provider.prepare();
                RewardMixturePerformanceStats.recordDependentPrepare(
                        RewardMixturePerformanceStats.elapsed(start));
            }
            if (operationCacheActive) {
                if (cacheAllBranchMessages) {
                    fillBranchMessageCache();
                } else if (cacheSingleBranchMessages) {
                    fillSingleBranchMessageCache(singleBranchCacheParameterIndex);
                }
                likelihoodMessagesFresh = true;
                advanceOperationCacheEpoch();
            }
        } finally {
            RewardMixturePerformanceStats.recordRewardProviderRefresh(
                    RewardMixturePerformanceStats.elapsed(refreshStart));
        }
    }

    public RewardsMixtureBranchResamplingHelper.BranchWeights computeBranchWeightsForParameterIndex(
            final int parameterIndex) {
        return computeBranchWeightsForNode(getNodeNumberForParameterIndex(parameterIndex));
    }

    public double computeContinuousLogWeightForParameterIndex(final int parameterIndex,
                                                             final double rawContinuousReward) {
        if (!operationCacheActive) {
            RewardMixturePerformanceStats.recordBranchWeightCacheBypass();
            refreshLikelihoodMessages();
        } else if (!likelihoodMessagesFresh) {
            RewardMixturePerformanceStats.recordLikelihoodMessageMiss();
            refreshLikelihoodMessages();
        }
        return computeContinuousLogWeightForNode(
                getNodeNumberForParameterIndex(parameterIndex),
                rawContinuousReward);
    }

    public RewardsMixtureBranchResamplingHelper.BranchWeights getOperationCachedBranchWeightsForParameterIndex(
            final int parameterIndex) {
        final int branchNodeNumber = getNodeNumberForParameterIndex(parameterIndex);

        if (!operationCacheActive) {
            RewardMixturePerformanceStats.recordBranchWeightCacheBypass();
            refreshLikelihoodMessages();
            return computeBranchWeightsForNode(branchNodeNumber);
        }

        if (!likelihoodMessagesFresh) {
            RewardMixturePerformanceStats.recordLikelihoodMessageMiss();
            refreshLikelihoodMessages();
        }
        if (cachedBranchWeightsEpoch[parameterIndex] != operationCacheEpoch) {
            RewardMixturePerformanceStats.recordBranchWeightCacheMiss();
            cachedBranchWeightsByParameterIndex[parameterIndex] =
                    computeBranchWeightsForNode(branchNodeNumber);
            cachedBranchWeightsEpoch[parameterIndex] = operationCacheEpoch;
        } else {
            RewardMixturePerformanceStats.recordBranchWeightCacheHit();
        }
        return cachedBranchWeightsByParameterIndex[parameterIndex];
    }

    private void advanceOperationCacheEpoch() {
        operationCacheEpoch++;
        if (operationCacheEpoch == Integer.MAX_VALUE) {
            Arrays.fill(cachedBranchWeightsEpoch, 0);
            operationCacheEpoch = 1;
        }
    }

    public RewardsMixtureBranchResamplingHelper.BranchWeights computeBranchWeightsForNode(
            final int branchNodeNumber) {
        final long start = RewardMixturePerformanceStats.startTimer();
        try {
            return computeBranchWeightsForNodeUnprofiled(branchNodeNumber);
        } finally {
            RewardMixturePerformanceStats.recordBranchWeightComputation(
                    RewardMixturePerformanceStats.elapsed(start));
        }
    }

    private RewardsMixtureBranchResamplingHelper.BranchWeights computeBranchWeightsForNodeUnprofiled(
            final int branchNodeNumber) {
        long messageStart = RewardMixturePerformanceStats.startTimer();
        final int parameterIndex = getParameterIndexForNode(branchNodeNumber);
        final double preScale;
        final double postScale;

        if (hasCachedBranchMessages(parameterIndex)) {
            loadCachedBranchMessages(parameterIndex, prePartial, postPartial);
            preScale = getCachedPreScale(parameterIndex);
            postScale = getCachedPostScale(parameterIndex);
        } else {
            loadBranchPartials(branchNodeNumber, prePartial, postPartial);

            long start = RewardMixturePerformanceStats.startTimer();
            discreteDelegate.getPreOrderBranchScalesInto(branchNodeNumber, preScales);
            RewardMixturePerformanceStats.recordBranchPreScaleLoad(
                    RewardMixturePerformanceStats.elapsed(start));

            start = RewardMixturePerformanceStats.startTimer();
            discreteDelegate.getPostOrderBranchScalesInto(branchNodeNumber, postScales);
            RewardMixturePerformanceStats.recordBranchPostScaleLoad(
                    RewardMixturePerformanceStats.elapsed(start));

            preScale = preScales[0];
            postScale = postScales[0];
        }

        RewardMixturePerformanceStats.recordBranchMessageLoad(
                RewardMixturePerformanceStats.elapsed(messageStart));

        final NodeRef node = tree.getNode(branchNodeNumber);
        final double branchLength = tree.getBranchLength(node);

        long start = RewardMixturePerformanceStats.startTimer();
        for (int j = 0; j < nstates; j++) {
            final double logAtomicLocalFactor =
                    rewardsAwareBranchModel.getAtomicLogScaleForState(j, branchLength);
            final double dependentLogEvidence =
                    getDependentLogEvidence(branchNodeNumber, rewardsAwareBranchModel.getRewardRateRawForState(j));
            logAtomicWeights[j] =
                    RewardsMixtureBranchResamplingHelper.logAtomicWeight(
                            prePartial[j],
                            postPartial[j],
                            logAtomicLocalFactor,
                            preScale,
                            postScale
                    ) + dependentLogEvidence;
        }

        final double logAtomicTotalWeight =
                RewardsMixtureBranchResamplingHelper.logSum(logAtomicWeights, nstates);
        RewardMixturePerformanceStats.recordBranchAtomicWeightLoop(
                RewardMixturePerformanceStats.elapsed(start));

        final double rawContinuousReward =
                rewardsAwareBranchModel.getContinuousRewardRawForBranch(branchNodeNumber);

        final double logCtsWeight;
        if (branchLength == 0.0) {
            // A mid-branch CTMC jump is impossible on a zero-length branch, so
            // the continuous category has no support there -- exclude it
            // outright rather than relying on a degenerate (identity) transition
            // matrix, which would leave the branch selectable as "continuous"
            // and then crash the cts-reward HMC gradient (Sericola derivative
            // code requires time > 0). See the ctmc_bm4d_timeseries scenario (c)
            // ladder-tree diagnosis in the project log for the full trace.
            logCtsWeight = Double.NEGATIVE_INFINITY;
        } else if (isContinuousRewardOutsideOpenSupport(rawContinuousReward)) {
            logCtsWeight = Double.NEGATIVE_INFINITY;
        } else {
            start = RewardMixturePerformanceStats.startTimer();
            rewardsAwareBranchModel.computeTransitionMatrixCtsForBranchInto(
                    branchNodeNumber,
                    continuousMatrixScratch,
                    singleCtsRewardScratch,
                    singleCtsTimeScratch,
                    singleCtsMatrixScratch);
            RewardMixturePerformanceStats.recordBranchContinuousMatrixAccess(
                    RewardMixturePerformanceStats.elapsed(start));

            start = RewardMixturePerformanceStats.startTimer();
            logCtsWeight =
                    RewardsMixtureBranchResamplingHelper.logContinuousWeight(
                            prePartial,
                            continuousMatrixScratch,
                            postPartial,
                            nstates,
                            preScale,
                            postScale
                    ) + getDependentLogEvidence(branchNodeNumber, rawContinuousReward);
            RewardMixturePerformanceStats.recordBranchContinuousWeight(
                    RewardMixturePerformanceStats.elapsed(start));
        }

        return new RewardsMixtureBranchResamplingHelper.BranchWeights(
                Arrays.copyOf(logAtomicWeights, nstates),
                logAtomicTotalWeight,
                logCtsWeight
        );
    }

    public double computeContinuousLogWeightForNode(final int branchNodeNumber,
                                                   final double rawContinuousReward) {
        final long start = RewardMixturePerformanceStats.startTimer();
        try {
            return computeContinuousLogWeightForNodeUnprofiled(branchNodeNumber, rawContinuousReward);
        } finally {
            RewardMixturePerformanceStats.recordBranchWeightComputation(
                    RewardMixturePerformanceStats.elapsed(start));
        }
    }

    private double computeContinuousLogWeightForNodeUnprofiled(final int branchNodeNumber,
                                                              final double rawContinuousReward) {
        long messageStart = RewardMixturePerformanceStats.startTimer();
        final int parameterIndex = getParameterIndexForNode(branchNodeNumber);
        final double preScale;
        final double postScale;

        if (hasCachedBranchMessages(parameterIndex)) {
            loadCachedBranchMessages(parameterIndex, prePartial, postPartial);
            preScale = getCachedPreScale(parameterIndex);
            postScale = getCachedPostScale(parameterIndex);
        } else {
            loadBranchPartials(branchNodeNumber, prePartial, postPartial);

            long start = RewardMixturePerformanceStats.startTimer();
            discreteDelegate.getPreOrderBranchScalesInto(branchNodeNumber, preScales);
            RewardMixturePerformanceStats.recordBranchPreScaleLoad(
                    RewardMixturePerformanceStats.elapsed(start));

            start = RewardMixturePerformanceStats.startTimer();
            discreteDelegate.getPostOrderBranchScalesInto(branchNodeNumber, postScales);
            RewardMixturePerformanceStats.recordBranchPostScaleLoad(
                    RewardMixturePerformanceStats.elapsed(start));

            preScale = preScales[0];
            postScale = postScales[0];
        }

        RewardMixturePerformanceStats.recordBranchMessageLoad(
                RewardMixturePerformanceStats.elapsed(messageStart));

        final NodeRef node = tree.getNode(branchNodeNumber);
        final double branchLength = tree.getBranchLength(node);

        if (branchLength == 0.0 || isContinuousRewardOutsideOpenSupport(rawContinuousReward)) {
            return Double.NEGATIVE_INFINITY;
        }

        long start = RewardMixturePerformanceStats.startTimer();
        rewardsAwareBranchModel.computeTransitionMatrixCtsForBranchInto(
                branchNodeNumber,
                rawContinuousReward,
                continuousMatrixScratch,
                singleCtsRewardScratch,
                singleCtsTimeScratch,
                singleCtsMatrixScratch);
        RewardMixturePerformanceStats.recordBranchContinuousMatrixAccess(
                RewardMixturePerformanceStats.elapsed(start));

        start = RewardMixturePerformanceStats.startTimer();
        final double logCtsWeight =
                RewardsMixtureBranchResamplingHelper.logContinuousWeight(
                        prePartial,
                        continuousMatrixScratch,
                        postPartial,
                        nstates,
                        preScale,
                        postScale
                ) + getDependentLogEvidence(branchNodeNumber, rawContinuousReward);
        RewardMixturePerformanceStats.recordBranchContinuousWeight(
                RewardMixturePerformanceStats.elapsed(start));

        return logCtsWeight;
    }

    public double getLogWeightForCategory(final RewardsMixtureBranchResamplingHelper.BranchWeights weights,
                                          final int category) {
        if (category == 0) {
            return weights.logCtsWeight;
        }
        final int atomState = category - 1;
        if (atomState < 0 || atomState >= nstates) {
            throw new IllegalArgumentException("Reward mixture category out of range: " + category);
        }
        return weights.logAtomicWeights[atomState];
    }

    public double computeLogWeightForCategoryAtValue(final int parameterIndex,
                                                     final int category,
                                                     final double rawContinuousReward) {
        if (category == 0) {
            return computeContinuousLogWeightForParameterIndex(parameterIndex, rawContinuousReward);
        }

        // Atomic categories do not depend on the candidate CTS reward value.
        // Their CTS contribution is carried by the optional atomic pseudo-prior.
        final int atomState = category - 1;
        if (atomState < 0 || atomState >= nstates) {
            throw new IllegalArgumentException("Reward mixture category out of range: " + category);
        }

        final RewardsMixtureBranchResamplingHelper.BranchWeights weights =
                getOperationCachedBranchWeightsForParameterIndex(parameterIndex);
        return weights.logAtomicWeights[atomState];
    }

    private void initializeBranchMappings() {
        Arrays.fill(nodeNumberByParameterIndex, -1);

        int observedBranches = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            final NodeRef node = tree.getNode(i);
            if (tree.isRoot(node)) {
                continue;
            }

            final int nodeNumber = node.getNumber();
            final int parameterIndex = rewardsAwareBranchModel.getParameterIndexForNode(nodeNumber);
            if (parameterIndex < 0 || parameterIndex >= branchCount) {
                throw new IllegalArgumentException(
                        "Invalid branch parameter index " + parameterIndex + " for node " + nodeNumber +
                                "; branch parameter dimension is " + branchCount);
            }
            if (nodeNumberByParameterIndex[parameterIndex] != -1) {
                throw new IllegalArgumentException(
                        "Multiple non-root nodes map to branch parameter index " + parameterIndex);
            }

            nodeNumberByParameterIndex[parameterIndex] = nodeNumber;
            observedBranches++;
        }

        if (observedBranches != branchCount) {
            throw new IllegalArgumentException(
                    "Observed " + observedBranches + " non-root branches, but branchCount is " + branchCount);
        }
        for (int i = 0; i < branchCount; i++) {
            if (nodeNumberByParameterIndex[i] < 0) {
                throw new IllegalArgumentException("No branch node maps to parameter index " + i);
            }
        }
    }

    private void loadBranchPartials(final int nodeNum,
                                    final double[] prePartialOut,
                                    final double[] postPartialOut) {
        long start = RewardMixturePerformanceStats.startTimer();
        Arrays.fill(prePartialOut, 0.0);
        discreteDelegate.getPreOrderAtBranchStartInto(nodeNum, prePartialOut);
        RewardMixturePerformanceStats.recordBranchPrePartialLoad(
                RewardMixturePerformanceStats.elapsed(start));

        start = RewardMixturePerformanceStats.startTimer();
        Arrays.fill(postPartialOut, 0.0);
        discreteDelegate.getPostOrderAtBranchEndInto(nodeNum, postPartialOut);
        RewardMixturePerformanceStats.recordBranchPostPartialLoad(
                RewardMixturePerformanceStats.elapsed(start));
    }

    private void fillBranchMessageCache() {
        final long start = RewardMixturePerformanceStats.startTimer();
        for (int parameterIndex = 0; parameterIndex < branchCount; parameterIndex++) {
            final int branchNodeNumber = nodeNumberByParameterIndex[parameterIndex];

            Arrays.fill(cachedPrePartialByParameterIndex[parameterIndex], 0.0);
            discreteDelegate.getPreOrderAtBranchStartInto(
                    branchNodeNumber, cachedPrePartialByParameterIndex[parameterIndex]);

            Arrays.fill(cachedPostPartialByParameterIndex[parameterIndex], 0.0);
            discreteDelegate.getPostOrderAtBranchEndInto(
                    branchNodeNumber, cachedPostPartialByParameterIndex[parameterIndex]);

            discreteDelegate.getPreOrderBranchScalesInto(branchNodeNumber, preScales);
            discreteDelegate.getPostOrderBranchScalesInto(branchNodeNumber, postScales);
            cachedPreScaleByParameterIndex[parameterIndex] = preScales[0];
            cachedPostScaleByParameterIndex[parameterIndex] = postScales[0];
        }
        RewardMixturePerformanceStats.recordBranchMessageCacheFill(
                RewardMixturePerformanceStats.elapsed(start));
    }

    private void fillSingleBranchMessageCache(final int parameterIndex) {
        final long start = RewardMixturePerformanceStats.startTimer();
        final int branchNodeNumber = nodeNumberByParameterIndex[parameterIndex];

        Arrays.fill(cachedSinglePrePartial, 0.0);
        discreteDelegate.getPreOrderAtBranchStartInto(branchNodeNumber, cachedSinglePrePartial);

        Arrays.fill(cachedSinglePostPartial, 0.0);
        discreteDelegate.getPostOrderAtBranchEndInto(branchNodeNumber, cachedSinglePostPartial);

        discreteDelegate.getPreOrderBranchScalesInto(branchNodeNumber, preScales);
        discreteDelegate.getPostOrderBranchScalesInto(branchNodeNumber, postScales);
        cachedSinglePreScale = preScales[0];
        cachedSinglePostScale = postScales[0];

        RewardMixturePerformanceStats.recordBranchMessageCacheFill(
                RewardMixturePerformanceStats.elapsed(start));
    }

    private void loadCachedBranchMessages(final int parameterIndex,
                                          final double[] prePartialOut,
                                          final double[] postPartialOut) {
        final long start = RewardMixturePerformanceStats.startTimer();
        final double[] cachedPre = cacheAllBranchMessages
                ? cachedPrePartialByParameterIndex[parameterIndex]
                : cachedSinglePrePartial;
        final double[] cachedPost = cacheAllBranchMessages
                ? cachedPostPartialByParameterIndex[parameterIndex]
                : cachedSinglePostPartial;

        System.arraycopy(cachedPre, 0,
                prePartialOut, 0, nstates);
        RewardMixturePerformanceStats.recordBranchPrePartialLoad(
                RewardMixturePerformanceStats.elapsed(start));

        final long postStart = RewardMixturePerformanceStats.startTimer();
        System.arraycopy(cachedPost, 0,
                postPartialOut, 0, nstates);
        RewardMixturePerformanceStats.recordBranchPostPartialLoad(
                RewardMixturePerformanceStats.elapsed(postStart));
    }

    private boolean hasCachedBranchMessages(final int parameterIndex) {
        return operationCacheActive && likelihoodMessagesFresh &&
                (cacheAllBranchMessages ||
                        (cacheSingleBranchMessages && parameterIndex == singleBranchCacheParameterIndex));
    }

    private double getCachedPreScale(final int parameterIndex) {
        return cacheAllBranchMessages
                ? cachedPreScaleByParameterIndex[parameterIndex]
                : cachedSinglePreScale;
    }

    private double getCachedPostScale(final int parameterIndex) {
        return cacheAllBranchMessages
                ? cachedPostScaleByParameterIndex[parameterIndex]
                : cachedSinglePostScale;
    }

    private double getDependentLogEvidence(final int branchNodeNumber, final double rawReward) {
        final long start = RewardMixturePerformanceStats.startTimer();
        double logEvidence = 0.0;
        try {
            for (RewardDependentEdgeEvidenceProvider provider : dependentEvidenceProviders) {
                final double contribution = provider.logEvidence(branchNodeNumber, rawReward);
                if (!Double.isFinite(contribution)) {
                    return Double.NEGATIVE_INFINITY;
                }
                logEvidence += contribution;
            }
            return logEvidence;
        } finally {
            RewardMixturePerformanceStats.recordDependentEvidence(
                    RewardMixturePerformanceStats.elapsed(start));
        }
    }

    private boolean isContinuousRewardOutsideOpenSupport(final double rawReward) {
        double minReward = Double.POSITIVE_INFINITY;
        double maxReward = Double.NEGATIVE_INFINITY;
        for (int state = 0; state < nstates; state++) {
            final double reward = rewardsAwareBranchModel.getRewardRateRawForState(state);
            minReward = Math.min(minReward, reward);
            maxReward = Math.max(maxReward, reward);
        }
        return RewardsMixtureBranchResamplingHelper.isContinuousRewardOutsideOpenSupport(
                rawReward, minReward, maxReward);
    }

    private static DiscreteDataLikelihoodDelegate requireDiscreteDelegate(final TreeDataLikelihood treeDataLikelihood) {
        if (treeDataLikelihood == null) {
            throw new IllegalArgumentException("treeDataLikelihood must be non-null");
        }
        final DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (!(delegate instanceof DiscreteDataLikelihoodDelegate)) {
            throw new IllegalArgumentException(
                    "Reward-mixture branch weights require TreeDataLikelihood to use DiscreteDataLikelihoodDelegate");
        }
        return (DiscreteDataLikelihoodDelegate) delegate;
    }

    private static RewardDependentEdgeEvidenceProvider[] createDependentEvidenceProviders(
            final RewardsAwareBranchModel rewardsAwareBranchModel,
            final TreeDataLikelihood[] dependentTreeDataLikelihoods,
            final TreeDataLikelihood[] dependentContinuousTreeDataLikelihoods) {
        if (rewardsAwareBranchModel == null) {
            throw new IllegalArgumentException("rewardsAwareBranchModel must be non-null");
        }

        final Tree tree = rewardsAwareBranchModel.getTree();
        final TreeDataLikelihood[] ctmcLikelihoods = dependentTreeDataLikelihoods == null
                ? new TreeDataLikelihood[0]
                : Arrays.copyOf(dependentTreeDataLikelihoods, dependentTreeDataLikelihoods.length);
        final TreeDataLikelihood[] continuousLikelihoods = dependentContinuousTreeDataLikelihoods == null
                ? new TreeDataLikelihood[0]
                : Arrays.copyOf(dependentContinuousTreeDataLikelihoods, dependentContinuousTreeDataLikelihoods.length);

        final RewardDependentEdgeEvidenceProvider[] providers =
                new RewardDependentEdgeEvidenceProvider[ctmcLikelihoods.length + continuousLikelihoods.length];

        for (int i = 0; i < ctmcLikelihoods.length; i++) {
            final TreeDataLikelihood dependent = ctmcLikelihoods[i];
            validateDependentTree(dependent, tree, "dependentTreeDataLikelihoods", i);
            providers[i] = new BeagleRewardDependentCtmcEdgeEvidenceProvider(dependent);
        }
        for (int i = 0; i < continuousLikelihoods.length; i++) {
            final TreeDataLikelihood dependent = continuousLikelihoods[i];
            validateDependentTree(dependent, tree, "dependentContinuousTreeDataLikelihoods", i);
            providers[ctmcLikelihoods.length + i] =
                    new BranchLocalContinuousRewardDependentEdgeEvidenceProvider(dependent);
        }

        return providers;
    }

    private static void validateDependentTree(final TreeDataLikelihood dependent,
                                              final Tree tree,
                                              final String label,
                                              final int index) {
        if (dependent == null) {
            throw new IllegalArgumentException(label + " contains null at index " + index);
        }
        if (dependent.getTree().getNodeCount() != tree.getNodeCount()) {
            throw new IllegalArgumentException(
                    label + " at index " + index +
                            " must use a tree with the same node count as RewardsAwareBranchModel.");
        }
    }
}
