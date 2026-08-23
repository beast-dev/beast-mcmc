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
    private final RewardsMixtureBranchResamplingHelper.BranchWeights[] cachedBranchWeightsByParameterIndex;
    private final int[] cachedBranchWeightsEpoch;

    private boolean operationCacheActive = false;
    private boolean likelihoodMessagesFresh = false;
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
        rewardsAwareBranchModel.refreshCategoryDecoderEmbedding();
        clearOperationCache();
    }

    public void beginOperationCache() {
        operationCacheActive = true;
        clearOperationCache();
    }

    public void clearOperationCache() {
        likelihoodMessagesFresh = false;
        advanceOperationCacheEpoch();
    }

    public void refreshLikelihoodMessages() {
        discreteDelegate.updatePostOrdersFromTreeDataLikelihood(treeDataLikelihood);
        discreteDelegate.ensurePreOrderComputed();
        for (RewardDependentEdgeEvidenceProvider provider : dependentEvidenceProviders) {
            provider.prepare();
        }
        if (operationCacheActive) {
            likelihoodMessagesFresh = true;
            advanceOperationCacheEpoch();
        }
    }

    public RewardsMixtureBranchResamplingHelper.BranchWeights computeBranchWeightsForParameterIndex(
            final int parameterIndex) {
        return computeBranchWeightsForNode(getNodeNumberForParameterIndex(parameterIndex));
    }

    public RewardsMixtureBranchResamplingHelper.BranchWeights getOperationCachedBranchWeightsForParameterIndex(
            final int parameterIndex) {
        final int branchNodeNumber = getNodeNumberForParameterIndex(parameterIndex);

        if (!operationCacheActive) {
            refreshLikelihoodMessages();
            return computeBranchWeightsForNode(branchNodeNumber);
        }

        if (!likelihoodMessagesFresh) {
            refreshLikelihoodMessages();
        }
        if (cachedBranchWeightsEpoch[parameterIndex] != operationCacheEpoch) {
            cachedBranchWeightsByParameterIndex[parameterIndex] =
                    computeBranchWeightsForNode(branchNodeNumber);
            cachedBranchWeightsEpoch[parameterIndex] = operationCacheEpoch;
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
        loadBranchPartials(branchNodeNumber, prePartial, postPartial);

        discreteDelegate.getPreOrderBranchScalesInto(branchNodeNumber, preScales);
        discreteDelegate.getPostOrderBranchScalesInto(branchNodeNumber, postScales);

        final double preScale = preScales[0];
        final double postScale = postScales[0];

        final NodeRef node = tree.getNode(branchNodeNumber);
        final double branchLength = tree.getBranchLength(node);

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

        final double rawContinuousReward =
                rewardsAwareBranchModel.getContinuousRewardRawForBranch(branchNodeNumber);

        final double logCtsWeight;
        if (isContinuousRewardOutsideOpenSupport(rawContinuousReward)) {
            logCtsWeight = Double.NEGATIVE_INFINITY;
        } else {
            final double[] continuousMatrix = rewardsAwareBranchModel.getTransitionMatrixCts(branchNodeNumber);
            logCtsWeight =
                    RewardsMixtureBranchResamplingHelper.logContinuousWeight(
                            prePartial,
                            continuousMatrix,
                            postPartial,
                            nstates,
                            preScale,
                            postScale
                    ) + getDependentLogEvidence(branchNodeNumber, rawContinuousReward);
        }

        return new RewardsMixtureBranchResamplingHelper.BranchWeights(
                Arrays.copyOf(logAtomicWeights, nstates),
                logAtomicTotalWeight,
                logCtsWeight
        );
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
        Arrays.fill(prePartialOut, 0.0);
        Arrays.fill(postPartialOut, 0.0);

        discreteDelegate.getPreOrderAtBranchStartInto(nodeNum, prePartialOut);
        discreteDelegate.getPostOrderAtBranchEndInto(nodeNum, postPartialOut);
    }

    private double getDependentLogEvidence(final int branchNodeNumber, final double rawReward) {
        double logEvidence = 0.0;
        for (RewardDependentEdgeEvidenceProvider provider : dependentEvidenceProviders) {
            final double contribution = provider.logEvidence(branchNodeNumber, rawReward);
            if (!Double.isFinite(contribution)) {
                return Double.NEGATIVE_INFINITY;
            }
            logEvidence += contribution;
        }
        return logEvidence;
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
