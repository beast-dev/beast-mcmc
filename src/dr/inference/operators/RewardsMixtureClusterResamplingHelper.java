package dr.inference.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.treedatalikelihood.DiscreteDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

import java.util.Arrays;

/**
 * Helper for blocked cluster proposals on the discrete mixture state:
 * - indicatorZ[b] in {0,1}
 * - atomIndex[b] in {0,...,nstates-1} when indicatorZ[b] = 1
 */
public final class RewardsMixtureClusterResamplingHelper {

    /**
     * Callback used to evaluate the exact current log target under the present parameter state.
     */
    public interface ExactLogTargetEvaluator {
        double computeCurrentLogTarget(int[] clusterBranchNodeNumbers, int clusterCount);
    }

    /**
     * Callback used after quiet state edits to refresh post-order/pre-order
     * messages before evaluating proposal probabilities under that state.
     */
    public interface LikelihoodStateRefresher {
        void refreshCurrentState();
    }

    /**
     * Result of a blocked proposal. The proposed new state remains applied after construction.
     */
    public static final class ClusterProposal {
        public final int[] branchNodeNumbers;
        public final int clusterCount;
        public final int[] parameterIndices;
        public final int[] oldIndicators;
        public final int[] oldAtoms;
        public final int[] newIndicators;
        public final int[] newAtoms;
        public final double logTargetOld;
        public final double logTargetNew;
        public final double logQForward;
        public final double logQReverse;

        private final Parameter indicatorZ;
        private final Parameter atomIndex;

        private ClusterProposal(
                final int[] branchNodeNumbers,
                final int clusterCount,
                final int[] parameterIndices,
                final int[] oldIndicators,
                final int[] oldAtoms,
                final int[] newIndicators,
                final int[] newAtoms,
                final double logTargetOld,
                final double logTargetNew,
                final double logQForward,
                final double logQReverse,
                final Parameter indicatorZ,
                final Parameter atomIndex
        ) {
            this.branchNodeNumbers = branchNodeNumbers;
            this.clusterCount = clusterCount;
            this.parameterIndices = parameterIndices;
            this.oldIndicators = oldIndicators;
            this.oldAtoms = oldAtoms;
            this.newIndicators = newIndicators;
            this.newAtoms = newAtoms;
            this.logTargetOld = logTargetOld;
            this.logTargetNew = logTargetNew;
            this.logQForward = logQForward;
            this.logQReverse = logQReverse;
            this.indicatorZ = indicatorZ;
            this.atomIndex = atomIndex;
        }

        public double getLogAcceptanceRatio() {
            return (logTargetNew - logTargetOld) + (logQReverse - logQForward);
        }

        public boolean acceptWithCurrentProposalApplied() {
            return Math.log(MathUtils.nextDouble()) < getLogAcceptanceRatio();
        }

        public void restoreOldState() {
            for (int i = 0; i < clusterCount; i++) {
                final int branchIndex = parameterIndices[i];
                indicatorZ.setParameterValueQuietly(branchIndex, oldIndicators[i]);
                atomIndex.setParameterValueQuietly(branchIndex, oldAtoms[i]);
            }
        }

        public void fireParameterEvents() {
            atomIndex.fireParameterChangedEvent();
            indicatorZ.fireParameterChangedEvent();
        }
    }

    private final Parameter indicatorZ;
    private final Parameter atomIndex;
    private final RewardsAwareBranchModel rewardsAwareBranchModel;
    private final TreeDataLikelihood treeDataLikelihood;
    private final DiscreteDataLikelihoodDelegate discreteDelegate;
    private final RewardDependentCtmcEdgeEvidenceProvider[] dependentEvidenceProviders;
    private final LikelihoodStateRefresher likelihoodStateRefresher;
    private final Tree tree;
    private final int nstates;
    private final int maxClusterSize;
    private final double borderBias;
    private final ExactLogTargetEvaluator logTargetEvaluator;

    private final double[] prePartial;
    private final double[] postPartial;
    private final double[] logAtomicWeights;
    private final double[] preScales;
    private final double[] postScales;

    private final int[] borderBuffer;
    private int borderCount = 0;

    private final int[] queueBuffer;
    private final boolean[] visited;

    public RewardsMixtureClusterResamplingHelper(
            final Parameter indicatorZ,
            final Parameter atomIndex,
            final RewardsAwareBranchModel rewardsAwareBranchModel,
            final TreeDataLikelihood treeDataLikelihood,
            final DiscreteDataLikelihoodDelegate discreteDelegate,
            final int maxClusterSize,
            final double borderBias,
            final ExactLogTargetEvaluator logTargetEvaluator
    ) {
        this(
                indicatorZ,
                atomIndex,
                rewardsAwareBranchModel,
                treeDataLikelihood,
                discreteDelegate,
                null,
                maxClusterSize,
                borderBias,
                logTargetEvaluator
        );
    }

    public RewardsMixtureClusterResamplingHelper(
            final Parameter indicatorZ,
            final Parameter atomIndex,
            final RewardsAwareBranchModel rewardsAwareBranchModel,
            final TreeDataLikelihood treeDataLikelihood,
            final DiscreteDataLikelihoodDelegate discreteDelegate,
            final RewardDependentCtmcEdgeEvidenceProvider[] dependentEvidenceProviders,
            final int maxClusterSize,
            final double borderBias,
            final ExactLogTargetEvaluator logTargetEvaluator
    ) {
        this(
                indicatorZ,
                atomIndex,
                rewardsAwareBranchModel,
                treeDataLikelihood,
                discreteDelegate,
                dependentEvidenceProviders,
                maxClusterSize,
                borderBias,
                logTargetEvaluator,
                null
        );
    }

    public RewardsMixtureClusterResamplingHelper(
            final Parameter indicatorZ,
            final Parameter atomIndex,
            final RewardsAwareBranchModel rewardsAwareBranchModel,
            final TreeDataLikelihood treeDataLikelihood,
            final DiscreteDataLikelihoodDelegate discreteDelegate,
            final RewardDependentCtmcEdgeEvidenceProvider[] dependentEvidenceProviders,
            final int maxClusterSize,
            final double borderBias,
            final ExactLogTargetEvaluator logTargetEvaluator,
            final LikelihoodStateRefresher likelihoodStateRefresher
    ) {
        if (indicatorZ == null) {
            throw new IllegalArgumentException("indicatorZ cannot be null.");
        }
        if (atomIndex == null) {
            throw new IllegalArgumentException("atomIndex cannot be null.");
        }
        if (rewardsAwareBranchModel == null) {
            throw new IllegalArgumentException("rewardsAwareBranchModel cannot be null.");
        }
        if (treeDataLikelihood == null) {
            throw new IllegalArgumentException("treeDataLikelihood cannot be null.");
        }
        if (discreteDelegate == null) {
            throw new IllegalArgumentException("discreteDelegate cannot be null.");
        }
        if (logTargetEvaluator == null) {
            throw new IllegalArgumentException("logTargetEvaluator cannot be null.");
        }
        if (indicatorZ.getDimension() != atomIndex.getDimension()) {
            throw new IllegalArgumentException("indicatorZ and atomIndex must have same dimension.");
        }

        this.indicatorZ = indicatorZ;
        this.atomIndex = atomIndex;
        this.rewardsAwareBranchModel = rewardsAwareBranchModel;
        this.treeDataLikelihood = treeDataLikelihood;
        this.discreteDelegate = discreteDelegate;
        this.dependentEvidenceProviders = dependentEvidenceProviders == null
                ? new RewardDependentCtmcEdgeEvidenceProvider[0]
                : Arrays.copyOf(dependentEvidenceProviders, dependentEvidenceProviders.length);
        this.likelihoodStateRefresher = likelihoodStateRefresher;
        this.tree = rewardsAwareBranchModel.getTree();
        this.nstates = rewardsAwareBranchModel.getStateCount();
        this.maxClusterSize = Math.max(1, maxClusterSize);
        this.borderBias = Math.max(0.0, Math.min(1.0, borderBias));
        this.logTargetEvaluator = logTargetEvaluator;

        this.prePartial = new double[nstates];
        this.postPartial = new double[nstates];
        this.logAtomicWeights = new double[nstates];
        this.preScales = new double[discreteDelegate.getPatternCount()];
        this.postScales = new double[discreteDelegate.getPatternCount()];

        this.borderBuffer = new int[tree.getNodeCount()];
        this.queueBuffer = new int[tree.getNodeCount()];
        this.visited = new boolean[tree.getNodeCount()];
    }

    public ClusterProposal proposeClusterMove(final int seedBranchNodeNumber) {

        final int[] clusterBranchNodeNumbers = new int[maxClusterSize];
        final int clusterCount = buildConnectedCluster(seedBranchNodeNumber, clusterBranchNodeNumbers);

        final int[] branchNodeNumbers = Arrays.copyOf(clusterBranchNodeNumbers, clusterCount);
        final int[] parameterIndices = new int[clusterCount];
        final int[] oldIndicators = new int[clusterCount];
        final int[] oldAtoms = new int[clusterCount];
        final int[] newIndicators = new int[clusterCount];
        final int[] newAtoms = new int[clusterCount];

        for (int i = 0; i < clusterCount; i++) {
            final int branchNodeNumber = branchNodeNumbers[i];
            final int branchIndex = rewardsAwareBranchModel.getParameterIndexForNode(branchNodeNumber);
            parameterIndices[i] = branchIndex;

            oldIndicators[i] = (int) Math.round(indicatorZ.getParameterValue(branchIndex));
            oldAtoms[i] = (int) Math.round(atomIndex.getParameterValue(branchIndex));
            newIndicators[i] = oldIndicators[i];
            newAtoms[i] = oldAtoms[i];
        }

        refreshLikelihoodState();

        final double logTargetOld =
                logTargetEvaluator.computeCurrentLogTarget(branchNodeNumbers, clusterCount);

        final double logQForward =
                proposeForward(branchNodeNumbers, clusterCount, newIndicators, newAtoms);

        applyState(parameterIndices, clusterCount, newIndicators, newAtoms);
        fireParameterEvents();
        refreshLikelihoodState();

        final double logTargetNew =
                logTargetEvaluator.computeCurrentLogTarget(branchNodeNumbers, clusterCount);

        final double logQReverse =
                computeReverseProposalProbability(branchNodeNumbers, clusterCount, oldIndicators, oldAtoms);

        return new ClusterProposal(
                branchNodeNumbers,
                clusterCount,
                parameterIndices,
                oldIndicators,
                oldAtoms,
                newIndicators,
                newAtoms,
                logTargetOld,
                logTargetNew,
                logQForward,
                logQReverse,
                indicatorZ,
                atomIndex
        );
    }

    public boolean doOneClusterMoveAndFireEvents() {
        final int seedBranchNodeNumber = chooseSeedBranch();
        final ClusterProposal proposal = proposeClusterMove(seedBranchNodeNumber);

        final boolean accept = proposal.acceptWithCurrentProposalApplied();
        if (!accept) {
            proposal.restoreOldState();
        }
        proposal.fireParameterEvents();
        return accept;
    }

    private void refreshLikelihoodState() {
        if (likelihoodStateRefresher != null) {
            likelihoodStateRefresher.refreshCurrentState();
        } else {
            treeDataLikelihood.makeDirty();
            discreteDelegate.updatePostOrdersFromTreeDataLikelihood(treeDataLikelihood);
            discreteDelegate.ensurePreOrderComputed();
        }
    }

    private void fireParameterEvents() {
        atomIndex.fireParameterChangedEvent();
        indicatorZ.fireParameterChangedEvent();
    }

    public int chooseSeedBranch() {
        collectBorderBranches();

        if (borderCount > 0 && MathUtils.nextDouble() < borderBias) {
            return borderBuffer[MathUtils.nextInt(borderCount)];
        }

        int branchNodeNumber;
        do {
            branchNodeNumber = MathUtils.nextInt(tree.getNodeCount());
        } while (tree.isRoot(tree.getNode(branchNodeNumber)));

        return branchNodeNumber;
    }

    private void collectBorderBranches() {
        borderCount = 0;

        for (int nodeNum = 0; nodeNum < tree.getNodeCount(); nodeNum++) {
            final NodeRef node = tree.getNode(nodeNum);

            if (tree.isRoot(node)) {
                continue;
            }

            final int branchIndex = rewardsAwareBranchModel.getParameterIndexForNode(nodeNum);
            final int branchIndicator = (int) Math.round(indicatorZ.getParameterValue(branchIndex));

            boolean isBorder = false;

            final NodeRef parent = tree.getParent(node);
            if (parent != null && !tree.isRoot(parent)) {
                final int parentNodeNum = parent.getNumber();
                final int parentBranchIndex = rewardsAwareBranchModel.getParameterIndexForNode(parentNodeNum);
                final int parentIndicator = (int) Math.round(indicatorZ.getParameterValue(parentBranchIndex));
                if (parentIndicator != branchIndicator) {
                    isBorder = true;
                }
            }

            if (!isBorder) {
                final int childCount = tree.getChildCount(node);
                for (int i = 0; i < childCount; i++) {
                    final NodeRef child = tree.getChild(node, i);
                    final int childNodeNum = child.getNumber();
                    final int childBranchIndex = rewardsAwareBranchModel.getParameterIndexForNode(childNodeNum);
                    final int childIndicator = (int) Math.round(indicatorZ.getParameterValue(childBranchIndex));
                    if (childIndicator != branchIndicator) {
                        isBorder = true;
                        break;
                    }
                }
            }

            if (isBorder) {
                borderBuffer[borderCount++] = nodeNum;
            }
        }
    }

    private int buildConnectedCluster(final int seedBranchNodeNumber, final int[] out) {
        Arrays.fill(visited, false);

        int queueHead = 0;
        int queueTail = 0;
        int count = 0;

        queueBuffer[queueTail++] = seedBranchNodeNumber;
        visited[seedBranchNodeNumber] = true;

        while (queueHead < queueTail && count < maxClusterSize) {
            final int branchNodeNumber = queueBuffer[queueHead++];
            out[count++] = branchNodeNumber;

            final NodeRef node = tree.getNode(branchNodeNumber);

            final NodeRef parent = tree.getParent(node);
            if (parent != null && !tree.isRoot(parent)) {
                final int parentBranchNodeNumber = parent.getNumber();
                if (!visited[parentBranchNodeNumber]) {
                    visited[parentBranchNodeNumber] = true;
                    queueBuffer[queueTail++] = parentBranchNodeNumber;
                }
            }

            final int childCount = tree.getChildCount(node);
            for (int i = 0; i < childCount; i++) {
                final NodeRef child = tree.getChild(node, i);
                final int childBranchNodeNumber = child.getNumber();
                if (!visited[childBranchNodeNumber]) {
                    visited[childBranchNodeNumber] = true;
                    queueBuffer[queueTail++] = childBranchNodeNumber;
                }
            }
        }

        return count;
    }

    private double proposeForward(
            final int[] branchNodeNumbers,
            final int clusterCount,
            final int[] newIndicators,
            final int[] newAtoms
    ) {
        double logQ = 0.0;

        for (int i = 0; i < clusterCount; i++) {
            final int branchNodeNumber = branchNodeNumbers[i];

            final RewardsMixtureBranchResamplingHelper.BranchWeights weights =
                    computeBranchWeights(branchNodeNumber);

            final double logAtomic = weights.logAtomicTotalWeight;
            final double logCts = weights.logCtsWeight;
            final double logNorm = logSumPair(logAtomic, logCts);

            final int newIndicator =
                    RewardsMixtureBranchResamplingHelper.sampleIndicatorFromLogs(logAtomic, logCts);

            newIndicators[i] = newIndicator;

            if (newIndicator == 1) {
                final int newAtom =
                        RewardsMixtureBranchResamplingHelper.sampleAtomFromLogs(weights.logAtomicWeights, nstates);
                newAtoms[i] = newAtom;

                logQ += logAtomic - logNorm;
                logQ += weights.logAtomicWeights[newAtom] - logAtomic;
            } else {
                logQ += logCts - logNorm;
                final int branchIndex = rewardsAwareBranchModel.getParameterIndexForNode(branchNodeNumber);
                newAtoms[i] = (int) Math.round(atomIndex.getParameterValue(branchIndex));
            }
        }

        return logQ;
    }

    private double computeReverseProposalProbability(
            final int[] branchNodeNumbers,
            final int clusterCount,
            final int[] oldIndicators,
            final int[] oldAtoms
    ) {
        double logQ = 0.0;

        for (int i = 0; i < clusterCount; i++) {
            final int branchNodeNumber = branchNodeNumbers[i];

            final RewardsMixtureBranchResamplingHelper.BranchWeights weights =
                    computeBranchWeights(branchNodeNumber);

            final double logAtomic = weights.logAtomicTotalWeight;
            final double logCts = weights.logCtsWeight;
            final double logNorm = logSumPair(logAtomic, logCts);

            if (oldIndicators[i] == 1) {
                final int oldAtom = oldAtoms[i];
                logQ += logAtomic - logNorm;
                logQ += weights.logAtomicWeights[oldAtom] - logAtomic;
            } else {
                logQ += logCts - logNorm;
            }
        }

        return logQ;
    }

    private void applyState(
            final int[] parameterIndices,
            final int clusterCount,
            final int[] indicators,
            final int[] atoms
    ) {
        for (int i = 0; i < clusterCount; i++) {
            final int branchIndex = parameterIndices[i];
            indicatorZ.setParameterValueQuietly(branchIndex, indicators[i]);
            if (indicators[i] == 1) {
                atomIndex.setParameterValueQuietly(branchIndex, atoms[i]);
            }
        }
    }

    private RewardsMixtureBranchResamplingHelper.BranchWeights computeBranchWeights(final int branchNodeNumber) {

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

        final double[] continuousMatrix = rewardsAwareBranchModel.getTransitionMatrixCts(branchNodeNumber);

        final double logCtsWeight =
                RewardsMixtureBranchResamplingHelper.logContinuousWeight(
                        prePartial,
                        continuousMatrix,
                        postPartial,
                        nstates,
                        preScale,
                        postScale
                ) + getDependentLogEvidence(
                        branchNodeNumber,
                        rewardsAwareBranchModel.getContinuousRewardRawForBranch(branchNodeNumber)
                );

        return new RewardsMixtureBranchResamplingHelper.BranchWeights(
                Arrays.copyOf(logAtomicWeights, nstates),
                logAtomicTotalWeight,
                logCtsWeight
        );
    }

    private double getDependentLogEvidence(final int branchNodeNumber, final double rawReward) {
        double logEvidence = 0.0;
        for (RewardDependentCtmcEdgeEvidenceProvider provider : dependentEvidenceProviders) {
            final double contribution = provider.logEvidence(branchNodeNumber, rawReward);
            if (!Double.isFinite(contribution)) {
                return Double.NEGATIVE_INFINITY;
            }
            logEvidence += contribution;
        }
        return logEvidence;
    }

    private void loadBranchPartials(final int nodeNum,
                                    final double[] prePartialOut,
                                    final double[] postPartialOut) {
        Arrays.fill(prePartialOut, 0.0);
        Arrays.fill(postPartialOut, 0.0);

        discreteDelegate.getPreOrderAtBranchStartInto(nodeNum, prePartialOut);
        discreteDelegate.getPostOrderAtBranchEndInto(nodeNum, postPartialOut);
    }

    private static double logSumPair(final double a, final double b) {
        if (Double.isInfinite(a) && a < 0.0) {
            return b;
        }
        if (Double.isInfinite(b) && b < 0.0) {
            return a;
        }
        final double m = Math.max(a, b);
        return m + Math.log(Math.exp(a - m) + Math.exp(b - m));
    }
}
