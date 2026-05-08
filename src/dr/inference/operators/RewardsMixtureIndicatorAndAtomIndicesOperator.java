package dr.inference.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DiscreteDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.Reportable;

import java.util.Arrays;

import static dr.inference.operators.AdaptationMode.ADAPTATION_OFF;

/*
 * @author Filippo Monti
 */
public final class RewardsMixtureIndicatorAndAtomIndicesOperator extends AbstractAdaptableOperator implements Reportable {

    private final Parameter indicatorZ;
    private final Parameter atomIndex;

    private double updateProportion;
    private final boolean adaptUpdateProportion;
    private double adaptableParameter;

    private final int branchCount;
    private final int[] nodeNumberByParameterIndex;
    private final int[] candidateParameterBuffer;
    private final int[] borderParameterBuffer;
    private int borderCount = 0;

    private final RewardsAwareBranchModel rewardsAwareBranchModel;
    private final TreeDataLikelihood treeDataLikelihood;
    private final DiscreteDataLikelihoodDelegate discreteDelegate;

    private final Tree tree;
    private final int nstates;

    private final double[] prePartial;
    private final double[] postPartial;
    private final double[] preScales = new double[1];
    private final double[] postScales = new double[1];
    private final double[] logAtomicWeights;

    private final boolean useClusterMoves;
    private final int clusterSize;
    private final double clusterBorderBias;
    private final RewardsMixtureClusterResamplingHelper clusterHelper;

    public RewardsMixtureIndicatorAndAtomIndicesOperator(
            final Parameter indicatorZ,
            final Parameter atomIndex,
            final RewardsAwareBranchModel rewardsAwareBranchModel,
            final TreeDataLikelihood treeDataLikelihood,
            final double updateProportion,
            final boolean adaptUpdateProportion,
            final boolean useClusterMoves,
            final int clusterSize,
            final double clusterBorderBias,
            final double weight
    ) {
        super(adaptUpdateProportion ? AdaptationMode.ADAPTATION_ON : ADAPTATION_OFF);

        if (indicatorZ == null) {
            throw new IllegalArgumentException("indicatorZ must be non-null");
        }
        if (atomIndex == null) {
            throw new IllegalArgumentException("atomIndex must be non-null");
        }
        if (rewardsAwareBranchModel == null) {
            throw new IllegalArgumentException("rewardsAwareBranchModel must be non-null");
        }
        if (treeDataLikelihood == null) {
            throw new IllegalArgumentException("treeDataLikelihood must be non-null");
        }

        this.indicatorZ = indicatorZ;
        this.atomIndex = atomIndex;
        this.rewardsAwareBranchModel = rewardsAwareBranchModel;
        this.treeDataLikelihood = treeDataLikelihood;

        final DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (!(delegate instanceof DiscreteDataLikelihoodDelegate)) {
            throw new IllegalArgumentException(
                    "RewardsMixtureIndicatorAndAtomIndicesOperator requires TreeDataLikelihood to use DiscreteDataLikelihoodDelegate"
            );
        }
        this.discreteDelegate = (DiscreteDataLikelihoodDelegate) delegate;

        this.tree = rewardsAwareBranchModel.getTree();
        this.nstates = rewardsAwareBranchModel.getStateCount();

        if (treeDataLikelihood.getTree().getNodeCount() != tree.getNodeCount()) {
            throw new IllegalArgumentException(
                    "TreeDataLikelihood and RewardsAwareBranchModel must use trees with the same node count."
            );
        }
        if (updateProportion <= 0.0 || updateProportion > 1.0) {
            throw new IllegalArgumentException("updateProportion must be in (0, 1]. Found: " + updateProportion);
        }

        this.updateProportion = updateProportion;
        this.adaptUpdateProportion = adaptUpdateProportion;
        this.adaptableParameter = logit(updateProportion);

        if (atomIndex.getDimension() != indicatorZ.getDimension()) {
            throw new IllegalArgumentException("atomIndex and indicatorZ must have the same dimension.");
        }

        this.branchCount = indicatorZ.getDimension();
        if (branchCount != tree.getNodeCount() - 1) {
            throw new IllegalArgumentException(
                    "indicatorZ dimension must equal the number of non-root branches. Found dimension=" +
                            branchCount + ", expected=" + (tree.getNodeCount() - 1)
            );
        }

        this.nodeNumberByParameterIndex = new int[branchCount];
        this.candidateParameterBuffer = new int[branchCount];
        this.borderParameterBuffer = new int[branchCount];
        initializeBranchMappings();

        this.prePartial = new double[nstates];
        this.postPartial = new double[nstates];
        this.logAtomicWeights = new double[nstates];

        this.useClusterMoves = useClusterMoves;
        this.clusterSize = Math.max(1, clusterSize);
        this.clusterBorderBias = Math.max(0.0, Math.min(1.0, clusterBorderBias));

        this.clusterHelper = useClusterMoves
                ? new RewardsMixtureClusterResamplingHelper(
                indicatorZ,
                atomIndex,
                rewardsAwareBranchModel,
                treeDataLikelihood,
                discreteDelegate,
                this.clusterSize,
                this.clusterBorderBias,
                new RewardsMixtureClusterResamplingHelper.ExactLogTargetEvaluator() {
                    @Override
                    public double computeCurrentLogTarget(final int[] clusterBranchNodeNumbers, final int clusterCount) {
                        return RewardsMixtureIndicatorAndAtomIndicesOperator.this.treeDataLikelihood.getLogLikelihood();
                    }
                }
        )
                : null;

        setWeight(weight);
    }

    @Override
    public String getOperatorName() {
        return "RewardsMixtureIndicatorAndAtomIndicesOperator";
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
                                "; branch parameter dimension is " + branchCount
                );
            }
            if (nodeNumberByParameterIndex[parameterIndex] != -1) {
                throw new IllegalArgumentException(
                        "Multiple non-root nodes map to branch parameter index " + parameterIndex
                );
            }

            nodeNumberByParameterIndex[parameterIndex] = nodeNumber;
            observedBranches++;
        }

        if (observedBranches != branchCount) {
            throw new IllegalArgumentException(
                    "Observed " + observedBranches + " non-root branches, but indicatorZ dimension is " + branchCount
            );
        }
        for (int i = 0; i < branchCount; i++) {
            if (nodeNumberByParameterIndex[i] < 0) {
                throw new IllegalArgumentException("No branch node maps to parameter index " + i);
            }
        }
    }

    @Override
    public double doOperation() {

        initializeBranchMappings();

        discreteDelegate.updatePostOrdersFromTreeDataLikelihood(treeDataLikelihood);
        discreteDelegate.ensurePreOrderComputed();

        if (useClusterMoves) {
            return doClusterOperation();
        }

        return doSingleBranchOperation();
    }

    private double doSingleBranchOperation() {

        final int nToUpdate = Math.max(1, Math.min(branchCount, (int) Math.round(updateProportion * branchCount)));

        if (nToUpdate >= branchCount) {
            for (int parameterIndex = 0; parameterIndex < branchCount; parameterIndex++) {
                resampleBranchQuietly(parameterIndex, nodeNumberByParameterIndex[parameterIndex]);
            }
            atomIndex.fireParameterChangedEvent();
            indicatorZ.fireParameterChangedEvent();
            return 0.0;
        }

        if (nToUpdate == 1) {
            final int parameterIndex = chooseCandidateParameterIndex();
            resampleBranchQuietly(parameterIndex, nodeNumberByParameterIndex[parameterIndex]);
        } else {
            for (int i = 0; i < branchCount; i++) {
                candidateParameterBuffer[i] = i;
            }
            for (int i = 0; i < nToUpdate; i++) {
                final int j = i + MathUtils.nextInt(branchCount - i);
                final int tmp = candidateParameterBuffer[i];
                candidateParameterBuffer[i] = candidateParameterBuffer[j];
                candidateParameterBuffer[j] = tmp;

                final int parameterIndex = candidateParameterBuffer[i];
                resampleBranchQuietly(parameterIndex, nodeNumberByParameterIndex[parameterIndex]);
            }
        }

        atomIndex.fireParameterChangedEvent();
        indicatorZ.fireParameterChangedEvent();

        return 0.0;
    }

    private double doClusterOperation() {
        clusterHelper.doOneClusterMoveAndFireEvents();
        return 0.0;
    }

    private void resampleBranchQuietly(final int branchIndex, final int branchNodeNumber) {

        if (branchIndex < 0 || branchIndex >= indicatorZ.getDimension()) {
            throw new IllegalArgumentException("Branch parameter index out of range for indicatorZ: " + branchIndex);
        }
        if (branchIndex < 0 || branchIndex >= atomIndex.getDimension()) {
            throw new IllegalArgumentException("Branch parameter index out of range for atomIndex: " + branchIndex);
        }

        final RewardsMixtureBranchResamplingHelper.BranchWeights weights =
                computeBranchWeights(branchNodeNumber);

        final int newIndicator =
                RewardsMixtureBranchResamplingHelper.sampleIndicatorFromLogs(
                        weights.logAtomicTotalWeight,
                        weights.logCtsWeight
                );

        final int currentIndicator = (int) Math.round(indicatorZ.getParameterValue(branchIndex));
        final int currentAtom = (int) Math.round(atomIndex.getParameterValue(branchIndex));

        if (newIndicator == 1) {
            final int newAtom =
                    RewardsMixtureBranchResamplingHelper.sampleAtomFromLogs(
                            weights.logAtomicWeights,
                            nstates
                    );

            if (currentAtom != newAtom) {
                atomIndex.setParameterValueQuietly(branchIndex, newAtom);
            }
            if (currentIndicator != 1) {
                indicatorZ.setParameterValueQuietly(branchIndex, 1.0);
            }

        } else if (currentIndicator != 0) {
            indicatorZ.setParameterValueQuietly(branchIndex, 0.0);
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
            logAtomicWeights[j] =
                    RewardsMixtureBranchResamplingHelper.logAtomicWeight(
                            prePartial[j],
                            postPartial[j],
                            logAtomicLocalFactor,
                            preScale,
                            postScale
                    );
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
                );

        if (!Double.isFinite(logAtomicTotalWeight) && !Double.isFinite(logCtsWeight)) {
            final StringBuilder sb = new StringBuilder();
            sb.append("Invalid total weight: atomic=").append(logAtomicTotalWeight)
                    .append(", continuous=").append(logCtsWeight)
                    .append(", branchNodeNumber=").append(branchNodeNumber)
                    .append(", branchLength=").append(branchLength)
                    .append(", preScale=").append(preScale)
                    .append(", postScale=").append(postScale)
                    .append(", prePartial=").append(Arrays.toString(prePartial))
                    .append(", postPartial=").append(Arrays.toString(postPartial))
                    .append(", logAtomicWeights=").append(Arrays.toString(logAtomicWeights));
            throw new IllegalStateException(sb.toString());
        }

        return new RewardsMixtureBranchResamplingHelper.BranchWeights(
                Arrays.copyOf(logAtomicWeights, nstates),
                logAtomicTotalWeight,
                logCtsWeight
        );
    }

    private void collectBorderBranches() {
        borderCount = 0;

        for (int parameterIndex = 0; parameterIndex < branchCount; parameterIndex++) {
            final int nodeNumber = nodeNumberByParameterIndex[parameterIndex];
            final NodeRef node = tree.getNode(nodeNumber);

            final int z = getIndicatorState(parameterIndex);
            boolean isBorder = false;

            final NodeRef parent = tree.getParent(node);
            if (parent != null && !tree.isRoot(parent)) {
                final int parentParameterIndex = rewardsAwareBranchModel.getParameterIndexForNode(parent.getNumber());
                if (getIndicatorState(parentParameterIndex) != z) {
                    isBorder = true;
                }
            }

            if (!isBorder) {
                final int childCount = tree.getChildCount(node);
                for (int i = 0; i < childCount; i++) {
                    final NodeRef child = tree.getChild(node, i);
                    final int childParameterIndex = rewardsAwareBranchModel.getParameterIndexForNode(child.getNumber());
                    if (getIndicatorState(childParameterIndex) != z) {
                        isBorder = true;
                        break;
                    }
                }
            }

            if (isBorder) {
                borderParameterBuffer[borderCount++] = parameterIndex;
            }
        }
    }

    private int chooseCandidateParameterIndex() {
        collectBorderBranches();

        if (borderCount > 0 && MathUtils.nextDouble() < 0.5) {
            return borderParameterBuffer[MathUtils.nextInt(borderCount)];
        }

        return MathUtils.nextInt(branchCount);
    }

    private int getIndicatorState(final int parameterIndex) {
        if (parameterIndex < 0 || parameterIndex >= indicatorZ.getDimension()) {
            throw new IllegalArgumentException("Indicator parameter index out of range: " + parameterIndex);
        }
        final double raw = indicatorZ.getParameterValue(parameterIndex);
        final int value = (int) Math.round(raw);
        if (Math.abs(raw - value) > 1.0e-9 || (value != 0 && value != 1)) {
            throw new IllegalArgumentException(
                    "indicatorZ must contain 0/1 values, found " + raw +
                            " at parameter index " + parameterIndex
            );
        }
        return value;
    }

    @Override
    protected void setAdaptableParameterValue(final double value) {
        if (!adaptUpdateProportion) {
            return;
        }
        adaptableParameter = value;
        updateProportion = logistic(value);

        final double eps = 1e-6;
        if (updateProportion < eps) {
            updateProportion = eps;
        }
        if (updateProportion > 1.0) {
            updateProportion = 1.0;
        }
    }

    @Override
    protected double getAdaptableParameterValue() {
        return adaptableParameter;
    }

    @Override
    public String getAdaptableParameterName() {
        return "updateProportion";
    }

    @Override
    public double getRawParameter() {
        return updateProportion;
    }

    private static double logistic(final double x) {
        if (x >= 0.0) {
            final double e = Math.exp(-x);
            return 1.0 / (1.0 + e);
        } else {
            final double e = Math.exp(x);
            return e / (1.0 + e);
        }
    }

    private static double logit(final double p) {
        final double eps = 1e-12;
        final double pp = Math.max(eps, Math.min(1.0 - eps, p));
        return Math.log(pp / (1.0 - pp));
    }

    @Override
    public String getReport() {
        initializeBranchMappings();

        final StringBuilder sb = new StringBuilder();

        final int dim = indicatorZ.getDimension();

        int activeCount = 0;
        final int inactiveCount;
        int sumAtomIndex = 0;
        int sumActiveAtomIndex = 0;

        sb.append("RewardsMixtureIndicatorAndAtomIndicesOperator\n");
        sb.append("dimension: ").append(dim).append("\n");
        sb.append("useClusterMoves: ").append(useClusterMoves).append("\n");
        sb.append("clusterSize: ").append(clusterSize).append("\n");
        sb.append("clusterBorderBias: ").append(clusterBorderBias).append("\n");

        sb.append("nodeNumberByParameterIndex: [");
        for (int i = 0; i < dim; i++) {
            if (i > 0) {
                sb.append(" ");
            }
            sb.append(nodeNumberByParameterIndex[i]);
        }
        sb.append("]\n");

        sb.append("indicatorZ: [");
        for (int i = 0; i < dim; i++) {
            if (i > 0) {
                sb.append(" ");
            }
            final int z = (int) indicatorZ.getParameterValue(i);
            sb.append(z);
            activeCount += z;
        }
        sb.append("]\n");

        inactiveCount = dim - activeCount;

        sb.append("atomIndex: [");
        for (int i = 0; i < dim; i++) {
            if (i > 0) {
                sb.append(" ");
            }
            final int a = (int) atomIndex.getParameterValue(i);
            sb.append(a);
            sumAtomIndex += a;
            if ((int) indicatorZ.getParameterValue(i) == 1) {
                sumActiveAtomIndex += a;
            }
        }
        sb.append("]\n");

        sb.append("activeCount: ").append(activeCount).append("\n");
        sb.append("inactiveCount: ").append(inactiveCount).append("\n");

        sb.append("activeBranches: [");
        boolean first = true;
        for (int i = 0; i < dim; i++) {
            if ((int) indicatorZ.getParameterValue(i) == 1) {
                if (!first) {
                    sb.append(" ");
                }
                sb.append(i);
                first = false;
            }
        }
        sb.append("]\n");

        sb.append("inactiveBranches: [");
        first = true;
        for (int i = 0; i < dim; i++) {
            if ((int) indicatorZ.getParameterValue(i) == 0) {
                if (!first) {
                    sb.append(" ");
                }
                sb.append(i);
                first = false;
            }
        }
        sb.append("]\n");

        sb.append("activeAtomIndex: [");
        first = true;
        for (int i = 0; i < dim; i++) {
            if ((int) indicatorZ.getParameterValue(i) == 1) {
                if (!first) {
                    sb.append(" ");
                }
                sb.append((int) atomIndex.getParameterValue(i));
                first = false;
            }
        }
        sb.append("]\n");

        sb.append("inactiveAtomIndex: [");
        first = true;
        for (int i = 0; i < dim; i++) {
            if ((int) indicatorZ.getParameterValue(i) == 0) {
                if (!first) {
                    sb.append(" ");
                }
                sb.append((int) atomIndex.getParameterValue(i));
                first = false;
            }
        }
        sb.append("]\n");

        sb.append("sumAtomIndex: ").append(sumAtomIndex).append("\n");
        sb.append("sumActiveAtomIndex: ").append(sumActiveAtomIndex).append("\n");

        sb.append("branchStates: [");
        for (int i = 0; i < dim; i++) {
            if (i > 0) {
                sb.append(" ");
            }
            sb.append(i)
                    .append(":")
                    .append((int) indicatorZ.getParameterValue(i))
                    .append("/")
                    .append((int) atomIndex.getParameterValue(i));
        }
        sb.append("]\n");

        sb.append("lastPrePartial: [");
        for (int i = 0; i < prePartial.length; i++) {
            if (i > 0) {
                sb.append(" ");
            }
            sb.append(prePartial[i]);
        }
        sb.append("]\n");

        sb.append("lastPostPartial: [");
        for (int i = 0; i < postPartial.length; i++) {
            if (i > 0) {
                sb.append(" ");
            }
            sb.append(postPartial[i]);
        }
        sb.append("]\n");

        sb.append("updateProportion: ").append(updateProportion).append("\n");
        sb.append("adaptableParameter: ").append(adaptableParameter).append("\n");

        return sb.toString();
    }
}
