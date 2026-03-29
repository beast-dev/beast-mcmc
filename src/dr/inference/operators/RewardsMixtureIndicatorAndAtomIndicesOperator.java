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

    private final int[] branchBuffer;

    private final RewardsAwareBranchModel rewardsAwareBranchModel;
    private final TreeDataLikelihood treeDataLikelihood;
    private final DiscreteDataLikelihoodDelegate discreteDelegate;

    private final Tree tree;
    private final int nstates;

    private final double[] prePartial;
    private final double[] postPartial;
    private final double[] atomicWeightsScratch;

    private final boolean DEBUG = false;

    private final double[] preScales = new double[1];
    private final double[] postScales = new double[1];
    private final double[] logAtomicWeights;

    public RewardsMixtureIndicatorAndAtomIndicesOperator(
            final Parameter indicatorZ,
            final Parameter atomIndex,
            final RewardsAwareBranchModel rewardsAwareBranchModel,
            final TreeDataLikelihood treeDataLikelihood,
            final double updateProportion,
            final boolean adaptUpdateProportion,
            final double weight
    ) {
        super(adaptUpdateProportion ? AdaptationMode.ADAPTATION_ON : ADAPTATION_OFF);

        this.indicatorZ = indicatorZ;
        this.atomIndex = atomIndex;
        this.rewardsAwareBranchModel = rewardsAwareBranchModel;
        this.treeDataLikelihood = treeDataLikelihood;
        this.borderBuffer = new int[treeDataLikelihood.getTree().getNodeCount()]; // max possible size

        final DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (!(delegate instanceof DiscreteDataLikelihoodDelegate)) {
            throw new IllegalArgumentException(
                    "RewardsMixtureIndicatorAndAtomIndicesOperator requires TreeDataLikelihood to use DiscreteDataLikelihoodDelegate"
            );
        }
        this.discreteDelegate = (DiscreteDataLikelihoodDelegate) delegate;

        this.tree = rewardsAwareBranchModel.getTree();
        this.nstates = rewardsAwareBranchModel.getStateCount();

        this.updateProportion = updateProportion;
        this.adaptUpdateProportion = adaptUpdateProportion;
        this.adaptableParameter = logit(updateProportion);

        if (atomIndex.getDimension() != indicatorZ.getDimension()) {
            throw new IllegalArgumentException("atomIndex and indicatorZ must have the same dimension.");
        }

        final int B = indicatorZ.getDimension();
        this.branchBuffer = new int[B];
        for (int i = 0; i < B; i++) {
            branchBuffer[i] = i;
        }

        this.prePartial = new double[nstates];
        this.postPartial = new double[nstates];
        this.atomicWeightsScratch = new double[nstates];

        this.logAtomicWeights = new double[nstates];

        setWeight(weight);
    }

    @Override
    public String getOperatorName() {
        return "RewardsMixtureIndicatorAndAtomIndicesOperator";
    }

    @Override
    public double doOperation() {

        treeDataLikelihood.calculatePostOrderStatistics();
        discreteDelegate.ensurePreOrderComputed();

        final int B = indicatorZ.getDimension();
        final int nToUpdate = 1; // TODO temporary
//                Math.max(1, (int) Math.round(updateProportion * B));

        if (nToUpdate >= B) {
            for (int b = 0; b < B; b++) {
                resampleBranch(b);
            }
            atomIndex.fireParameterChangedEvent();
            indicatorZ.fireParameterChangedEvent();
            return 0.0;
        }

        for (int i = 0; i < B; i++) {
            branchBuffer[i] = i;
        }
        int branchNumber = chooseCandidateBranch(tree, indicatorZ);
        int b = rewardsAwareBranchModel.getParameterIndexForNode(branchNumber);
        branchBuffer[b] = branchNumber;
        resampleBranch(b);

//        for (int i = 0; i < nToUpdate; i++) {
//            final int j = i + MathUtils.nextInt(B - i);
//            final int tmp = branchBuffer[i];
//            branchBuffer[i] = branchBuffer[j];
//            branchBuffer[j] = tmp;
//
//            resampleBranch(branchBuffer[i]);
//        }

        atomIndex.fireParameterChangedEvent();
        indicatorZ.fireParameterChangedEvent();

        return 0.0;
    }
    private void resampleBranch(final int b) {

        if (b < 0 || b >= branchBuffer.length) {
            throw new IllegalArgumentException("Branch index out of range: " + b);
        }
        if (b >= indicatorZ.getDimension()) {
            throw new IllegalArgumentException("Branch index out of range for indicatorZ: " + b);
        }
        if (b >= atomIndex.getDimension()) {
            throw new IllegalArgumentException("Branch index out of range for atomIndex: " + b);
        }

        final int branchNodeNumber = branchBuffer[b];

        final RewardsMixtureBranchResamplingHelper.BranchWeights weights =
                computeBranchWeights(branchNodeNumber);

        final int newIndicator =
                RewardsMixtureBranchResamplingHelper.sampleIndicatorFromLogs(
                        weights.logAtomicTotalWeight,
                        weights.logCtsWeight
                );

        final int currentIndicator = (int) Math.round(indicatorZ.getParameterValue(b));
        final int currentAtom = (int) Math.round(atomIndex.getParameterValue(b));

        if (newIndicator == 1) {
            final int newAtom =
                    RewardsMixtureBranchResamplingHelper.sampleAtomFromLogs(
                            weights.logAtomicWeights,
                            nstates
                    );

            // update atom index first, then indicator
            // so that when the branch becomes atomic it already has a valid state
            if (currentAtom != newAtom) {
                atomIndex.setParameterValueQuietly(b, newAtom);
            }
            if (currentIndicator != 1) {
                indicatorZ.setParameterValueQuietly(b, 1.0);
            }

        } else {
            // continuous branch: keep atomIndex unchanged
            if (currentIndicator != 0) {
                indicatorZ.setParameterValueQuietly(b, 0.0);
            }
        }

        // fire one event per parameter at the end
        atomIndex.fireParameterChangedEvent();
        indicatorZ.fireParameterChangedEvent();
    }


    private void loadBranchPartials(final int nodeNum,
                                    final double[] prePartialOut,
                                    final double[] postPartialOut) {
        Arrays.fill(prePartialOut, 0.0);
        Arrays.fill(postPartialOut, 0.0);

        discreteDelegate.getPreOrderAtBranchStartInto(nodeNum, prePartialOut);
        discreteDelegate.getPostOrderAtBranchEndInto(nodeNum, postPartialOut);
    }
    private RewardsMixtureBranchResamplingHelper.BranchWeights computeBranchWeights(int branchNodeNumber) {

        // Load normalized preorder/postorder branch messages
        loadBranchPartials(branchNodeNumber, prePartial, postPartial);

        // Load associated log-scales
        discreteDelegate.getPreOrderBranchScalesInto(branchNodeNumber, preScales);
        discreteDelegate.getPostOrderBranchScalesInto(branchNodeNumber, postScales);

        final double preScale = preScales[0];
        final double postScale = postScales[0];

        // Atomic local factor: no-jump mass = exp(-lambda * t)
        final NodeRef node = tree.getNode(branchNodeNumber);
        final double branchLength = tree.getBranchLength(node);
        final double lambda = rewardsAwareBranchModel.getUniformizationRate();
        final double logAtomicLocalFactor = -lambda * branchLength;

        // Atomic weights by state
        for (int j = 0; j < nstates; j++) {
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

        // Continuous weight
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
                    .append(", lambda=").append(lambda)
                    .append(", preScale=").append(preScale)
                    .append(", postScale=").append(postScale)
                    .append(", prePartial=").append(java.util.Arrays.toString(prePartial))
                    .append(", postPartial=").append(java.util.Arrays.toString(postPartial))
                    .append(", logAtomicWeights=").append(java.util.Arrays.toString(logAtomicWeights));
            throw new IllegalStateException(sb.toString());
        }

        return new RewardsMixtureBranchResamplingHelper.BranchWeights(
                java.util.Arrays.copyOf(logAtomicWeights, nstates),
                logAtomicTotalWeight,
                logCtsWeight
        );
    }

    private final int[] borderBuffer;
    private int borderCount = 0;

    private void collectBorderBranches(Tree tree, Parameter indicatorZ) {
        borderCount = 0;

        for (int b = 0; b < tree.getNodeCount(); b++) {
            NodeRef node = tree.getNode(b);

            if (tree.isRoot(node)) {
                continue; // no branch above the root
            }

            final int zb = (int) indicatorZ.getParameterValue(b);
            boolean isBorder = false;

            // check parent-adjacent branch
            NodeRef parent = tree.getParent(node);
            if (parent != null && !tree.isRoot(parent)) {
                int parentBranch = parent.getNumber();
                int zParent = (int) indicatorZ.getParameterValue(parentBranch);
                if (zParent != zb) {
                    isBorder = true;
                }
            }

            // check child-adjacent branches
            if (!isBorder) {
                int childCount = tree.getChildCount(node);
                for (int i = 0; i < childCount; i++) {
                    NodeRef child = tree.getChild(node, i);
                    int childBranch = child.getNumber();
                    int zChild = (int) indicatorZ.getParameterValue(childBranch);
                    if (zChild != zb) {
                        isBorder = true;
                        break;
                    }
                }
            }

            if (isBorder) {
                borderBuffer[borderCount++] = b;
            }
        }
    }
    private int chooseCandidateBranch(Tree tree, Parameter indicatorZ) {
        collectBorderBranches(tree, indicatorZ);

        // mostly choose from border, occasionally from all branches
        if (borderCount > 0 && MathUtils.nextDouble() < 0.5) {
            return borderBuffer[MathUtils.nextInt(borderCount)];
        }

        // fallback global draw for irreducibility
        int b;
        do {
            b = MathUtils.nextInt(tree.getNodeCount());
        } while (tree.isRoot(tree.getNode(b)));

        return b;
    }

//    private RewardsMixtureBranchResamplingHelper.BranchWeights computeBranchWeights(
//            final NodeRef node,
//            final int nodeNum,
//            final double[] prePartial,
//            final double[] postPartial,
//            final double[] atomicWeightsOut
//    ) {
//        Arrays.fill(atomicWeightsOut, 0.0);
//
//        computeAtomicWeightsForBranchInto(node, prePartial, postPartial, atomicWeightsOut);
//
//        double atomicTotal = 0.0;
//        for (int j = 0; j < nstates; j++) {
//            final double w = atomicWeightsOut[j];
//            if (Double.isNaN(w) || Double.isInfinite(w) || w < 0.0) {
//                throw new IllegalStateException("Invalid atomic weight at state " + j + ": " + w);
//            }
//            atomicTotal += w;
//        }
//
//        final double ctsWeight = computeCtsWeightForBranch(nodeNum, prePartial, postPartial);
//        if (Double.isNaN(ctsWeight) || Double.isInfinite(ctsWeight) || ctsWeight < 0.0) {
//            throw new IllegalStateException("Invalid continuous weight: " + ctsWeight);
//        }
//        if (DEBUG) {
//            System.out.println("Atomic weights: " + Arrays.toString(atomicWeightsScratch));
//            System.out.println("Atomic weight: " + atomicTotal);
//            System.out.println("Cts weight: " + ctsWeight);
//        }
//
//        final double denom = atomicTotal + ctsWeight;
//        if (!(denom > 0.0) || Double.isInfinite(denom)) {
//            throw new IllegalStateException(
//                    "Invalid total weight: atomic=" + atomicTotal + ", continuous=" + ctsWeight
//            );
//        }
//
//        return new RewardsMixtureBranchResamplingHelper.BranchWeights(
//                atomicWeightsOut,
//                atomicTotal,
//                ctsWeight
//        );
//    }

    private void computeAtomicWeightsForBranchInto(final NodeRef node,
                                                   double[] prePartial,
                                                   double[] postPartial,
                                                   double[] atomicWeightsOut) {
        final double uniformizationRate = rewardsAwareBranchModel.getUniformizationRate();
        final double branchLength = tree.getBranchLength(node);
        final double scale = Math.exp(-uniformizationRate * branchLength);

        for (int j = 0; j < nstates; j++) {
            atomicWeightsOut[j] = scale * prePartial[j] * postPartial[j];
        }
    }

    private double computeCtsWeightForBranch(final int nodeNum,
                                             double[] prePartial,
                                             double[] postPartial) {
        final double[] W = rewardsAwareBranchModel.getTransitionMatrixCts(nodeNum);
        return RewardsMixtureBranchResamplingHelper.bilinearFormStable(prePartial, W, postPartial, nstates);
    }

    @Override
    protected void setAdaptableParameterValue(final double value) {
        if (!adaptUpdateProportion) return;
        adaptableParameter = value;
        updateProportion = logistic(value);

        final double eps = 1e-6;
        if (updateProportion < eps) updateProportion = eps;
        if (updateProportion > 1.0) updateProportion = 1.0;
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
        StringBuilder sb = new StringBuilder();

        final int dim = indicatorZ.getDimension();

        int activeCount = 0;
        int inactiveCount;
        int sumAtomIndex = 0;
        int sumActiveAtomIndex = 0;

        sb.append("RewardsMixtureIndicatorAndAtomIndicesOperator\n");

        sb.append("dimension: ").append(dim).append("\n");

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

//        sb.append("lastBranchIndex: ").append(lastBranchIndex).append("\n");
//        sb.append("lastNodeNumber: ").append(lastNodeNumber).append("\n");

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