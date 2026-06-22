package dr.evomodel.treedatalikelihood.continuous.canonical.adapter;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.treedatalikelihood.ProcessOnTreeDelegate.BranchOperation;
import dr.evomodel.treedatalikelihood.continuous.ContinuousBranchTimeModel;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;

/**
 * Computes effective branch lengths and caches them for diagnostic/traversal
 * phases that repeatedly request the same child-node values.
 */
final class CanonicalBranchLengthSnapshot implements ContinuousBranchTimeModel {

    private final Tree tree;
    private final BranchRateModel rateModel;
    private final ContinuousRateTransformation rateTransformation;
    private final double[] phaseCache;
    private volatile int phaseDepth;
    private double normalization;

    CanonicalBranchLengthSnapshot(final Tree tree,
                                  final BranchRateModel rateModel,
                                  final ContinuousRateTransformation rateTransformation) {
        this.tree = tree;
        this.rateModel = rateModel;
        this.rateTransformation = rateTransformation;
        this.phaseCache = new double[tree.getNodeCount()];
        this.phaseDepth = 0;
        this.normalization = 1.0;
    }

    void beginPhase() {
        if (phaseDepth > 0) {
            phaseDepth++;
            return;
        }
        updateNormalization();
        final int rootIndex = tree.getRoot().getNumber();
        for (int childNodeIndex = 0; childNodeIndex < tree.getNodeCount(); childNodeIndex++) {
            if (childNodeIndex == rootIndex) {
                phaseCache[childNodeIndex] = 0.0;
                continue;
            }
            phaseCache[childNodeIndex] = compute(childNodeIndex);
        }
        phaseDepth = 1;
    }

    void endPhase() {
        if (phaseDepth > 0) {
            phaseDepth--;
        }
    }

    double getEffectiveBranchLength(final int childNodeIndex) {
        return getEffectiveBranchTime(childNodeIndex);
    }

    @Override
    public double updateNormalization() {
        normalization = rateTransformation == null ? 1.0 : rateTransformation.getNormalization();
        return normalization;
    }

    @Override
    public double getNormalization() {
        return normalization;
    }

    @Override
    public double getEffectiveBranchTime(final BranchOperation branchOperation) {
        return branchOperation.getBranchLength() * normalization;
    }

    @Override
    public double getEffectiveBranchTime(final int childNodeIndex) {
        if (phaseDepth > 0) {
            return phaseCache[childNodeIndex];
        }
        updateNormalization();
        return compute(childNodeIndex);
    }

    @Override
    public int getSegmentCount(final int childNodeIndex) {
        return tree.isRoot(tree.getNode(childNodeIndex)) ? 0 : 1;
    }

    @Override
    public int getSegmentProcessIndex(final int childNodeIndex, final int segmentIndex) {
        requireSingleSegment(childNodeIndex, segmentIndex);
        return 0;
    }

    @Override
    public double getSegmentEffectiveTime(final int childNodeIndex, final int segmentIndex) {
        requireSingleSegment(childNodeIndex, segmentIndex);
        return getEffectiveBranchTime(childNodeIndex);
    }

    private double compute(final int childNodeIndex) {
        final NodeRef node = tree.getNode(childNodeIndex);
        if (tree.isRoot(node)) {
            return 0.0;
        }
        final double rawLength = tree.getBranchLength(node);
        if (rateModel == null) {
            return rawLength * normalization;
        }
        return rawLength * rateModel.getBranchRate(tree, node) * normalization;
    }

    private void requireSingleSegment(final int childNodeIndex, final int segmentIndex) {
        if (segmentIndex != 0 || getSegmentCount(childNodeIndex) != 1) {
            throw new IllegalArgumentException("Invalid one-segment branch time request for node "
                    + childNodeIndex + " segment " + segmentIndex);
        }
    }
}
