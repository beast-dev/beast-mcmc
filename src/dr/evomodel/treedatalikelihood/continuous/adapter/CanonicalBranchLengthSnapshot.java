package dr.evomodel.treedatalikelihood.continuous.adapter;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;

/**
 * Computes effective branch lengths and caches them for diagnostic/traversal
 * phases that repeatedly request the same child-node values.
 */
final class CanonicalBranchLengthSnapshot {

    private final Tree tree;
    private final BranchRateModel rateModel;
    private final ContinuousRateTransformation rateTransformation;
    private final double[] phaseCache;
    private int phaseDepth;

    CanonicalBranchLengthSnapshot(final Tree tree,
                                  final BranchRateModel rateModel,
                                  final ContinuousRateTransformation rateTransformation) {
        this.tree = tree;
        this.rateModel = rateModel;
        this.rateTransformation = rateTransformation;
        this.phaseCache = new double[tree.getNodeCount()];
        this.phaseDepth = 0;
    }

    void beginPhase() {
        if (phaseDepth++ > 0) {
            return;
        }
        final int rootIndex = tree.getRoot().getNumber();
        for (int childNodeIndex = 0; childNodeIndex < tree.getNodeCount(); childNodeIndex++) {
            if (childNodeIndex == rootIndex) {
                phaseCache[childNodeIndex] = 0.0;
                continue;
            }
            phaseCache[childNodeIndex] = compute(childNodeIndex);
        }
    }

    void endPhase() {
        if (phaseDepth > 0) {
            phaseDepth--;
        }
    }

    double getEffectiveBranchLength(final int childNodeIndex) {
        if (phaseDepth > 0) {
            return phaseCache[childNodeIndex];
        }
        return compute(childNodeIndex);
    }

    private double compute(final int childNodeIndex) {
        final NodeRef node = tree.getNode(childNodeIndex);
        final double rawLength = tree.getBranchLength(node);
        final double normalization = rateTransformation == null ? 1.0 : rateTransformation.getNormalization();
        if (rateModel == null) {
            return rawLength * normalization;
        }
        return rawLength * rateModel.getBranchRate(tree, node) * normalization;
    }
}
