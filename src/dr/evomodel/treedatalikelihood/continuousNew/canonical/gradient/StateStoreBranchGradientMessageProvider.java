package dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.traversal.CanonicalTreeStateStore;
import dr.evomodel.treedatalikelihood.continuousNew.observationmodel.CanonicalTipObservationModel;

/**
 * First branch-gradient message provider: a thin semantic view over the current state store.
 */
public final class StateStoreBranchGradientMessageProvider implements BranchGradientMessageProvider {

    private final Tree tree;
    private final CanonicalTreeStateStore stateStore;

    public StateStoreBranchGradientMessageProvider(final Tree tree,
                                                   final CanonicalTreeStateStore stateStore) {
        this.tree = tree;
        this.stateStore = stateStore;
    }

    @Override
    public void fillBranchMessages(final int childNodeIndex,
                                   final BranchGradientMessageBundle out) {
        final boolean externalChild = tree.isExternal(tree.getNode(childNodeIndex));
        final boolean fixedParent =
                stateStore.hasFixedRootValue && tree.isRoot(tree.getParent(tree.getNode(childNodeIndex)));
        final CanonicalTipObservationModel observationModel =
                externalChild ? stateStore.tipObservationModels[childNodeIndex] : null;
        out.set(
                childNodeIndex,
                externalChild,
                fixedParent,
                stateStore.branchAboveParent[childNodeIndex],
                externalChild ? null : stateStore.postOrder[childNodeIndex],
                observationModel,
                fixedParent ? stateStore.fixedRootValue : null);
    }

    @Override
    public double getRootDiffusionScale() {
        return stateStore.lastRootDiffusionScale;
    }

    @Override
    public void fillRootMessages(final BranchGradientInputs out) {
        final int rootIndex = tree.getRoot().getNumber();
        out.setRoot(
                stateStore.lastRootDiffusionScale,
                stateStore.branchAboveParent[rootIndex],
                stateStore.postOrder[rootIndex]);
    }
}
