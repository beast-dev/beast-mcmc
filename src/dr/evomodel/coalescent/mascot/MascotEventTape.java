/*
 * MascotEventTape.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.evolution.coalescent.IntervalType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.coalescent.StructuredTipStates;
import dr.evomodel.tree.TreeModel;

/**
 * Builds a compact backward-time event list from a shared {@link
 * BigFastTreeIntervals} (the same tree-interval walk BASTA's {@code
 * CoalescentIntervalTraversal} consumes) and an attributePatterns-style
 * PatternList of tip states, resolved by taxon name (not by tree traversal
 * order — see {@link #fromTree}).
 */
public final class MascotEventTape {

    private final MascotCore.Event[] events;
    private final MascotCore.PreparedEvents preparedEvents;

    private MascotEventTape(MascotCore.Event[] events) {
        this.events = events;
        this.preparedEvents = MascotCore.prepareEvents(events);
    }

    /**
     * @param tipPartials one partial vector per external node, indexed by
     *                    {@code NodeRef.getNumber()} -- see {@link
     *                    StructuredTipStates#buildPartialsCache}. Built once by
     *                    the caller (tip states never change during a run), not
     *                    re-derived from {@code tipPatterns} here on every
     *                    rebuild -- this method used to call {@link
     *                    StructuredTipStates#getPartials} per tip per call,
     *                    which re-did an O(tipCount) taxon-name lookup and a
     *                    fresh allocation for every tip on every tree-changing
     *                    operator, even though the result never differs.
     */
    public static MascotEventTape fromTree(BigFastTreeIntervals treeIntervals, double[][] tipPartials, int stateCount) {
        Tree tree = treeIntervals.getTree();
        int intervalCount = treeIntervals.getIntervalCount();
        // BigFastTreeIntervals always has one more node/event than interval (see
        // its own class doc): the 0th node is the start of the 0th interval, and
        // is reached here via the interval=-1 convention that
        // getSamplingNode/getCoalescentNode already use elsewhere in the codebase
        // (e.g. CoalescentIntervalTraversal). That first node is always a sample
        // event by construction -- BigFastTreeIntervals.calculateIntervals()
        // itself throws if it isn't.
        MascotCore.Event[] events = new MascotCore.Event[intervalCount + 1];
        for (int interval = -1; interval < intervalCount; interval++) {
            boolean isSample = interval == -1 || treeIntervals.getIntervalType(interval) == IntervalType.SAMPLE;
            events[interval + 1] = isSample
                    ? sampleEvent(tree, treeIntervals.getSamplingNode(interval), tipPartials)
                    : coalescentEvent(tree, treeIntervals.getCoalescentNode(interval));
        }

        return new MascotEventTape(events);
    }

    private static MascotCore.Event sampleEvent(Tree tree, NodeRef node, double[][] tipPartials) {
        return MascotCore.Event.sample(tree.getNodeHeight(node), node.getNumber(), tipPartials[node.getNumber()]);
    }

    private static MascotCore.Event coalescentEvent(Tree tree, NodeRef node) {
        int childCount = tree.getChildCount(node);
        if (childCount != 2) {
            throw new IllegalArgumentException("MASCOT currently requires binary trees; node " +
                    node.getNumber() + " has " + childCount + " children");
        }
        NodeRef child1 = tree.getChild(node, 0);
        NodeRef child2 = tree.getChild(node, 1);
        return MascotCore.Event.coalescence(
                tree.getNodeHeight(node),
                child1.getNumber(),
                child2.getNumber(),
                node.getNumber()
        );
    }

    /**
     * Writes a {@code MascotCore}-compatible {@code branchRates} array (indexed
     * by lineage id, i.e. {@code NodeRef.getNumber()} -- see {@link #fromTree}),
     * from any {@link BranchRateModel}, into a caller-owned {@code destination}
     * of length {@code tree.getNodeCount()}. The root's slot is explicitly set
     * to {@code 0.0} (never read by {@code MascotCore}, but every non-root slot
     * must be freshly overwritten on every call regardless -- a topology change
     * can move which node is root, which would otherwise leave a stale nonzero
     * value in a reused buffer's old root slot). Rebuilt fresh on every call,
     * deliberately uncached, matching {@link MascotDynamics#getThetaValues()}'s
     * existing always-rebuild pattern: this is an O(nodeCount) traversal, so
     * caching would only add an invalidation-tracking axis (tree topology
     * change vs. clock-parameter change) for no real benefit.
     */
    public static void writeBranchRates(TreeModel tree, BranchRateModel branchRateModel, double[] destination) {
        if (destination.length != tree.getNodeCount()) {
            throw new IllegalArgumentException("destination dimension " + destination.length +
                    " does not match tree node count " + tree.getNodeCount());
        }
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            destination[node.getNumber()] = tree.isRoot(node) ? 0.0 : branchRateModel.getBranchRate(tree, node);
        }
    }

    /** Allocating compatibility wrapper over {@link #writeBranchRates}; prefer the in-place form in hot paths. */
    public static double[] buildBranchRates(TreeModel tree, BranchRateModel branchRateModel) {
        double[] rates = new double[tree.getNodeCount()];
        writeBranchRates(tree, branchRateModel, rates);
        return rates;
    }

    public MascotCore.Event[] getEvents() {
        return events.clone();
    }

    /**
     * Already sorted and validated; safe to reuse across many evaluations of the
     * same fixed tree without re-cloning or re-sorting.
     */
    public MascotCore.PreparedEvents getPreparedEvents() {
        return preparedEvents;
    }
}
