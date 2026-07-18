/*
 * MascotEventTape.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;

/**
 * Builds a compact backward-time event list from a BEAST-X TreeModel and an
 * attributePatterns-style PatternList of tip states, resolved by taxon name
 * (not by tree traversal order — see {@link #fromTree}).
 */
public final class MascotEventTape {

    private final MascotCore.Event[] events;
    private final MascotCore.PreparedEvents preparedEvents;

    private MascotEventTape(MascotCore.Event[] events) {
        this.events = events;
        this.preparedEvents = MascotCore.prepareEvents(events);
    }

    public static MascotEventTape fromTree(TreeModel tree, PatternList tipPatterns, int stateCount) {
        int tipCount = tree.getExternalNodeCount();
        int internalCount = tree.getInternalNodeCount();
        if (tipPatterns.getPatternCount() != 1) {
            throw new IllegalArgumentException("tip-state attributePatterns must contain exactly one pattern, found " +
                    tipPatterns.getPatternCount());
        }
        if (tipPatterns.getStateCount() != stateCount) {
            throw new IllegalArgumentException("tip-state attributePatterns dataType has " +
                    tipPatterns.getStateCount() + " states, which does not match mascotLikelihood's stateCount=" +
                    stateCount);
        }

        MascotCore.Event[] events = new MascotCore.Event[tipCount + internalCount];
        int index = 0;
        for (int i = 0; i < tipCount; i++) {
            NodeRef node = tree.getExternalNode(i);
            String taxonId = tree.getTaxonId(i);
            int taxonIndex = tipPatterns.getTaxonIndex(taxonId);
            if (taxonIndex < 0) {
                throw new IllegalArgumentException("tip-state attributePatterns has no entry for " +
                        tipDescription(tree, node, i));
            }
            int state = tipPatterns.getPatternState(taxonIndex, 0);
            if (state < 0 || state >= stateCount) {
                throw new IllegalArgumentException("tip state out of range for " +
                        tipDescription(tree, node, i) + ": " + state +
                        " (stateCount=" + stateCount + ")");
            }
            events[index++] = MascotCore.Event.sample(tree.getNodeHeight(node), node.getNumber(), state);
        }

        for (int i = 0; i < internalCount; i++) {
            NodeRef node = tree.getInternalNode(i);
            int childCount = tree.getChildCount(node);
            if (childCount != 2) {
                throw new IllegalArgumentException("MASCOT currently requires binary trees; node " +
                        node.getNumber() + " has " + childCount + " children");
            }
            NodeRef child1 = tree.getChild(node, 0);
            NodeRef child2 = tree.getChild(node, 1);
            events[index++] = MascotCore.Event.coalescence(
                    tree.getNodeHeight(node),
                    child1.getNumber(),
                    child2.getNumber(),
                    node.getNumber()
            );
        }

        return new MascotEventTape(events);
    }

    private static String tipDescription(TreeModel tree, NodeRef node, int externalNodeIndex) {
        String taxonId = tree.getTaxonId(externalNodeIndex);
        if (taxonId == null) {
            return "external node " + externalNodeIndex + " (node " + node.getNumber() + ")";
        }
        return "external node " + externalNodeIndex + " (node " + node.getNumber() +
                ", taxon " + taxonId + ")";
    }

    /**
     * Builds a {@code MascotCore}-compatible {@code branchRates} array (indexed
     * by lineage id, i.e. {@code NodeRef.getNumber()} -- see {@link #fromTree}),
     * from any {@link BranchRateModel}. The root's slot is left {@code 0.0} and
     * is never read by {@code MascotCore} (the root is never an active lineage).
     * Rebuilt fresh on every call, deliberately uncached, matching {@link
     * MascotDynamics#getThetaValues()}'s existing always-rebuild pattern: this is
     * an O(nodeCount) traversal with no allocation beyond the returned array, so
     * caching would only add an invalidation-tracking axis (tree topology change
     * vs. clock-parameter change) for no real benefit.
     */
    public static double[] buildBranchRates(TreeModel tree, BranchRateModel branchRateModel) {
        double[] rates = new double[tree.getNodeCount()];
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                rates[node.getNumber()] = branchRateModel.getBranchRate(tree, node);
            }
        }
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
