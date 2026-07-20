/*
 * StructuredTipStates.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;

/**
 * Shared tip-state lookup for structured-coalescent likelihoods that read one
 * attributePatterns-style {@link PatternList} by tree taxon id.
 */
public final class StructuredTipStates {

    private StructuredTipStates() {
    }

    public static void validateSinglePattern(PatternList tipPatterns, int stateCount, String context) {
        if (tipPatterns.getPatternCount() != 1) {
            throw new IllegalArgumentException(context + " must contain exactly one pattern, found " +
                    tipPatterns.getPatternCount());
        }
        if (tipPatterns.getStateCount() != stateCount) {
            throw new IllegalArgumentException(context + " dataType has " + tipPatterns.getStateCount() +
                    " states, which does not match stateCount=" + stateCount);
        }
    }

    /**
     * Builds every external node's partial vector once, indexed by {@code
     * NodeRef.getNumber()} (0..externalNodeCount-1 by BEAST's tree-numbering
     * convention). Tip states are fixed input data -- like the sequence
     * alignment, never an estimated model variable -- so unlike the tree
     * itself, this never needs to be rebuilt after construction: callers that
     * currently re-derive a tip's partial vector on every tree-topology change
     * (e.g. {@code MascotEventTape.fromTree}, invoked on every tree-changing
     * operator) should build this once and index into it instead.
     */
    public static double[][] buildPartialsCache(Tree tree, PatternList tipPatterns, int stateCount,
                                                boolean allowAmbiguities, String context) {
        double[][] cache = new double[tree.getExternalNodeCount()][];
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef node = tree.getExternalNode(i);
            cache[node.getNumber()] = getPartials(tree, node, tipPatterns, stateCount, allowAmbiguities, context);
        }
        return cache;
    }

    public static double[] getPartials(Tree tree, NodeRef node, PatternList tipPatterns, int stateCount,
                                       boolean allowAmbiguities, String context) {
        int state = getState(tree, node, tipPatterns, context);
        DataType dataType = tipPatterns.getDataType();
        double[] partials = new double[stateCount];

        if (dataType.isAmbiguousState(state)) {
            if (!allowAmbiguities) {
                throw new IllegalArgumentException("ambiguous tip state for " +
                        tipDescription(node, taxonId(tree, node)) + " is not allowed in " + context + ": " + state);
            }
            boolean[] stateSet = dataType.getStateSet(state);
            int possibleStateCount = 0;
            for (int i = 0; i < stateCount; i++) {
                if (stateSet[i]) {
                    partials[i] = 1.0;
                    possibleStateCount++;
                }
            }
            if (possibleStateCount == 0) {
                throw new IllegalArgumentException("ambiguous tip state for " +
                        tipDescription(node, taxonId(tree, node)) + " maps to no states in " + context + ": " + state);
            }
            for (int i = 0; i < stateCount; i++) {
                partials[i] /= possibleStateCount;
            }
        } else {
            if (state < 0 || state >= stateCount) {
                throw new IllegalArgumentException("tip state out of range for " +
                        tipDescription(node, taxonId(tree, node)) + " in " + context + ": " + state +
                        " (stateCount=" + stateCount + ")");
            }
            partials[state] = 1.0;
        }

        return partials;
    }

    public static int getObservedState(Tree tree, NodeRef node, PatternList tipPatterns, int stateCount,
                                       String context) {
        int state = getState(tree, node, tipPatterns, context);
        if (tipPatterns.getDataType().isAmbiguousState(state) || state < 0 || state >= stateCount) {
            throw new IllegalArgumentException("tip state must be a single observed state for " +
                    tipDescription(node, taxonId(tree, node)) + " in " + context + ": " + state +
                    " (stateCount=" + stateCount + ")");
        }
        return state;
    }

    private static int getState(Tree tree, NodeRef node, PatternList tipPatterns, String context) {
        int taxonIndex = getTaxonIndex(tree, node, tipPatterns, context);
        return tipPatterns.getPatternState(taxonIndex, 0);
    }

    private static int getTaxonIndex(Tree tree, NodeRef node, PatternList tipPatterns, String context) {
        String taxonId = taxonId(tree, node);
        int taxonIndex = taxonId == null ? -1 : tipPatterns.getTaxonIndex(taxonId);
        if (taxonIndex < 0) {
            throw new IllegalArgumentException(context + " has no entry for " + tipDescription(node, taxonId));
        }
        return taxonIndex;
    }

    private static String taxonId(Tree tree, NodeRef node) {
        Taxon taxon = tree.getNodeTaxon(node);
        return taxon == null ? null : taxon.getId();
    }

    private static String tipDescription(NodeRef node, String taxonId) {
        if (taxonId == null) {
            return "node " + node.getNumber();
        }
        return "node " + node.getNumber() + " (taxon " + taxonId + ")";
    }
}
