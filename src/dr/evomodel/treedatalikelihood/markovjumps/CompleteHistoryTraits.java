/*
 * CompleteHistoryTraits.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.treedatalikelihood.markovjumps;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treelikelihood.MarkovJumpsBeagleTreeLikelihood;

/**
 * The untagged traits {@code CompleteHistoryLogger} looks up by name, built from an engine's own single-pattern
 * traits: {@code history}, one branch's complete history as {site, time, from, to} events, and
 * {@code allTransitions}, the number of jumps summed over every branch of the tree. They are not loggable, so tree
 * files carry only the engine's tagged traits; CompleteHistoryLogger reads them directly.
 */
public final class CompleteHistoryTraits {

    private CompleteHistoryTraits() {
    }

    /**
     * @param historyTrait per-branch history strings, one per pattern (only the first is used)
     * @param countTrait   per-branch jump counts, one per pattern (only the first is used)
     */
    public static TreeTrait[] createLoggerTraits(final TreeTrait historyTrait, final TreeTrait countTrait) {
        if (historyTrait == null || countTrait == null) {
            throw new IllegalArgumentException("history and count traits are both required");
        }

        TreeTrait<String> history = new TreeTrait.S() {
            public String getTraitName() {
                return MarkovJumpsBeagleTreeLikelihood.HISTORY;
            }

            public Intent getIntent() {
                return Intent.BRANCH;
            }

            public String getTrait(Tree tree, NodeRef node) {
                String[] values = (String[]) historyTrait.getTrait(tree, node);
                return values == null ? null : values[0];
            }

            public boolean getLoggable() {
                return false;
            }
        };

        TreeTrait<Double> total = new TreeTrait.D() {
            public String getTraitName() {
                return MarkovJumpsBeagleTreeLikelihood.TOTAL_COUNTS;
            }

            public Intent getIntent() {
                return Intent.WHOLE_TREE;
            }

            public Double getTrait(Tree tree, NodeRef ignored) {
                double sum = 0.0;
                for (int i = 0; i < tree.getNodeCount(); i++) {
                    NodeRef node = tree.getNode(i);
                    if (!tree.isRoot(node)) {
                        sum += ((double[]) countTrait.getTrait(tree, node))[0];
                    }
                }
                return sum;
            }

            public boolean getLoggable() {
                return false;
            }
        };

        return new TreeTrait[]{history, total};
    }
}
