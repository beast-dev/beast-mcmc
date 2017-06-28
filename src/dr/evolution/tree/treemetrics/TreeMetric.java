package dr.evolution.tree.treemetrics;

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;

import java.io.IOException;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface TreeMetric {
    double getMetric(Tree tree1, Tree tree2);

    enum Types {
        ROBINSON_FOULDS("Robinson-Foulds"),
        BRANCH_SCORE("branch score"),
        ROOTED_BRANCH_SCORE("rooted branch score"),
        CLADE_HEIGHT("clade height"),
        KENDALL_COLIJN("Kendall-Colijn path difference"),
        STEEL_PENNY("Steel-Penny path difference");

        Types(String name) {
            this.name = name;
        }

        @Override
        public String toString() {
            return name;
        }

        private final String name;
    }

    class Utils {
        static void checkTreeTaxa(Tree tree1, Tree tree2) {
            //check if taxon lists are in the same order!!
            if (tree1.getExternalNodeCount() != tree2.getExternalNodeCount()) {
                throw new RuntimeException("Different number of taxa in both trees.");
            } else {
                for (int i = 0; i < tree1.getExternalNodeCount(); i++) {
                    if (!tree1.getNodeTaxon(tree1.getExternalNode(i)).getId().equals(tree2.getNodeTaxon(tree2.getExternalNode(i)).getId())) {
                        throw new RuntimeException("Mismatch between taxa in both trees: " + tree1.getNodeTaxon(tree1.getExternalNode(i)).getId() + " vs. " + tree2.getNodeTaxon(tree2.getExternalNode(i)).getId());
                    }
                }
            }
        }
    }
}
