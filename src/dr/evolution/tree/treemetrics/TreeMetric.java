/*
 * TreeMetric.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.evolution.tree.treemetrics;

import dr.evolution.tree.Tree;

/**
 * @author Andrew Rambaut
 */
public interface TreeMetric {

    double getMetric(Tree tree1, Tree tree2);

    Type getType();

    enum Type {
        ROBINSON_FOULDS("Robinson-Foulds", "rf"),
        BRANCH_SCORE("branch score", "branchscore"),
        ROOTED_BRANCH_SCORE("rooted branch score", "branch"),
        CLADE_HEIGHT("clade height", "clade"),
        KENDALL_COLIJN("Kendall-Colijn path difference", "kc"),
        STEEL_PENNY("Steel-Penny path difference", "sp");

        Type(String name, String stortName) {
            this.name = name;
            this.stortName = stortName;
        }

        public String getName() {
            return name;
        }

        public String getShortName() {
            return stortName;
        }

        @Override
        public String toString() {
            return getName();
        }

        private final String name;
        private final String stortName;
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
