/*
 * InternalBranchLengths.java
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

package dr.app.treestat.statistics;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;

/**
 * @author Alexei Drummond
 */
public class InternalBranchLengths implements TreeSummaryStatistic {

    private InternalBranchLengths() {
    }

    public int getStatisticDimensions(Tree tree) {
        return tree.getInternalNodeCount() - 1;
    }

    public String getStatisticLabel(Tree tree, int i) {
        return "Branch " + Integer.toString(i + 1);
    }

    public double[] getSummaryStatistic(Tree tree) {

        int internalNodeCount = tree.getInternalNodeCount();
        double[] stats = new double[internalNodeCount - 1];
        int count = 0;
        for (int i = 0; i < internalNodeCount; i++) {
            NodeRef node = tree.getInternalNode(i);
            if (!tree.isRoot(node)) {
                stats[count++] =
                        tree.getNodeHeight(tree.getParent(node)) - tree.getNodeHeight(node);
            }
        }
        //Arrays.sort(stats);

        return stats;
    }

    public void setTaxonList(TaxonList taxonList) {
        throw new UnsupportedOperationException("not implemented in this statistic");
    }

    public void setInteger(int value) {
        throw new UnsupportedOperationException("not implemented in this statistic");
    }

    public void setDouble(double value) {
        throw new UnsupportedOperationException("not implemented in this statistic");
    }

    public void setString(String value) {
        throw new UnsupportedOperationException("not implemented in this statistic");
    }

    public String getSummaryStatisticName() {
        return FACTORY.getSummaryStatisticName();
    }

    public String getSummaryStatisticDescription() {
        return FACTORY.getSummaryStatisticDescription();
    }

    public String getSummaryStatisticReference() {
        return FACTORY.getSummaryStatisticReference();
    }

    public boolean allowsPolytomies() {
        return FACTORY.allowsPolytomies();
    }

    public boolean allowsNonultrametricTrees() {
        return FACTORY.allowsNonultrametricTrees();
    }

    public boolean allowsUnrootedTrees() {
        return FACTORY.allowsUnrootedTrees();
    }

    public Category getCategory() {
        return FACTORY.getCategory();
    }

    public static final Factory FACTORY = new Factory() {

        public TreeSummaryStatistic createStatistic() {
            return new InternalBranchLengths();
        }

        public String getSummaryStatisticName() {
            return "Internal Branch Lengths";
        }

        public String getSummaryStatisticDescription() {

            return "The branch lengths of the internal edges in the tree.";
        }

        public String getSummaryStatisticReference() {

            return "-";
        }

        public boolean allowsPolytomies() {
            return true;
        }

        public boolean allowsNonultrametricTrees() {
            return true;
        }

        public boolean allowsUnrootedTrees() {
            return false;
        }

        public Category getCategory() {
            return Category.GENERAL;
        }
    };
}