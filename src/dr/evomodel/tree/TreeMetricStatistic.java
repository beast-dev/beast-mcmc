/*
 * TreeMetricStatistic.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.evomodel.tree;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.tree.treemetrics.TreeMetric;

/**
 * A statistic that returns the distance between two trees.
 *
 * @author Andrew Rambaut
 */
public class TreeMetricStatistic extends TreeStatistic {

    /**
     * Constructor which creates statistic that just says whether two trees have the same topology
     * @param name
     * @param referenceTree
     */
    public TreeMetricStatistic(String name, Tree referenceTree, Tree targetTree) {
        this(name, referenceTree, targetTree, null);
    }

    public TreeMetricStatistic(String name, Tree referenceTree, Tree targetTree, TreeMetric treeMetric) {
        super(name);

        this.referenceTree = referenceTree;
        this.targetTree = targetTree;
        this.treeMetric = treeMetric;
        this.focalNewick = TreeUtils.uniqueNewick(referenceTree, referenceTree.getRoot());
    }

    public void setTree(Tree tree) {
        this.targetTree = tree;
    }

    public Tree getTree() {
        return targetTree;
    }

    public int getDimension() {
        return 1;
    }

    /**
     * @return value.
     */
    public double getStatisticValue(int dim) {

        if (treeMetric == null) {
            // simply return if the two trees have the same topology
            return compareTreesByTopology();
        }

        return treeMetric.getMetric(referenceTree, targetTree);
    }

    private double compareTreesByTopology() {
        final String targetNewick = TreeUtils.uniqueNewick(targetTree, targetTree.getRoot());
        return targetNewick.equals(focalNewick) ? 1.0 : 0.0;
    }

    private Tree targetTree;
    private final Tree referenceTree;

    private final String focalNewick;

    private final TreeMetric treeMetric;
}
