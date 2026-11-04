/*
 * SubtreeRootHeightStatistic.java
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

package dr.evomodel.bigfasttree.thorney;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeStatistic;
import dr.inference.model.Statistic;

import java.util.List;

/**
 * This is statistic that returns the subtree root heights from a
 * constrained tree. It can be used to log the tmrca of the clades in the
 * constraints tree.
 */
public class SubtreeRootHeightStatistic extends TreeStatistic {
    private ConstrainedTreeModel tree;
    private int dimension;
    private int subtree;
    private final double mostRecentTipTime;
    private final boolean isBackwards;
    public SubtreeRootHeightStatistic(String name, ConstrainedTreeModel tree, TaxonList taxa,boolean isAbsolute) {
        super(name);

        if(taxa==null){
            subtree=-1;
            dimension=tree.getSubtreeCount();
        }else{
            int[] nodes=new int[taxa.getTaxonCount()];
            for (int i = 0; i < nodes.length; i++) {
                nodes[i]=tree.getExternalNode(tree.getTaxonIndex(taxa.getTaxon(i))).getNumber();
            }
            NodeRef tmrca = TreeUtils.getCommonAncestor(tree,nodes);
            subtree=tree.getSubtreeIndex(tmrca);
            dimension=1;
        }

        if (isAbsolute && Taxon.getMostRecentDate() != null) {
            isBackwards = Taxon.getMostRecentDate().isBackwards();
            mostRecentTipTime = Taxon.getMostRecentDate().getAbsoluteTimeValue();
        } else {
            // give node heights or taxa don't have dates
            mostRecentTipTime = Double.NaN;
            isBackwards = false;
        }
        setTree(tree);
    }


        /**
         * @return the number of dimensions that this statistic has.
         */
    @Override
    public int getDimension() {
        return dimension;
    }

    /**
     * @param dim the dimension to return value of
     * @return the statistic's scalar value in the given dimension
     */
    @Override
    public double getStatisticValue(int dim) {
        Tree stree;
        if(subtree==-1){
            stree = tree.getSubtree(dim);
        }else{
            stree = tree.getSubtree(subtree);
        }

        if (!Double.isNaN(mostRecentTipTime)) {
            if (isBackwards) {
                return mostRecentTipTime + stree.getNodeHeight(stree.getRoot());
            } else {
                return mostRecentTipTime - stree.getNodeHeight(stree.getRoot());
            }
        } else {
            return stree.getNodeHeight(stree.getRoot());
        }
    }

    @Override
    public void setTree(Tree tree) {
        this.tree=(ConstrainedTreeModel) tree;
    }

    @Override
    public Tree getTree() {
        return this.tree;
    }
}