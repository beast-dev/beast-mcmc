/*
 * terminalBranchStatistic.java
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Statistic;

/**
 * A statistic that reports the branch lengths of external branches in the tree
 *
 * @author Luiz Carvalho
 */
public class TerminalBranchStatistic extends TreeStatistic {

    public TerminalBranchStatistic(String name, Tree tree) {
        super(name);
        this.tree = tree;
    }

    public void setTree(Tree tree) {
        this.tree = tree;
    }

    public Tree getTree() {
        return tree;
    }

    public int getDimension() {
        return tree.getExternalNodeCount();
    }

    /**
     * @returns external branch lengths
     */
    public double getStatisticValue(int dim) {
        final double[] terminalBranches = new double [tree.getExternalNodeCount()] ;
        int k = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if(tree.isExternal(node)){
            	terminalBranches[k] = tree.getBranchLength(node) ;
            	k = k + 1 ;
            }
        }
        return terminalBranches[dim];
    }
    
    private Tree tree = null;
}