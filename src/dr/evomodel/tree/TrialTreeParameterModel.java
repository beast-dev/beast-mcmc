/*
 * TrialTreeParameterModel.java
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

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;

/**
 * This class maintains a parameter of length equal to the number of nodes in the tree.
 * It can optionally include the root node. If the root node is not included then this
 * class will track tree changes that change the root node number and swap the parameter
 * values so that the parameter values stay with their node when changes to the tree
 * topology occur.
 *
 * @author Alexei Drummond
 * @author Marc A. Suchard
 */
public class TrialTreeParameterModel extends TreeParameterModel {

//    public TrialTreeParameterModel(TreeModel tree, Parameter parameter, boolean includeRoot, boolean includeTips) {
//        this(tree, parameter, includeRoot, includeTips, Intent.NODE);
//    }

    public TrialTreeParameterModel(MutableTreeModel tree, Parameter parameter, boolean includeRoot, boolean includeTips,
                                   Intent intent) {
        super(tree, parameter, includeRoot, intent);
        this.includeTips = includeTips;
    }

    public int getParameterSize() {
        int treeSize = super.getParameterSize();
        if (!includeTips) {
            treeSize -= tree.getExternalNodeCount();
        }
        return treeSize;
    }

    public int getNodeNumberFromParameterIndex(int parameterIndex) {
        int number = super.getNodeNumberFromParameterIndex(parameterIndex);
        if (!includeTips) {
            number += tree.getExternalNodeCount();
        }
        return number;
    }

    public int getParameterIndexFromNodeNumber(int nodeNumber) {
        int number = super.getParameterIndexFromNodeNumber(nodeNumber);
        if (!includeTips) {
            number -= tree.getExternalNodeCount();
        }
        return number;
    }

    public double getNodeValue(Tree tree, NodeRef node) {
        assert (!tree.isExternal(node) && !includeTips) : "tip nodes do not have parameter values!";
        return super.getNodeValue(tree, node);
    }

    public void setNodeValue(Tree tree, NodeRef node, double value) {
        assert (!tree.isExternal(node) && !includeTips) : "tip nodes do not have parameter values!";
        super.setNodeValue(tree, node, value);
    }

    private boolean includeTips = true;
//    protected final Tree tree;

}
