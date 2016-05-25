/*
 * NewTreeParameterModel.java
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
import dr.evolution.tree.TreeDoubleTraitProvider;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */
public class NewTreeParameterModel extends AbstractTreeParameterModel<Double> implements TreeDoubleTraitProvider {

    /**
     * This class constructs a tree parameter, and will set the dimension of the parameter
     * to match the appropriate number of nodes if necessary.
     *
     * @param tree        the tree that this parameter corresponds to
     * @param parameter   the parameter to keep in sync with tree topology moves.
     * @param includeRoot tree if the parameter includes a value associated with the root node.
     * @param includeTips
     */
    public NewTreeParameterModel(TreeModel tree, Parameter parameter, boolean includeRoot, boolean includeTips) {
        super(tree, parameter, includeRoot, includeTips);
    }

    @Override
    public int getParameterSize() {
        int treeSize = tree.getNodeCount();
        if (!doesIncludeRoot()) {
            treeSize -= 1;
        }
        if (!doesIncludeTips()) {
            treeSize -= tree.getExternalNodeCount();
        }
        return treeSize;
    }

    public double getNodeDoubleValue(Tree tree, NodeRef node) {
        return getNodeValue(tree, node);
    }

    @Override
    public Double getNodeValue(Tree tree, NodeRef node) {

        assert checkNode(tree, node);

        int nodeNumber = node.getNumber();
        int index = getParameterIndexFromNodeNumber(nodeNumber);
        return parameter.getParameterValue(index);
    }

    @Override
    public void setNodeValue(Tree tree, NodeRef node, Double value) {

        assert checkNode(tree, node);

        int nodeNumber = node.getNumber();
        int index = getParameterIndexFromNodeNumber(nodeNumber);
        parameter.setParameterValue(index, value);
    }

    @Override
    public Double getTrait(Tree tree, NodeRef node) {
        return getNodeValue(tree, node);
    }

    @Override
    public String getTraitString(Tree tree, NodeRef node) {
        return Double.toString(getNodeValue(tree, node));
    }
}
