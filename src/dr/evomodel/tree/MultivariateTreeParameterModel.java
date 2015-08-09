/*
 * MultivariateTreeParameterModel.java
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
import dr.evolution.tree.TreeTrait;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */
public class MultivariateTreeParameterModel extends AbstractTreeParameterModel<double[]> {

    /**
     * This class constructs a tree parameter, and will set the dimension of the parameter
     * to match the appropriate number of nodes if necessary.
     *
     * @param tree        the tree that this parameter corresponds to
     * @param parameter   the parameter to keep in sync with tree topology moves.
     * @param includeRoot true if the parameter includes values associated with the root node.
     * @param includeTips true if the parameter includes values associated with the tip nodes.
     * @param dim         dimension of multivariate trait
     */
    public MultivariateTreeParameterModel(TreeModel tree, Parameter parameter,
                                          boolean includeRoot, boolean includeTips, int dim) {
        super(tree, parameter, includeRoot, includeTips);
        this.dim = dim;
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
        return treeSize * dim;
    }

    @Override
    public double[] getNodeValue(Tree tree, NodeRef node) {

        assert checkNode(tree, node);

        int nodeNumber = node.getNumber();
        int index = getParameterIndexFromNodeNumber(nodeNumber);

        double[] result = new double[dim];
        for (int i = 0; i < dim; ++i) {
            result[i] = parameter.getParameterValue(index * dim + i);
        }
        return result;
    }

    @Override
    public void setNodeValue(Tree tree, NodeRef node, double[] value) {

        assert checkNode(tree, node);

        int nodeNumber = node.getNumber();
        int index = getParameterIndexFromNodeNumber(nodeNumber);

        for (int i = 0; i < dim; ++i) {
            parameter.setParameterValue(index * dim + i, value[i]); // TODO Fire change event once
        }
    }

    @Override
    public double[] getTrait(Tree tree, NodeRef node) {
        return getNodeValue(tree, node);
    }

    @Override
    public String getTraitString(Tree tree, NodeRef node) {
        return TreeTrait.DA.formatTrait(getNodeValue(tree, node));
    }

    private final int dim;
}
