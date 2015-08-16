/*
 * SlicedMultivariateTreeParameterModel.java
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
public class SlicedMultivariateTreeParameterModel extends MultivariateTreeParameterModel implements TreeDoubleTraitProvider {

    /**
     * This class constructs a tree parameter, and will set the dimension of the parameter
     * to match the appropriate number of nodes if necessary.
     *
     * @param tree        the tree that this parameter corresponds to
     * @param parameter   the parameter to keep in sync with tree topology moves.
     * @param includeRoot true if the parameter includes values associated with the root node.
     * @param includeTips true if the parameter includes values associated with the tip nodes.
     * @param dim         dimension of multivariate trait
     * @param slice       specify which entry in multivariate trait to return
     */
    public SlicedMultivariateTreeParameterModel(TreeModel tree, Parameter parameter,
                                                boolean includeRoot, boolean includeTips, int dim, int slice) {
        super(tree, parameter, includeRoot, includeTips, dim);
        this.slice = slice;
    }

    @Override
    public double[] getTrait(Tree tree, NodeRef node) {
        return getNodeValue(tree, node);
    }

    @Override
    public String getTraitString(Tree tree, NodeRef node) {
        return DA.formatTrait(getNodeValue(tree, node));
    }

    public double getNodeDoubleValue(Tree tree, NodeRef node) {
        return getNodeValue(tree, node)[slice];
    }

    private final int slice;
}
