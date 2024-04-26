/*
 * SubtreeLeapOperator.java
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

package dr.evomodel.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.operators.SubtreeLeapOperatorParser;
import dr.evomodelxml.operators.TipLeapOperatorParser;
import dr.inference.distribution.CauchyDistribution;
import dr.inference.model.Bounds;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.RandomWalkOperator;
import dr.inference.operators.Scalable;
import dr.inferencexml.operators.ScaleOperatorParser;
import dr.math.MathUtils;
import dr.math.distributions.Distribution;
import dr.util.Transform;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

/**
 * Implements moves that changes node heights - Uniform version
 *
 * Created to allow these operations on TreeModels that can't expose node heights as parameters for efficiency
 * reasons such as BigFastTreeModel.
 *
 * @author Andrew Rambaut
 * @version $Id$
 */
public class UniformNodeHeightOperator extends AbstractTreeOperator {

    private final TreeModel tree;

    /**
     * Constructor
     *
     * @param tree   the tree
     * @param weight the weight
     */
    public UniformNodeHeightOperator(TreeModel tree, double weight) {
        super();

        this.tree = tree;
        setWeight(weight);
    }

    /**
     * Do a subtree leap move.
     *
     * @return the log-transformed hastings ratio
     */
    public double doOperation() {
        final NodeRef root = tree.getRoot();

        NodeRef node;
        // Pick a node (but not the root)
        do {
            // choose a internal node avoiding root
            node = tree.getInternalNode(MathUtils.nextInt(tree.getInternalNodeCount()));
        } while (node == root);

        final double upperHeight = tree.getNodeHeight(tree.getParent(node));
        final double lowerHeight = Math.max(
                tree.getNodeHeight(tree.getChild(node, 0)),
                tree.getNodeHeight(tree.getChild(node, 1)));

        final double oldHeight = tree.getNodeHeight(node);

        return doUniform(node, oldHeight, upperHeight, lowerHeight);
    }

    private double doUniform(NodeRef node, double oldValue, double upper, double lower) {
        tree.setNodeHeight(node, (MathUtils.nextDouble() * (upper - lower)) + lower);
        return 0.0;
    }

    public String getOperatorName() {
        return "uniform(" + tree.getId() + " internal nodes)";
    }
}
