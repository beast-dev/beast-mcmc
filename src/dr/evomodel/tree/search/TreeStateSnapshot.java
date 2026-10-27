/*
 * TreeStateSnapshot.java
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

package dr.evomodel.tree.search;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * A saved state of a time tree and the parameters fitted on it, for the topology search (plan 7 §4.4.2): the best
 * tree, the candidate set, the state before a batch of moves, and reverts after a full re-fit, all of which outlive
 * the one state TreeModel's own store/restore keeps.
 * <p>
 * The snapshot is keyed by node number, as {@code DefaultTreeModel.storeState} is: each node's parent and children,
 * the root number, and the internal node heights, plus the values of the given parameters. This is exact as long as
 * the root keeps its node number, which every move before the root-changing ones does (rooted NNIs, root-child NNIs
 * and SPRs with a non-root parent keep the root; TreeModel's own store and restore swap node objects but keep the
 * numbers): rates and node-height proxies index nodes by number, skipping the root's, so their indices keep their
 * meaning. A restore under a different root number would permute them, so it fails instead.
 * <p>
 * Tip heights are never stored or written: they come from the dates.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class TreeStateSnapshot {

    private static final double TOLERANCE = 1E-6;

    private final TreeModel tree;
    private final int rootNumber;
    private final int[][] children; // by node number; empty for tips
    private final double[] heights; // by node number; NaN for tips
    private final List<Parameter> parameters;
    private final double[][] values;
    private final Likelihood objective;
    private final double logDensity;

    /**
     * Saves the current state.
     *
     * @param parameters the parameters to save with the tree (rates, model parameters, or a maximizer's parameter);
     *                   restored by value, so each must keep its meaning while the root is unchanged
     * @param objective  if not null, evaluated now and checked against the restored state
     */
    public TreeStateSnapshot(TreeModel tree, List<? extends Parameter> parameters, Likelihood objective) {
        this.tree = tree;
        this.rootNumber = tree.getRoot().getNumber();
        final int nodeCount = tree.getNodeCount();
        this.children = new int[nodeCount][];
        this.heights = new double[nodeCount];
        for (int i = 0; i < nodeCount; ++i) {
            final NodeRef node = tree.getNode(i);
            children[i] = new int[tree.getChildCount(node)];
            for (int c = 0; c < children[i].length; ++c) {
                children[i][c] = tree.getChild(node, c).getNumber();
            }
            heights[i] = tree.isExternal(node) ? Double.NaN : tree.getNodeHeight(node);
        }
        this.parameters = parameters == null ? Collections.<Parameter>emptyList() :
                Collections.unmodifiableList(new ArrayList<Parameter>(parameters));
        this.values = new double[this.parameters.size()][];
        for (int p = 0; p < values.length; ++p) {
            values[p] = this.parameters.get(p).getParameterValues();
        }
        this.objective = objective;
        this.logDensity = objective != null ? objective.getLogLikelihood() : Double.NaN;
    }

    /**
     * Puts the tree and the parameters back as they were saved, then (with an objective) evaluates it from scratch
     * and checks that it equals the saved value.
     *
     * @throws IllegalStateException when the root node number has changed since the save, or the objective differs
     */
    public void restore() {
        final int currentRoot = tree.getRoot().getNumber();
        if (currentRoot != rootNumber) {
            throw new IllegalStateException("The root is node " + currentRoot + ", not node " + rootNumber +
                    " as at the save: a snapshot keyed by node number cannot be restored across a root change");
        }

        tree.beginTreeEdit();
        for (int i = 0; i < children.length; ++i) {
            final NodeRef node = tree.getNode(i);
            for (int c = tree.getChildCount(node) - 1; c >= 0; --c) {
                tree.removeChild(node, tree.getChild(node, c));
            }
        }
        for (int i = 0; i < children.length; ++i) {
            final NodeRef node = tree.getNode(i);
            for (int child : children[i]) {
                tree.addChild(node, tree.getNode(child));
            }
            if (!Double.isNaN(heights[i])) {
                tree.setNodeHeightQuietly(node, heights[i]);
            }
        }
        tree.endTreeEdit();
        tree.pushTreeChangedEvent();

        for (int p = 0; p < values.length; ++p) {
            final Parameter parameter = parameters.get(p);
            for (int i = 0; i < values[p].length; ++i) {
                if (parameter.getParameterValue(i) != values[p][i]) {
                    parameter.setParameterValue(i, values[p][i]);
                }
            }
        }

        if (objective != null) {
            objective.makeDirty();
            final double restored = objective.getLogLikelihood();
            if (!(Math.abs(restored - logDensity) <= TOLERANCE * Math.max(1.0, Math.abs(logDensity)))) {
                throw new IllegalStateException("The restored state has log density " + restored + ", not " +
                        logDensity + " as at the save");
            }
        }
    }

    public int getRootNumber() {
        return rootNumber;
    }

    /**
     * @return the objective's value at the save (NaN without an objective)
     */
    public double getLogDensity() {
        return logDensity;
    }

    public TreeModel getTree() {
        return tree;
    }
}
