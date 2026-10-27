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
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;

/**
 * A saved state of a time tree and the parameters fitted on it, for the topology search (plan 7 §4.4.2): the best
 * tree, the candidate set, the state before a batch of moves, and reverts after a full re-fit, all of which outlive
 * the one state TreeModel's own store/restore keeps.
 * <p>
 * The snapshot is keyed by node number, as {@code DefaultTreeModel.storeState} is: each node's parent and children,
 * the root number, and the internal node heights, plus the values of the given parameters. While the root keeps its
 * node number this is exact: rates and node-height proxies index nodes by number, skipping the root's, so their
 * indices keep their meaning, and the restore puts every node number back as it was.
 * <p>
 * A root-changing move (plan 7 §5.4) gives the root another node, which shifts those indices. Restoring across a
 * root change is therefore by clade (M6): the saved topology is rebuilt on the existing node objects, with the
 * current root node kept as the root, each other saved clade on the node that holds that clade now if there is one,
 * and the rest on the remaining internal nodes, lowest number first. Heights then follow the clades, as do the values
 * of branch-indexed parameters (those of the objective's {@link TreeParameterModel}s on this tree), which are set
 * at their new indices; other parameters are restored by value, and node-height proxies through the heights.
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

    // for a restore by clade: the leaf parameters restored by value, and the branch-indexed ones by node number
    private final List<Parameter> byValue = new ArrayList<Parameter>();
    private final double[][] byValueValues;
    private final List<TreeParameterModel> branchModels = new ArrayList<TreeParameterModel>();
    private final List<Parameter> branchParameters = new ArrayList<Parameter>();
    private final double[][] branchValues; // by node number; NaN where the node has no value (the root)

    /**
     * Saves the current state.
     *
     * @param parameters the parameters to save with the tree (rates, model parameters, or a maximizer's parameter)
     * @param objective  if not null, evaluated now and checked against the restored state; its models tell which
     *                   parameters are indexed by branch, which a restore across a root change needs
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

        // the leaves, each once, split into branch-indexed parameters and the rest
        final List<Parameter> leaves = new ArrayList<Parameter>();
        final Map<Parameter, Boolean> seen = new IdentityHashMap<Parameter, Boolean>();
        for (Parameter parameter : this.parameters) {
            final List<Parameter> flat = new ArrayList<Parameter>();
            TrialEvaluator.flatten(parameter, flat);
            for (Parameter leaf : flat) {
                if (seen.put(leaf, Boolean.TRUE) == null) {
                    leaves.add(leaf);
                }
            }
        }
        final Map<Parameter, TreeParameterModel> branchIndexed = objective != null ?
                branchIndexed(tree, objective.getModel()) : Collections.<Parameter, TreeParameterModel>emptyMap();
        for (Parameter leaf : leaves) {
            final TreeParameterModel model = branchIndexed.get(leaf);
            if (model != null) {
                branchModels.add(model);
                branchParameters.add(leaf);
            } else {
                byValue.add(leaf);
            }
        }
        this.byValueValues = new double[byValue.size()][];
        for (int p = 0; p < byValueValues.length; ++p) {
            byValueValues[p] = byValue.get(p).getParameterValues();
        }
        this.branchValues = new double[branchParameters.size()][nodeCount];
        for (int p = 0; p < branchValues.length; ++p) {
            final TreeParameterModel model = branchModels.get(p);
            final Parameter parameter = branchParameters.get(p);
            final boolean includeRoot = model.getParameterSize() == nodeCount;
            for (int i = 0; i < nodeCount; ++i) {
                branchValues[p][i] = !includeRoot && i == rootNumber ? Double.NaN :
                        parameter.getParameterValue(model.getParameterIndexFromNodeNumber(i));
            }
        }
    }

    // a topology with heights and no parameters, to adopt by clade
    private TreeStateSnapshot(TreeModel tree, int rootNumber, int[][] children, double[] heights) {
        this.tree = tree;
        this.rootNumber = rootNumber;
        this.children = children;
        this.heights = heights;
        this.parameters = Collections.emptyList();
        this.values = new double[0][];
        this.objective = null;
        this.logDensity = Double.NaN;
        this.byValueValues = new double[0][];
        this.branchValues = new double[0][];
    }

    /**
     * A topology to adopt by clade (an NJ start tree, plan 7 §7.1): node numbers index the arrays, tips by their
     * numbers in the tree; heights are those of the internal nodes (tips are ignored). Restoring it leaves every
     * parameter as it is, so per-branch rates must be set afterwards.
     *
     * @param children for each node number, its children's numbers (empty for tips)
     */
    public static TreeStateSnapshot ofTopology(TreeModel tree, int rootNumber, int[][] children, double[] heights) {
        if (children.length != tree.getNodeCount() || heights.length != children.length) {
            throw new IllegalArgumentException("The topology has " + children.length + " nodes, the tree " +
                    tree.getNodeCount());
        }
        return new TreeStateSnapshot(tree, rootNumber, children.clone(), heights.clone());
    }

    // the parameters of the TreeParameterModels on this tree among the models, with their models
    private static Map<Parameter, TreeParameterModel> branchIndexed(TreeModel tree, Model root) {
        final Map<Parameter, TreeParameterModel> result = new IdentityHashMap<Parameter, TreeParameterModel>();
        if (root == null) {
            return result;
        }
        for (Model model : TrialEvaluator.models(root)) {
            if (model instanceof TreeParameterModel && ((TreeParameterModel) model).getTreeModel() == tree) {
                final TreeParameterModel treeParameter = (TreeParameterModel) model;
                for (int v = 0; v < treeParameter.getVariableCount(); ++v) {
                    if (treeParameter.getVariable(v) instanceof Parameter
                            && treeParameter.getVariable(v).getSize() == treeParameter.getParameterSize()) {
                        result.put((Parameter) treeParameter.getVariable(v), treeParameter);
                    }
                }
            }
        }
        return result;
    }

    /**
     * Puts the tree and the parameters back as they were saved, then (with an objective) evaluates it from scratch
     * and checks that it equals the saved value. Under the root number of the save, node numbers are restored as
     * they were; under another, the restore is by clade and keeps the current root node.
     *
     * @throws IllegalStateException when the objective differs
     */
    public void restore() {
        if (tree.getRoot().getNumber() != rootNumber) {
            restoreByClade(tree.getRoot().getNumber());
            return;
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
            setValues(parameters.get(p), values[p]);
        }
        check();
    }

    /**
     * Puts the saved topology, heights and parameters back by clade (plan 7 §4.4.2), with the given internal node as
     * the root: the current root to keep it, as restore does across a root change, or another node to move the root
     * there. Each other saved clade goes to the node that holds that clade now, if there is one, and the rest to the
     * remaining internal nodes in the saved post-order, lowest number first, so the result is deterministic.
     *
     * @throws IllegalStateException when the objective differs, or without an objective (which tells the
     *                               branch-indexed parameters apart)
     */
    public void restoreByClade(int newRootNumber) {
        if (objective == null && values.length > 0) {
            throw new IllegalStateException("A restore by clade needs the objective to tell branch-indexed " +
                    "parameters apart");
        }
        final int nodeCount = children.length;
        final NodeRef newRoot = tree.getNode(newRootNumber);
        if (tree.isExternal(newRoot)) {
            throw new IllegalArgumentException("Node " + newRootNumber + " is a tip");
        }

        // the saved nodes in post-order, with their clades
        final List<Integer> postOrder = new ArrayList<Integer>();
        final BitSet[] saved = new BitSet[nodeCount];
        postOrder(rootNumber, postOrder, saved);
        final BitSet[] current = TopologyConstraints.clades(tree);

        // saved node number -> node number it is rebuilt on
        final int[] on = new int[nodeCount];
        final boolean[] used = new boolean[nodeCount];
        for (int i = 0; i < nodeCount; ++i) {
            on[i] = -1;
            if (children[i].length == 0) {
                on[i] = i; // tips keep their numbers
                used[i] = true;
            }
        }
        on[rootNumber] = newRootNumber;
        used[newRootNumber] = true;
        final Map<BitSet, Integer> holders = new HashMap<BitSet, Integer>();
        for (int m = 0; m < nodeCount; ++m) {
            if (!tree.isExternal(tree.getNode(m)) && m != newRootNumber) {
                holders.put(current[m], m);
            }
        }
        for (int s : postOrder) {
            if (on[s] < 0) {
                final Integer holder = holders.get(saved[s]);
                if (holder != null && !used[holder]) {
                    on[s] = holder;
                    used[holder] = true;
                }
            }
        }
        int free = 0;
        for (int s : postOrder) {
            if (on[s] < 0) {
                while (used[free] || tree.isExternal(tree.getNode(free))) {
                    ++free;
                }
                on[s] = free;
                used[free] = true;
            }
        }

        tree.beginTreeEdit();
        for (int m = 0; m < nodeCount; ++m) {
            final NodeRef node = tree.getNode(m);
            for (int c = tree.getChildCount(node) - 1; c >= 0; --c) {
                tree.removeChild(node, tree.getChild(node, c));
            }
        }
        for (int s : postOrder) {
            final NodeRef node = tree.getNode(on[s]);
            for (int child : children[s]) {
                tree.addChild(node, tree.getNode(on[child]));
            }
            if (!Double.isNaN(heights[s])) {
                tree.setNodeHeightQuietly(node, heights[s]);
            }
        }
        if (tree.getRoot() != newRoot) {
            tree.setRoot(newRoot);
        }
        tree.endTreeEdit(); // a new root object: TreeModel swaps its parameter objects, rate models shift indices
        tree.pushTreeChangedEvent();

        // every branch value at its new index: none may stay where it was
        for (int p = 0; p < branchValues.length; ++p) {
            final TreeParameterModel model = branchModels.get(p);
            final Parameter parameter = branchParameters.get(p);
            final double[] target = parameter.getParameterValues();
            for (int s = 0; s < nodeCount; ++s) {
                if (!Double.isNaN(branchValues[p][s])) {
                    target[model.getParameterIndexFromNodeNumber(on[s])] = branchValues[p][s];
                }
            }
            setValues(parameter, target);
        }
        for (int p = 0; p < byValueValues.length; ++p) {
            setValues(byValue.get(p), byValueValues[p]);
        }
        check();
    }

    private void postOrder(int node, List<Integer> order, BitSet[] clades) {
        final BitSet bits = new BitSet(tree.getExternalNodeCount());
        if (children[node].length == 0) {
            bits.set(node);
        } else {
            for (int child : children[node]) {
                postOrder(child, order, clades);
                bits.or(clades[child]);
            }
            order.add(node);
        }
        clades[node] = bits;
    }

    private static void setValues(Parameter parameter, double[] values) {
        for (int i = 0; i < values.length; ++i) {
            if (parameter.getParameterValue(i) != values[i]) {
                parameter.setParameterValue(i, values[i]);
            }
        }
    }

    private void check() {
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
