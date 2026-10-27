/*
 * TimeTreeInitializer.java
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

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.distance.DistanceMatrix;
import dr.evolution.distance.F84DistanceMatrix;
import dr.evolution.distance.JukesCantorDistanceMatrix;
import dr.evolution.tree.NeighborJoiningTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A neighbour-joining start for the topology search (plan 7 §7.1, optional): the NJ tree of the alignment's
 * distances, rooted, dated and adopted by clade, keeping the root's node.
 * <ol>
 * <li><b>Rooting.</b> With tip dates, the root position (an edge and a point on it) whose root-to-tip distances fit
 * a line in the dates best, with a positive slope, the rate r̂ (root-to-tip regression; the residual sum of squares
 * is quadratic in the position on an edge, so each edge's best point is exact). Without dates, or without a positive
 * slope, the position with the least variance of root-to-tip distances, and r̂ the mean root-to-tip distance over
 * the current root height.</li>
 * <li><b>Heights.</b> Tips keep their heights (the dates); in post-order h(v) = max over children c of
 * h(c) + max(b_c / r̂, δ), which is always a valid time tree.</li>
 * <li><b>Rates.</b> Per-branch rates b_c / (h(p) - h(c)), within [r̂ / 100, 100 r̂]; a strict clock's rate is left
 * to the first re-fit.</li>
 * <li><b>Constraints.</b> An NJ tree that breaks one is not adopted (the plan's repair by SPR is not built).</li>
 * </ol>
 * The distances are F84 for nucleotides, Jukes-Cantor otherwise, from the patterns of the first tree likelihood on
 * the tree.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class TimeTreeInitializer {

    private static final double RATE_RANGE = 100.0;

    private final TreeModel tree;
    private final ArbitraryBranchRates rates;

    // the NJ tree as an undirected graph over its nodes; tips mapped to the tree's tip numbers
    private int nodeCount;
    private List<int[]> edges; // {a, b}
    private double[] lengths;
    private List<List<Integer>> adjacent; // edge indices by node
    private int[] tipNumber; // NJ node -> tree tip number, -1 for internal nodes

    private double rate = Double.NaN;
    private boolean dated = false;

    public TimeTreeInitializer(TreeModel tree, ArbitraryBranchRates rates) {
        this.tree = tree;
        this.rates = rates;
    }

    /**
     * @return the patterns of the first BEAGLE tree likelihood on the tree among the objective's models, or null
     */
    public static PatternList patterns(TreeModel tree, Likelihood objective) {
        for (Model model : TrialEvaluator.models(objective.getModel())) {
            if (model instanceof TreeDataLikelihood && ((TreeDataLikelihood) model).getTree() == tree
                    && ((TreeDataLikelihood) model).getDataLikelihoodDelegate() instanceof BeagleDataLikelihoodDelegate) {
                return ((BeagleDataLikelihoodDelegate) ((TreeDataLikelihood) model).getDataLikelihoodDelegate())
                        .getPatternList();
            }
        }
        return null;
    }

    /**
     * Builds, roots and dates the NJ tree of the patterns and adopts it.
     *
     * @return false (and the tree as it was) when the NJ tree breaks a constraint
     */
    public boolean adopt(PatternList patterns, TopologyConstraints constraints) {
        final DistanceMatrix distances = patterns.getDataType() instanceof Nucleotides ?
                new F84DistanceMatrix(patterns) : new JukesCantorDistanceMatrix(patterns);
        return adopt(new NeighborJoiningTree(distances), constraints);
    }

    /**
     * Roots and dates an unrooted (or arbitrarily rooted) tree with branch lengths in substitutions and adopts it.
     */
    public boolean adopt(Tree donor, TopologyConstraints constraints) {
        graph(donor);
        final double[] position = root();
        final int[][] children = new int[tree.getNodeCount()][];
        final double[] heights = new double[tree.getNodeCount()];
        final double[] branch = new double[tree.getNodeCount()]; // substitutions above each node
        Arrays.fill(heights, Double.NaN);
        final int rootNumber = build(position, children, heights, branch);

        final TreeStateSnapshot topology = TreeStateSnapshot.ofTopology(tree, rootNumber, children, heights);
        final TreeStateSnapshot before = new TreeStateSnapshot(tree, null, null);
        topology.restoreByClade(tree.getRoot().getNumber());
        if (constraints != null && !constraints.areSatisfied(tree)) {
            before.restoreByClade(tree.getRoot().getNumber());
            return false;
        }
        if (rates != null) {
            setRates(children, heights, branch, rootNumber);
        }
        return true;
    }

    /**
     * @return the rate r̂ of the last adoption: the regression slope with dates, else the mean root-to-tip distance
     * over the root height
     */
    public double getRate() {
        return rate;
    }

    /**
     * @return whether the last adoption was rooted by the dates (root-to-tip regression)
     */
    public boolean isDated() {
        return dated;
    }

    private void graph(Tree donor) {
        nodeCount = donor.getNodeCount();
        edges = new ArrayList<int[]>();
        adjacent = new ArrayList<List<Integer>>();
        for (int n = 0; n < nodeCount; ++n) {
            adjacent.add(new ArrayList<Integer>());
        }
        final List<Double> length = new ArrayList<Double>();
        for (int n = 0; n < nodeCount; ++n) {
            final NodeRef node = donor.getNode(n);
            if (!donor.isRoot(node)) {
                edges.add(new int[]{n, donor.getParent(node).getNumber()});
                length.add(Math.max(0.0, donor.getNodeHeight(donor.getParent(node)) - donor.getNodeHeight(node)));
            }
        }
        // the donor's root has two children: merge its two edges into one, so the root is not a node of the graph
        final NodeRef donorRoot = donor.getRoot();
        if (donor.getChildCount(donorRoot) == 2) {
            final int r = donorRoot.getNumber();
            final int a = donor.getChild(donorRoot, 0).getNumber();
            final int b = donor.getChild(donorRoot, 1).getNumber();
            double merged = 0.0;
            for (int e = edges.size() - 1; e >= 0; --e) {
                if (edges.get(e)[1] == r) {
                    merged += length.get(e);
                    edges.remove(e);
                    length.remove(e);
                }
            }
            edges.add(new int[]{a, b});
            length.add(merged);
        }
        lengths = new double[edges.size()];
        for (int e = 0; e < lengths.length; ++e) {
            lengths[e] = length.get(e);
            adjacent.get(edges.get(e)[0]).add(e);
            adjacent.get(edges.get(e)[1]).add(e);
        }
        tipNumber = new int[nodeCount];
        Arrays.fill(tipNumber, -1);
        for (int n = 0; n < nodeCount; ++n) {
            final NodeRef node = donor.getNode(n);
            if (donor.isExternal(node)) {
                final String id = donor.getNodeTaxon(node).getId();
                for (int t = 0; t < tree.getExternalNodeCount(); ++t) {
                    if (tree.getNodeTaxon(tree.getExternalNode(t)).getId().equals(id)) {
                        tipNumber[n] = tree.getExternalNode(t).getNumber();
                    }
                }
                if (tipNumber[n] < 0) {
                    throw new IllegalArgumentException("Taxon " + id + " is not in the tree");
                }
            }
        }
        if (edges.size() != 2 * tree.getExternalNodeCount() - 3) {
            throw new IllegalArgumentException("The start tree is not binary");
        }
    }

    // the distances from a node to every NJ tip, by tree tip number
    private double[] tipDistances(int from) {
        final double[] result = new double[tree.getExternalNodeCount()];
        final double[] distance = new double[nodeCount];
        Arrays.fill(distance, -1.0);
        distance[from] = 0.0;
        final List<Integer> stack = new ArrayList<Integer>();
        stack.add(from);
        while (!stack.isEmpty()) {
            final int x = stack.remove(stack.size() - 1);
            if (tipNumber[x] >= 0) {
                result[tipNumber[x]] = distance[x];
            }
            for (int e : adjacent.get(x)) {
                final int y = edges.get(e)[0] == x ? edges.get(e)[1] : edges.get(e)[0];
                if (distance[y] < 0) {
                    distance[y] = distance[x] + lengths[e];
                    stack.add(y);
                }
            }
        }
        return result;
    }

    // the best root position {edge, distance from its first node}; sets rate and dated
    private double[] root() {
        final int tips = tree.getExternalNodeCount();
        final double[] dates = new double[tips];
        double meanDate = 0.0;
        for (int t = 0; t < tips; ++t) {
            dates[t] = -tree.getNodeHeight(tree.getExternalNode(t));
            meanDate += dates[t] / tips;
        }
        double sxx = 0.0;
        for (int t = 0; t < tips; ++t) {
            sxx += (dates[t] - meanDate) * (dates[t] - meanDate);
        }
        double[] best = null;
        double bestValue = Double.POSITIVE_INFINITY;
        double bestSlope = Double.NaN;
        double[] fallback = null;
        double fallbackValue = Double.POSITIVE_INFINITY;
        for (int e = 0; e < edges.size(); ++e) {
            // root-to-tip distance at x from a: c_t + s_t x, s_t = +1 on a's side, -1 on b's side
            final double[] fromA = tipDistances(edges.get(e)[0]);
            final double[] fromB = tipDistances(edges.get(e)[1]);
            final double[] c = new double[tips];
            final double[] s = new double[tips];
            for (int t = 0; t < tips; ++t) {
                // on a's side the path to b runs through a
                if (Math.abs(fromA[t] + lengths[e] - fromB[t]) <= Math.abs(fromB[t] + lengths[e] - fromA[t])) {
                    c[t] = fromA[t];
                    s[t] = 1.0;
                } else {
                    c[t] = fromB[t] + lengths[e];
                    s[t] = -1.0;
                }
            }
            // variance of the distances, and the regression's residual sum of squares, as quadratics in x
            double mc = 0.0, ms = 0.0;
            for (int t = 0; t < tips; ++t) {
                mc += c[t] / tips;
                ms += s[t] / tips;
            }
            double scc = 0.0, scs = 0.0, sss = 0.0, sxc = 0.0, sxs = 0.0;
            for (int t = 0; t < tips; ++t) {
                final double dc = c[t] - mc, ds = s[t] - ms, dx = dates[t] - meanDate;
                scc += dc * dc;
                scs += dc * ds;
                sss += ds * ds;
                sxc += dx * dc;
                sxs += dx * ds;
            }
            // Syy(x) = scc + 2 scs x + sss x^2; Sxy(x) = sxc + sxs x; RSS(x) = Syy - Sxy^2 / Sxx
            final double qa = sxx > 0 ? sss - sxs * sxs / sxx : sss;
            final double qb = sxx > 0 ? 2 * scs - 2 * sxc * sxs / sxx : 2 * scs;
            final double qc = sxx > 0 ? scc - sxc * sxc / sxx : scc;
            for (double x : candidates(qa, qb, lengths[e])) {
                final double value = qc + qb * x + qa * x * x;
                final double slope = sxx > 0 ? (sxc + sxs * x) / sxx : Double.NaN;
                if (slope > 0 && value < bestValue) {
                    bestValue = value;
                    best = new double[]{e, x};
                    bestSlope = slope;
                }
                final double variance = scc + 2 * scs * x + sss * x * x;
                if (variance < fallbackValue) {
                    fallbackValue = variance;
                    fallback = new double[]{e, x};
                }
            }
        }
        dated = best != null;
        if (dated) {
            rate = bestSlope;
            return best;
        }
        rate = Double.NaN; // set from the root height once the root is placed
        return fallback;
    }

    // the ends of [0, length] and the quadratic's stationary point inside it
    private static double[] candidates(double a, double b, double length) {
        if (a > 0) {
            final double x = -b / (2 * a);
            if (x > 0 && x < length) {
                return new double[]{0.0, length, x};
            }
        }
        return new double[]{0.0, length};
    }

    // the rooted tree in the tree's node numbers: tips keep theirs, internal nodes take the free numbers; returns the
    // root's number
    private int build(double[] position, int[][] children, double[] heights, double[] branch) {
        final int e = (int) position[0];
        final double x = position[1];
        final int a = edges.get(e)[0];
        final int b = edges.get(e)[1];

        // number the NJ internal nodes, then a new root on edge e
        final int[] number = new int[nodeCount];
        int next = tree.getExternalNodeCount();
        for (int n = 0; n < nodeCount; ++n) {
            number[n] = tipNumber[n] >= 0 ? tipNumber[n] : -1;
        }
        for (int n = 0; n < nodeCount; ++n) {
            if (number[n] < 0 && !adjacent.get(n).isEmpty()) {
                number[n] = next++;
            }
        }
        final int root = next;
        if (root != tree.getNodeCount() - 1) {
            throw new IllegalStateException("The start tree has " + (root + 1) + " nodes, the tree " +
                    tree.getNodeCount());
        }
        Arrays.fill(children, new int[0]);
        children[root] = new int[]{number[a], number[b]};
        branch[number[a]] = x;
        branch[number[b]] = lengths[e] - x;
        orient(a, e, number, children, branch);
        orient(b, e, number, children, branch);

        // the rate: from the regression, else the mean root-to-tip distance over the current root height
        final double[] depth = new double[tree.getNodeCount()];
        fillDepth(root, 0.0, children, branch, depth);
        if (!dated) {
            double mean = 0.0;
            double tipMean = 0.0;
            final int tips = tree.getExternalNodeCount();
            for (int t = 0; t < tips; ++t) {
                mean += depth[tree.getExternalNode(t).getNumber()] / tips;
                tipMean += tree.getNodeHeight(tree.getExternalNode(t)) / tips;
            }
            final double age = tree.getNodeHeight(tree.getRoot()) - tipMean;
            rate = mean > 0 && age > 0 ? mean / age : 1.0;
        }

        // heights in post-order: h(v) = max_c h(c) + max(b_c / r, delta)
        double totalDepth = 0.0;
        for (double d : depth) {
            totalDepth = Math.max(totalDepth, d);
        }
        final double delta = Math.max(1E-12, 1E-3 * totalDepth / rate / tree.getExternalNodeCount());
        heightsBelow(root, children, branch, heights, delta);
        return root;
    }

    private void orient(int node, int from, int[] number, int[][] children, double[] branch) {
        final List<Integer> below = new ArrayList<Integer>();
        for (int e : adjacent.get(node)) {
            if (e != from) {
                final int y = edges.get(e)[0] == node ? edges.get(e)[1] : edges.get(e)[0];
                below.add(number[y]);
                branch[number[y]] = lengths[e];
                orient(y, e, number, children, branch);
            }
        }
        children[number[node]] = new int[below.size()];
        for (int k = 0; k < below.size(); ++k) {
            children[number[node]][k] = below.get(k);
        }
    }

    private static void fillDepth(int node, double depth, int[][] children, double[] branch, double[] result) {
        result[node] = depth;
        for (int child : children[node]) {
            fillDepth(child, depth + branch[child], children, branch, result);
        }
    }

    private double heightsBelow(int node, int[][] children, double[] branch, double[] heights, double delta) {
        if (children[node].length == 0) {
            return tree.getNodeHeight(tree.getNode(node));
        }
        double height = Double.NEGATIVE_INFINITY;
        for (int child : children[node]) {
            height = Math.max(height, heightsBelow(child, children, branch, heights, delta)
                    + Math.max(branch[child] / rate, delta));
        }
        heights[node] = height;
        return height;
    }

    // each branch's rate b_c / (h(p) - h(c)), within [r / 100, 100 r], set on the node that now holds its clade
    private void setRates(int[][] children, double[] heights, double[] branch, int rootNumber) {
        final BranchRateWarmStart setter = new BranchRateWarmStart(tree, rates, BranchRateWarmStart.Mode.KEEP);
        final java.util.BitSet[] clades = TopologyConstraints.clades(tree);
        final java.util.Map<java.util.BitSet, Integer> holder = new java.util.HashMap<java.util.BitSet, Integer>();
        for (int n = 0; n < clades.length; ++n) {
            holder.put(clades[n], n);
        }
        final java.util.BitSet[] built = new java.util.BitSet[children.length];
        cladesBelow(rootNumber, children, built);
        for (int s = 0; s < children.length; ++s) {
            if (s == rootNumber) {
                continue;
            }
            final double parentHeight = heights[parentOf(s, children)];
            final double height = children[s].length == 0 ? tree.getNodeHeight(tree.getNode(s)) : heights[s];
            final double r = Math.min(RATE_RANGE * rate, Math.max(rate / RATE_RANGE,
                    branch[s] / (parentHeight - height)));
            setter.setRate(tree.getNode(holder.get(built[s])), r);
        }
    }

    private java.util.BitSet cladesBelow(int node, int[][] children, java.util.BitSet[] clades) {
        final java.util.BitSet bits = new java.util.BitSet();
        if (children[node].length == 0) {
            bits.set(node);
        }
        for (int child : children[node]) {
            bits.or(cladesBelow(child, children, clades));
        }
        clades[node] = bits;
        return bits;
    }

    private static int parentOf(int node, int[][] children) {
        for (int p = 0; p < children.length; ++p) {
            for (int c : children[p]) {
                if (c == node) {
                    return p;
                }
            }
        }
        throw new IllegalStateException();
    }
}
