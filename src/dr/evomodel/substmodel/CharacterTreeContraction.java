/*
 * CharacterTreeContraction.java
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

package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.ModelListener;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

/**
 * A multifurcating character tree held as a binary skeleton (a TreeModel) and a collapse state. A collapsed internal
 * node has clade rate exactly 0, which leaves the generator unchanged (docs/generator-constraints), so it is
 * transparent: its children belong to its parent in the contracted tree. The collapse state is either fixed (one flag
 * per node, set where a multifurcating Newick was resolved) or a keep Parameter by internal node index (node number
 * minus the state count; 0 = collapsed). The root is never collapsed.
 * <p>
 * The BEAGLE upload is one canonical binary resolution of the contracted tree, independent of how the skeleton
 * resolves a multifurcation and of child order: a node with three or more contracted children is resolved lightest
 * first (fewest states, ties by smallest state), the inserted nodes have rate 0, the two children of every node are
 * ordered by smallest state, and internal nodes are numbered S + post-order rank.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class CharacterTreeContraction {

    private final TreeModel tree;
    private final int stateCount;
    private final int[] stateOfNode;
    private final boolean[] fixedCollapsed;
    private final Parameter keep;

    private final int[] minState;
    private final double[] triples;
    private final int[] lambdaSource;
    private final double[] lastKeep;
    private int[] keptInternalNodes;
    private int contractedInternalCount;
    private double logResolutionCount;
    private int rank;
    private boolean topologyKnown = false;

    public CharacterTreeContraction(TreeModel tree, DataType dataType, boolean[] fixedCollapsed, Parameter keep) {
        this.tree = tree;
        this.stateCount = tree.getExternalNodeCount();
        this.stateOfNode = new int[stateCount];
        for (int i = 0; i < stateCount; i++) {
            stateOfNode[i] = AbstractHierarchicalRateProvider.stateOfCode(dataType, tree.getTaxonId(i));
        }
        if (fixedCollapsed != null && fixedCollapsed.length != tree.getNodeCount()) {
            throw new IllegalArgumentException("collapse flags need one entry per node of the character tree ("
                    + tree.getNodeCount() + "); found " + fixedCollapsed.length);
        }
        if (keep != null && keep.getDimension() != tree.getInternalNodeCount()) {
            throw new IllegalArgumentException("keep needs one entry per internal node of the character tree ("
                    + tree.getInternalNodeCount() + "); found " + keep.getDimension());
        }
        if (fixedCollapsed != null && keep != null) {
            throw new IllegalArgumentException("a character tree has fixed collapse flags or a keep parameter, not both");
        }
        if (fixedCollapsed != null && fixedCollapsed[tree.getRoot().getNumber()]) {
            throw new IllegalArgumentException("the root of the character tree cannot be collapsed");
        }
        this.fixedCollapsed = fixedCollapsed == null ? null : fixedCollapsed.clone();
        this.keep = keep;
        this.minState = new int[tree.getNodeCount()];
        this.triples = new double[Math.max(1, 3 * (stateCount - 1))];
        this.lambdaSource = new int[Math.max(1, stateCount - 1)];
        this.lastKeep = keep == null ? null : new double[keep.getDimension()];
        if (lastKeep != null) {
            Arrays.fill(lastKeep, Double.NaN);
        }
        final ModelListener listener = new ModelListener() {
            public void modelChangedEvent(Model model, Object object, int index) {
                if (HierarchicalSubstitutionModel.changesTopology(object)) {
                    topologyKnown = false;
                }
            }

            public void modelRestored(Model model) {
                topologyKnown = false;
            }
        };
        tree.addModelListener(listener);
        tree.addModelRestoreListener(listener);
    }

    public TreeModel getTree() {
        return tree;
    }

    /**
     * Marks the canonical resolution stale after a change of the skeleton's topology. The contraction also listens to
     * the skeleton; this makes the order of the listeners irrelevant.
     */
    public void topologyChanged() {
        topologyKnown = false;
    }

    public Parameter getKeep() {
        return keep;
    }

    public boolean hasFixedFlags() {
        return fixedCollapsed != null;
    }

    public int getStateCount() {
        return stateCount;
    }

    public boolean isCollapsed(int nodeNumber) {
        if (nodeNumber < stateCount || nodeNumber == tree.getRoot().getNumber()) {
            return false;
        }
        if (fixedCollapsed != null) {
            return fixedCollapsed[nodeNumber];
        }
        return keep != null && keep.getParameterValue(nodeNumber - stateCount) == 0.0;
    }

    private void ensure() {
        boolean known = topologyKnown;
        if (keep != null) {
            for (int i = 0; i < lastKeep.length; i++) {
                final double k = keep.getParameterValue(i);
                if (!(k == lastKeep[i])) {
                    lastKeep[i] = k;
                    known = false;
                }
            }
        }
        if (!known) {
            compute();
        }
    }

    private void compute() {
        fillMinState(tree.getRoot());
        rank = 0;
        contractedInternalCount = 0;
        logResolutionCount = 0.0;
        emit(resolve(tree.getRoot().getNumber()));
        if (rank != stateCount - 1) {
            throw new IllegalStateException("the canonical resolution has " + rank + " internal nodes, expected "
                    + (stateCount - 1));
        }
        final int[] kept = new int[contractedInternalCount];
        int k = 0;
        for (int r = 0; r < stateCount - 1; r++) {
            if (lambdaSource[r] >= 0) {
                kept[k++] = lambdaSource[r];
            }
        }
        keptInternalNodes = kept;
        topologyKnown = true;
    }

    /**
     * The canonical post-order triples (parent, left, right) of the upload, 3 (S - 1) values.
     */
    public double[] getTriples() {
        ensure();
        return triples;
    }

    /**
     * For each internal node of the upload, by post-order rank: the skeleton node whose rate it carries, or -1 for
     * a node inserted to resolve a multifurcation (rate 0).
     */
    public int[] getLambdaSource() {
        ensure();
        return lambdaSource;
    }

    /**
     * The kept internal nodes of the skeleton (the internal nodes of the contracted tree) in the canonical order of
     * the upload.
     */
    public int[] getKeptInternalNodes() {
        ensure();
        return keptInternalNodes;
    }

    public int getContractedInternalNodeCount() {
        ensure();
        return contractedInternalCount;
    }

    /**
     * The log of the number of binary trees that resolve the contracted tree: the sum over its internal nodes v of
     * log (2 d_v - 3)!!, with d_v the number of contracted children.
     */
    public double getLogResolutionCount() {
        ensure();
        return logResolutionCount;
    }

    public static double logDoubleFactorialOfOddChildren(int d) {
        double sum = 0.0;
        for (int k = 3; k <= 2 * d - 3; k += 2) {
            sum += Math.log(k);
        }
        return sum;
    }

    private int fillMinState(NodeRef node) {
        final int n = node.getNumber();
        if (tree.isExternal(node)) {
            minState[n] = stateOfNode[n];
        } else {
            int m = Integer.MAX_VALUE;
            for (int c = 0; c < tree.getChildCount(node); c++) {
                m = Math.min(m, fillMinState(tree.getChild(node, c)));
            }
            minState[n] = m;
        }
        return minState[n];
    }

    /**
     * The children of a node in the contracted tree (skeleton node numbers), found by descending through collapsed
     * skeleton children, in skeleton order.
     */
    public List<Integer> contractedChildren(int nodeNumber) {
        final List<Integer> kids = new ArrayList<Integer>();
        collect(tree.getNode(nodeNumber), kids);
        return kids;
    }

    /**
     * The contracted children ordered by their smallest state.
     */
    public List<Integer> sortedContractedChildren(int nodeNumber) {
        ensure();
        final List<Integer> kids = contractedChildren(nodeNumber);
        Collections.sort(kids, new Comparator<Integer>() {
            public int compare(Integer a, Integer b) {
                return Integer.compare(minState[a], minState[b]);
            }
        });
        return kids;
    }

    /**
     * The nearest kept ancestor of a node (its parent in the contracted tree), or -1 at the root.
     */
    public int contractedParent(int nodeNumber) {
        NodeRef node = tree.getParent(tree.getNode(nodeNumber));
        while (node != null && isCollapsed(node.getNumber())) {
            node = tree.getParent(node);
        }
        return node == null ? -1 : node.getNumber();
    }

    public int getMinState(int nodeNumber) {
        ensure();
        return minState[nodeNumber];
    }

    private void collect(NodeRef node, List<Integer> kids) {
        for (int c = 0; c < tree.getChildCount(node); c++) {
            final NodeRef child = tree.getChild(node, c);
            if (!tree.isExternal(child) && isCollapsed(child.getNumber())) {
                collect(child, kids);
            } else {
                kids.add(child.getNumber());
            }
        }
    }

    /**
     * The states below a skeleton node in ascending order.
     */
    public int[] sortedStatesBelow(int nodeNumber) {
        final List<Integer> states = new ArrayList<Integer>();
        collectStates(tree.getNode(nodeNumber), states);
        final int[] out = new int[states.size()];
        for (int i = 0; i < out.length; i++) {
            out[i] = states.get(i);
        }
        Arrays.sort(out);
        return out;
    }

    private void collectStates(NodeRef node, List<Integer> states) {
        if (tree.isExternal(node)) {
            states.add(stateOfNode[node.getNumber()]);
        } else {
            for (int c = 0; c < tree.getChildCount(node); c++) {
                collectStates(tree.getChild(node, c), states);
            }
        }
    }

    private static final class Resolved {
        final Resolved left, right;
        final int state, min, size, source;

        Resolved(int state) {
            this.left = null;
            this.right = null;
            this.state = state;
            this.min = state;
            this.size = 1;
            this.source = -1;
        }

        Resolved(Resolved a, Resolved b, int source) {
            final boolean swap = b.min < a.min;
            this.left = swap ? b : a;
            this.right = swap ? a : b;
            this.state = -1;
            this.min = Math.min(a.min, b.min);
            this.size = a.size + b.size;
            this.source = source;
        }
    }

    private static final Comparator<Resolved> LIGHTEST_FIRST = new Comparator<Resolved>() {
        public int compare(Resolved a, Resolved b) {
            return a.size != b.size ? Integer.compare(a.size, b.size) : Integer.compare(a.min, b.min);
        }
    };

    private Resolved resolve(int nodeNumber) {
        if (nodeNumber < stateCount) {
            return new Resolved(stateOfNode[nodeNumber]);
        }
        contractedInternalCount++;
        final List<Integer> kids = contractedChildren(nodeNumber);
        logResolutionCount += logDoubleFactorialOfOddChildren(kids.size());
        final List<Resolved> pieces = new ArrayList<Resolved>();
        for (int kid : kids) {
            pieces.add(resolve(kid));
        }
        while (pieces.size() > 2) {
            Collections.sort(pieces, LIGHTEST_FIRST);
            final Resolved a = pieces.remove(0);
            final Resolved b = pieces.remove(0);
            pieces.add(new Resolved(a, b, -1));
        }
        return new Resolved(pieces.get(0), pieces.get(1), nodeNumber);
    }

    private int emit(Resolved r) {
        if (r.state >= 0) {
            return r.state;
        }
        final int left = emit(r.left);
        final int right = emit(r.right);
        final int id = stateCount + rank;
        triples[3 * rank] = id;
        triples[3 * rank + 1] = left;
        triples[3 * rank + 2] = right;
        lambdaSource[rank] = r.source;
        rank++;
        return id;
    }

    public interface NodeLabeler {
        String label(int nodeNumber);
    }

    /**
     * The contracted tree in Newick, children ordered by smallest state, so equal trees give equal strings.
     */
    public String contractedNewick(NodeLabeler labeler) {
        ensure();
        final StringBuilder sb = new StringBuilder();
        writeContracted(tree.getRoot().getNumber(), labeler, sb);
        return sb.append(';').toString();
    }

    private void writeContracted(int nodeNumber, NodeLabeler labeler, StringBuilder sb) {
        if (nodeNumber < stateCount) {
            sb.append(tree.getTaxonId(nodeNumber));
            return;
        }
        final List<Integer> kids = sortedContractedChildren(nodeNumber);
        sb.append('(');
        for (int k = 0; k < kids.size(); k++) {
            if (k > 0) {
                sb.append(',');
            }
            writeContracted(kids.get(k), labeler, sb);
        }
        sb.append(')');
        if (labeler != null) {
            sb.append(labeler.label(nodeNumber));
        }
    }

    /**
     * A skeleton built from a possibly multifurcating Newick: its fixed collapse flags (true at the nodes inserted to
     * resolve a multifurcation) and, for each internal node of the skeleton, the post-order index of its clade in
     * the Newick as written (-1 for an inserted node).
     */
    public static final class Built {
        public final DefaultTreeModel skeleton;
        public final boolean[] collapsed;
        public final int[] slotOfInternal;
        public final int contractedInternalCount;

        Built(DefaultTreeModel skeleton, boolean[] collapsed, int[] slotOfInternal, int contractedInternalCount) {
            this.skeleton = skeleton;
            this.collapsed = collapsed;
            this.slotOfInternal = slotOfInternal;
            this.contractedInternalCount = contractedInternalCount;
        }

        public int collapsedCount() {
            int n = 0;
            for (boolean c : collapsed) {
                n += c ? 1 : 0;
            }
            return n;
        }

        public boolean isMultifurcating() {
            return collapsedCount() > 0;
        }
    }

    private static final class Piece {
        final String newick;
        final BitSet clade;
        final int min;

        Piece(String newick, BitSet clade, int min) {
            this.newick = newick;
            this.clade = clade;
            this.min = min;
        }
    }

    /**
     * Builds the binary skeleton of a possibly multifurcating Newick. Binary nodes keep the Newick's order, so a
     * binary Newick gives exactly the DefaultTreeModel that AbstractHierarchicalRateProvider.characterTree builds. A
     * node with three or more children is resolved lightest first (random == null) or by a uniformly random pairing
     * (random != null, for tests); the inserted nodes are collapsed. Branch lengths are set to 1.
     */
    public static Built build(String id, String newick, DataType dataType, boolean fixTree, Random random) {
        final Tree imported;
        try {
            imported = new NewickImporter(newick).importTree(null);
        } catch (Exception e) {
            throw new IllegalArgumentException("cannot read the character tree " + newick + ": " + e.getMessage());
        }
        final List<BitSet> inserted = new ArrayList<BitSet>();
        final List<BitSet> postOrder = new ArrayList<BitSet>();
        final Piece root = resolve(imported, imported.getRoot(), dataType, inserted, random, postOrder);
        final DefaultTreeModel skeleton;
        if (inserted.isEmpty()) {
            skeleton = new DefaultTreeModel(id, imported, true, fixTree);
        } else {
            try {
                skeleton = new DefaultTreeModel(id, new NewickImporter(root.newick + ";").importTree(null), true, fixTree);
            } catch (Exception e) {
                throw new IllegalStateException("cannot read the resolved character tree " + root.newick + ": " + e.getMessage());
            }
            AbstractHierarchicalRateProvider.setUnitBranchLengths(skeleton);
        }
        final Set<BitSet> insertedSet = new HashSet<BitSet>(inserted);
        final boolean[] collapsed = new boolean[skeleton.getNodeCount()];
        final BitSet[] clades = new BitSet[skeleton.getNodeCount()];
        cladesOf(skeleton, skeleton.getRoot(), dataType, clades);
        final int tips = skeleton.getExternalNodeCount();
        final int[] slotOfInternal = new int[skeleton.getInternalNodeCount()];
        for (int i = tips; i < skeleton.getNodeCount(); i++) {
            collapsed[i] = insertedSet.contains(clades[i]);
            slotOfInternal[i - tips] = collapsed[i] ? -1 : postOrder.indexOf(clades[i]);
            if (!collapsed[i] && slotOfInternal[i - tips] < 0) {
                throw new IllegalStateException("clade " + clades[i] + " of the skeleton is not in the Newick");
            }
        }
        return new Built(skeleton, collapsed, slotOfInternal, postOrder.size());
    }

    /**
     * The star tree on the states of a data type, tips in state order.
     */
    public static String starNewick(DataType dataType) {
        final StringBuilder sb = new StringBuilder("(");
        for (int i = 0; i < dataType.getStateCount(); i++) {
            final String code = AbstractHierarchicalRateProvider.codeOfState(dataType, i);
            sb.append(i > 0 ? "," : "").append(code.matches("[A-Za-z0-9_.\\-]+") ? code : "'" + code + "'");
        }
        return sb.append(");").toString();
    }

    private static BitSet cladesOf(Tree tree, NodeRef node, DataType dataType, BitSet[] clades) {
        final BitSet set = new BitSet();
        if (tree.isExternal(node)) {
            set.set(AbstractHierarchicalRateProvider.stateOfCode(dataType, tree.getTaxonId(node.getNumber())));
        } else {
            for (int c = 0; c < tree.getChildCount(node); c++) {
                set.or(cladesOf(tree, tree.getChild(node, c), dataType, clades));
            }
        }
        clades[node.getNumber()] = set;
        return set;
    }

    public static boolean isMultifurcating(String newick) {
        final Tree imported;
        try {
            imported = new NewickImporter(newick).importTree(null);
        } catch (Exception e) {
            throw new IllegalArgumentException("cannot read the character tree " + newick + ": " + e.getMessage());
        }
        for (int i = 0; i < imported.getInternalNodeCount(); i++) {
            if (imported.getChildCount(imported.getInternalNode(i)) > 2) {
                return true;
            }
        }
        return false;
    }

    private static Piece resolve(Tree tree, NodeRef node, DataType dataType, List<BitSet> inserted, Random random,
                                 List<BitSet> postOrder) {
        if (tree.isExternal(node)) {
            final String id = tree.getTaxonId(node.getNumber());
            final int state = AbstractHierarchicalRateProvider.stateOfCode(dataType, id);
            if (state < 0 || state >= dataType.getStateCount()) {
                throw new IllegalArgumentException("character tree tip " + id + " is not a state of "
                        + dataType.getDescription());
            }
            final BitSet set = new BitSet();
            set.set(state);
            return new Piece(id.matches("[A-Za-z0-9_.\\-]+") ? id : "'" + id + "'", set, state);
        }
        final List<Piece> kids = new ArrayList<Piece>();
        for (int c = 0; c < tree.getChildCount(node); c++) {
            kids.add(resolve(tree, tree.getChild(node, c), dataType, inserted, random, postOrder));
        }
        if (kids.size() == 1) {
            return kids.get(0);
        }
        final BitSet clade = new BitSet();
        for (Piece kid : kids) {
            clade.or(kid.clade);
        }
        postOrder.add(clade);
        if (kids.size() == 2) {
            return join(kids.get(0), kids.get(1));
        }
        final List<Piece> pool = new ArrayList<Piece>(kids);
        if (random == null) {
            while (pool.size() > 2) {
                Collections.sort(pool, new Comparator<Piece>() {
                    public int compare(Piece a, Piece b) {
                        final int sa = a.clade.cardinality(), sb = b.clade.cardinality();
                        return sa != sb ? Integer.compare(sa, sb) : Integer.compare(a.min, b.min);
                    }
                });
                final Piece a = pool.remove(0);
                final Piece b = pool.remove(0);
                final Piece joined = a.min < b.min ? join(a, b) : join(b, a);
                inserted.add(joined.clade);
                pool.add(joined);
            }
            return pool.get(0).min < pool.get(1).min ? join(pool.get(0), pool.get(1)) : join(pool.get(1), pool.get(0));
        }
        while (pool.size() > 2) {
            final Piece a = pool.remove(random.nextInt(pool.size()));
            final Piece b = pool.remove(random.nextInt(pool.size()));
            final Piece joined = random.nextBoolean() ? join(a, b) : join(b, a);
            inserted.add(joined.clade);
            pool.add(random.nextInt(pool.size() + 1), joined);
        }
        return join(pool.get(0), pool.get(1));
    }

    private static Piece join(Piece a, Piece b) {
        final BitSet set = (BitSet) a.clade.clone();
        set.or(b.clade);
        return new Piece("(" + a.newick + ":1," + b.newick + ":1)", set, Math.min(a.min, b.min));
    }
}
