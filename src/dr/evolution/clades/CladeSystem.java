/*
 * CladeSystem.java
 *
 * Copyright © 2002-2026, the BEAST Development Team.
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
 */

package dr.evolution.clades;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.stats.DiscreteStatistics;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * The clades of a set of trees on the same taxa, with how many of the trees each is in.
 *
 * Trees can be added ({@link #add}) and their heights and attributes collected
 * ({@link #collectCladeHeights}, {@link #traverseTree}) from any number of threads at once.
 * The summaries (credibilities, the top clades and so on) are for when that has finished.
 *
 * Clades are identified by the keys of {@link CladeKeys}, by default hashes of the taxon
 * names, so the same clade has the same key in every clade system.
 *
 * @author Andrew Rambaut
 * @version $
 */
public final class CladeSystem {

    private final boolean keepSubClades;
    private final boolean keepParents;
    private final AtomicInteger treeCount = new AtomicInteger(0);

    /**
     * Constructor starting with an empty clade system
     *
     * @param keepSubClades whether to keep all subtrees in each clade
     */
    public CladeSystem(boolean keepSubClades, boolean keepParents) {
        this.keepSubClades = keepSubClades;
        this.keepParents = keepParents;
    }

    /**
     * Constructor adding a single target tree
     */
    public CladeSystem(Tree targetTree, CladeSystem sourceCladeSystem, boolean storeBitsets) {
        this.keepSubClades = false;
        this.keepParents = false;
        add(targetTree, storeBitsets);
        for (Clade clade : this.getClades()) {
            Clade sourceClade = sourceCladeSystem.getClade(clade.getKey());
            if (sourceClade != null) {
                clade.setCredibility(sourceClade.getCredibility());
                clade.setCount(sourceClade.getCount());
            } else {
                // a user-supplied target tree may contain clades not present in the sampled trees
                clade.setCredibility(0.0);
                clade.setCount(0);
            }
        }
    }

    /**
     * adds all the clades in the tree
     */
    public void add(Tree tree) {
        add(tree, false);
    }

    /**
     * adds all the clades in the tree
     */
    public void add(Tree tree, boolean storeBitsets) {
        // the first tree gives the taxa if they have not been set
        final Tips tips = getTips(tree);

        // Recurse over the tree and add all the clades (or increment their
        // frequency if already present). The root clade is added too (for
        // annotation purposes).
        Clade rootClade = addClades(tree, tree.getRoot(), tips, storeBitsets);

        if (this.rootClade == null) {
            // every thread would set the same clade
            this.rootClade = rootClade;
        }
        assert rootClade == this.rootClade;

        assert rootClade.getSize() == tree.getExternalNodeCount();

        treeCount.incrementAndGet();
    }

    /**
     * Sets the taxa, which gives each its index. If this is not called the taxa are taken from
     * the first tree to be added, which with several threads could be any of them, so call
     * this first if the indices should not depend on that.
     *
     * @throws IllegalStateException if the taxa have already been set
     */
    public void setTaxonList(TaxonList taxonList) {
        synchronized (tipsLock) {
            if (tips != null) {
                throw new IllegalStateException("The taxa of the clade system have already been set");
            }
            tips = new Tips(taxonList);
        }
    }

    /**
     * @return the number of trees added
     */
    public int getTreeCount() {
        return treeCount.get();
    }

    public Clade getRootClade() {
        return rootClade;
    }

    private Tips getTips(TaxonList defaultTaxonList) {
        Tips tips = this.tips;
        if (tips == null) {
            synchronized (tipsLock) {
                if (this.tips == null) {
                    this.tips = new Tips(defaultTaxonList);
                }
                tips = this.tips;
            }
        }
        return tips;
    }

    /**
     * The taxa and their tip clades: made once and never changed.
     */
    private static final class Tips {
        private final Map<Taxon, Integer> taxonNumberMap = new HashMap<>();
        // the tip clades by taxon index
        private final BiClade[] tipClades;
        private final Map<Object, BiClade> tipCladeMap = new HashMap<>();

        private Tips(TaxonList taxonList) {
            tipClades = new BiClade[taxonList.getTaxonCount()];
            for (int i = 0; i < tipClades.length; i++) {
                Taxon taxon = taxonList.getTaxon(i);
                taxonNumberMap.put(taxon, i);
                BiClade clade = new BiClade(i, taxon);
                tipClades[i] = clade;
                BiClade sameKey = tipCladeMap.put(clade.getKey(), clade);
                if (sameKey != null) {
                    // the keys are hashes of the names so there is nothing else to try
                    throw new IllegalArgumentException("Taxa '" + sameKey.getTaxon().getId() + "' and '" +
                            taxon.getId() + "' have the same clade key");
                }
            }
        }

        private BiClade getTipClade(Tree tree, NodeRef node) {
            Integer index = taxonNumberMap.get(tree.getNodeTaxon(node));
            if (index == null) {
                throw new IllegalArgumentException("Taxon '" + tree.getNodeTaxon(node).getId() +
                        "' is not among the taxa of the clade system");
            }
            return tipClades[index];
        }
    }

    /**
     * recursively add all the clades in a tree
     */
    private BiClade addClades(Tree tree, NodeRef node, Tips tips, boolean storeBitsets) {

        BiClade clade;
        if (tree.isExternal(node)) {
            clade = tips.getTipClade(tree, node);

            if (storeBitsets && !bitsetKeyMap.containsKey(clade.getKey())) {
                BitsetKey bitset = new BitsetKey(tree.getExternalNodeCount());
                bitset.set(clade.getIndex());
                bitsetKeyMap.putIfAbsent(clade.getKey(), bitset);
            }
        } else {
            assert tree.getChildCount(node) == 2 : "requires a strictly bifurcating tree";

            BiClade clade1 = addClades(tree, tree.getChild(node, 0), tips, storeBitsets);
            BiClade clade2 = addClades(tree, tree.getChild(node, 1), tips, storeBitsets);

            clade = getOrAddClade(clade1, clade2);

            if (keepParents) {
                clade1.addParent(clade);
                clade2.addParent(clade);
            }

            if (storeBitsets) {
                getOrStoreBitset(clade, clade1, clade2);
            }
        }
        assert clade != null;

        clade.incrementCount();

        return clade;
    }

    /**
     * Store the bitset for the parent of two clades
     *
     * @param clade
     * @param clade1
     * @param clade2
     */
    private BitsetKey getOrStoreBitset(BiClade clade, BiClade clade1, BiClade clade2) {
        return getOrStoreBitset(clade.getKey(), clade1.getKey(), clade2.getKey());
    }

    /**
     * Store the bitset for the parent of two clades
     *
     * @param cladeKey
     * @param key1
     * @param key2
     */
    private BitsetKey getOrStoreBitset(Object cladeKey, Object key1, Object key2) {
        BitsetKey bitset = bitsetKeyMap.get(cladeKey);
        if (bitset != null) {
            return bitset;
        }

        return bitsetKeyMap.computeIfAbsent(cladeKey, k -> {
            BitsetKey bitset1 = bitsetKeyMap.get(key1);
            BitsetKey bitset2 = bitsetKeyMap.get(key2);
            assert bitset1 != null && bitset2 != null;

            BitsetKey combined = new BitsetKey(bitset1);
            combined.or(bitset2);
            return combined;
        });
    }

    private final Map<Object, BitsetKey> bitsetKeyMap = new ConcurrentHashMap<>();

    /**
     * see if a clade exists otherwise create it
     */
    private BiClade getOrAddClade(BiClade child1, BiClade child2) {
        Object key = BiClade.getParentKey(child1.getKey(), child2.getKey());
        BiClade clade = cladeMap.get(key);
        if (clade == null) {
            BiClade newClade;
            if (keepSubClades) {
                newClade = new BiClade(child1, child2);
            } else {
                newClade = new BiClade(key, child1.getSize() + child2.getSize());
            }
            clade = cladeMap.putIfAbsent(key, newClade);
            if (clade == null) {
                // this thread's clade is the one that was added
                requiresMrcaKeyList = null;
                return newClade;
            }
        }
        if (keepSubClades) {
            clade.addSubClades(child1, child2);
        }
        return clade;
    }

    /**
     * @return the clade with this key (a tip clade or not), or null if there is none
     */
    public Clade getClade(Object key) {
        Clade clade = cladeMap.get(key);
        if (clade == null) {
            Tips tips = this.tips;
            if (tips != null) {
                clade = tips.tipCladeMap.get(key);
            }
        }
        return clade;
    }

    /**
     * Calls the action for each node of the tree whose clade is in this clade system. The
     * action is called on this thread, so it must be safe to call from several threads if
     * this method is.
     */
    public void traverseTree(Tree tree, CladeAction action) {
        traverseTree(tree, tree.getRoot(), getTips(tree), action);
    }

    private Object traverseTree(Tree tree, NodeRef node, Tips tips, CladeAction action) {

        Object key;

        if (tree.isExternal(node)) {
            key = tips.getTipClade(tree, node).getKey();
        } else {
            assert tree.getChildCount(node) == 2;

            Object key1 = traverseTree(tree, tree.getChild(node, 0), tips, action);
            Object key2 = traverseTree(tree, tree.getChild(node, 1), tips, action);

            key = BiClade.getParentKey(key1, key2);
        }

        Clade clade = getClade(key);
        if (clade != null) {
            action.actOnClade(clade, tree, node);
        } else {
            assert action.expectAllClades();
        }

        return key;
    }

    // use the more efficient algorithm
    private static final boolean MRCA_COLLECTION_BY_NODE = false;

    /**
     * This method traverses the provided tree and then adds the node heights
     * to the respective clades list (for later summarizing and annotating on
     * the target tree). The simple algorithm simply sees if the clade is
     * in the clade system and adds the height to that. This can result in very
     * few height values for unsupported clades and produce negative branch lengths.
     *
     * An alternative is provided when mrcaCladeHeights = true, this is the CA
     * (clade ancestor) method of:
     * Heled and Bouckaert, (2013) Looking for trees in the forest:
     * summary tree from posterior samples'. BMC Evolutionary Biology 13:221;
     *
     * For this, if the clade represented by the node is not present then the
     * node representing the most recent common ancestor of the tips in the clade
     * is found in the tree and that node's height is added to the clade heights.
     * This means that all clades in the clade system get a height for every tree
     * even if it is just the root height.
     *
     * MRCA_COLLECTION_BY_NODE = true
     * This method's algorithm keeps an ordered set of clades that haven't been
     * matched yet. As each node is visited it check to see if this matches any
     * of the clades or is the MRCA of a clade that hasn't been match - in which
     * case the height is added and the clade removed from the `unresolved` set.
     * As this is post-order traversal, the nodes of the tree are larger and
     * larger clades (and ultimately the root) so will eventually be the MRCA of
     * all clades.
     *
     * MRCA_COLLECTION_BY_NODE = false
     * A second, potentially more efficient algorithm `collectCladeHeightsByClade`
     * collects all the clades from the tree (with their heights) once and then
     * iterates through the clades and finds the MRCA height.
     *
     * @param tree
     * @param mrcaCladeHeights
     * @return
     */
    public void collectCladeHeights(Tree tree, boolean mrcaCladeHeights) {
        if (!mrcaCladeHeights) {
            collectCladeHeightsByNode(tree, tree.getRoot(), false, null);
            return;
        }

        // the list is shared by every thread: it is made once and then only read
        List<BiClade> clades = requiresMrcaKeyList;
        if (clades == null) {
            synchronized (mrcaLock) {
                clades = requiresMrcaKeyList;
                if (clades == null) {
                    clades = new ArrayList<>(cladeMap.values());
                    // sort in descending size
                    clades.sort(Comparator.comparingInt((BiClade clade) -> clade.getSize()).reversed());
                    requiresMrcaKeyList = clades;
                }
            }
        }

        if (MRCA_COLLECTION_BY_NODE) {
            // older slower algorithm
            final Set<BiClade> unresolved = new LinkedHashSet<>(clades);
            collectCladeHeightsByNode(tree, tree.getRoot(), true, unresolved);
            assert unresolved.isEmpty();
            return;
        }

        collectCladeHeightsByClade(tree, clades);
    }


    /**
     * The version of the algorithm that iterates through all clades for
     * all nodes in the traversal. Simpler but less efficient.
     *
     * @param tree
     * @param node
     * @param mrcaCladeHeights
     * @param unresolved
     * @return
     */
    private Object collectCladeHeightsByNode(Tree tree,
                                             NodeRef node,
                                             boolean mrcaCladeHeights,
                                             Set<BiClade> unresolved) {

        if (tree.isExternal(node)) {
            return collectTipHeight(tree, node);
        }

        assert tree.getChildCount(node) == 2;

        Object key1 = collectCladeHeightsByNode(tree, tree.getChild(node, 0), mrcaCladeHeights, unresolved);
        Object key2 = collectCladeHeightsByNode(tree, tree.getChild(node, 1), mrcaCladeHeights, unresolved);

        Object key = BiClade.getParentKey(key1, key2);

        BiClade nodeClade = (BiClade) getClade(key);

        if (!mrcaCladeHeights) {
            // only add the node height to the height list of the matching clade
            if (nodeClade != null) {
                // parent is in the Clade set
                nodeClade.addHeightValue(tree.getNodeHeight(node));
            }
            return key;
        }

        final double nodeHeight = tree.getNodeHeight(node);

        // get the bitset for clade and store it in the cache if not there
        BitsetKey nodeBitset = getOrStoreBitset(key, key1, key2);

        if (nodeClade != null) {
            // parent is in the Clade set...
            nodeClade.addHeightValue(nodeHeight);
            unresolved.remove(nodeClade);
        }

        // but this nodeClade may still be the MRCA of other unmatched clades...
        int nodeCladeSize = nodeBitset.cardinality();

        // then check if this node is the MRCA of one of the
        // clades in the unresolved list. If it is then remove
        // it from the list.

        for (Iterator<BiClade> iterator = unresolved.iterator(); iterator.hasNext(); ) {
            BiClade clade = iterator.next();
            if (clade.getSize() < nodeCladeSize) {
                BitsetKey cladeBitset = bitsetKeyMap.get(clade.getKey());
                if (nodeBitset.isSubset(cladeBitset)) {
                    clade.addHeightValue(nodeHeight);
                    iterator.remove();
                }
            }
        }
        return key;
    }

    /**
     * The version of the algorithm that traverses the tree once, and then iterates
     * over the clades once.
     *
     * @param tree
     */
    private void collectCladeHeightsByClade(Tree tree, List<BiClade> clades) {
        final List<NodeSummary> nodeSummaries = new ArrayList<>(tree.getInternalNodeCount());

        collectCladeHeightsByClade(tree, tree.getRoot(), nodeSummaries);

        nodeSummaries.sort(Comparator.comparingInt(summary -> summary.cladeSize));

        for (BiClade clade : clades) {
            final BitsetKey cladeBitset = bitsetKeyMap.get(clade.getKey());
            findMrcaHeight(nodeSummaries, clade, cladeBitset);

        }
    }

    private void findMrcaHeight(List<NodeSummary> nodeSummaries, BiClade clade, BitsetKey cladeBitset) {
        for (NodeSummary nodeSummary : nodeSummaries) {
            if (nodeSummary.cladeSize < clade.getSize()) {
                continue;
            }
            if (nodeSummary.bitset.isSubset(cladeBitset)) {
                clade.addHeightValue(nodeSummary.height);
                return;
            }
        }
        assert false;
    }

    /**
     * Traverse the tree and collect node summaries for each node
     * @param tree
     * @param node
     * @param nodeSummaries
     * @return
     */
    private Object collectCladeHeightsByClade(Tree tree,
                                              NodeRef node,
                                              List<NodeSummary> nodeSummaries) {

        if (tree.isExternal(node)) {
            return collectTipHeight(tree, node);
        }

        assert tree.getChildCount(node) == 2;

        Object key1 = collectCladeHeightsByClade(tree, tree.getChild(node, 0), nodeSummaries);
        Object key2 = collectCladeHeightsByClade(tree, tree.getChild(node, 1), nodeSummaries);

        Object key = BiClade.getParentKey(key1, key2);

        final BitsetKey nodeBitset = getOrStoreBitset(key, key1, key2);
        final int cladeSize = nodeBitset.cardinality();
        final double nodeHeight = tree.getNodeHeight(node);
        nodeSummaries.add(new NodeSummary(nodeBitset, cladeSize, nodeHeight));

        return key;
}

    private Object collectTipHeight(Tree tree, NodeRef node) {
        BiClade tip = getTips(tree).getTipClade(tree, node);
        tip.addHeightValue(tree.getNodeHeight(node));

        return tip.getKey();
    }

private static final class NodeSummary {
    private final BitsetKey bitset;
    private final int cladeSize;
    private final double height;

    private NodeSummary(BitsetKey bitset, int cladeSize, double height) {
        this.bitset = bitset;
        this.cladeSize = cladeSize;
        this.height = height;
    }
}

    // the clades in descending order of size, made when first needed and discarded when a
    // clade is added
    private final Object mrcaLock = new Object();
    private volatile List<BiClade> requiresMrcaKeyList = null;


public void calculateCladeCredibilities(int totalTreesUsed) {
    for (Clade clade : cladeMap.values()) {
        assert clade.getCount() <= totalTreesUsed : "clade.getCount=(" + clade.getCount() +
                ") should be <= totalTreesUsed = (" + totalTreesUsed + ")";

        clade.setCredibility(((double) clade.getCount()) / (double) totalTreesUsed);
    }
}

public double getLogCladeCredibility(Tree tree) {
    final double[] logCladeCredibility = {0.0};
    traverseTree(tree, new CladeAction() {
        @Override
        public void actOnClade(Clade clade, Tree tree, NodeRef node) {
            logCladeCredibility[0] += Math.log(clade.getCredibility());
        }

        @Override
        public boolean expectAllClades() {
            return true;
        }
    });
    return logCladeCredibility[0];
}

public double getMinimumCladeCredibility(Tree tree) {
    final double[] minCladeCredibility = {Double.MAX_VALUE};
    traverseTree(tree, new CladeAction() {
        @Override
        public void actOnClade(Clade clade, Tree tree, NodeRef node) {
            if (clade.getCredibility() < minCladeCredibility[0]) {
                minCladeCredibility[0] = clade.getCredibility();
            }
        }

        @Override
        public boolean expectAllClades() {
            return true;
        }
    });
    return minCladeCredibility[0];
}

public double getMeanCladeCredibility(Tree tree) {
    final double[] minCladeCredibility = {0.0};
    traverseTree(tree, new CladeAction() {
        @Override
        public void actOnClade(Clade clade, Tree tree, NodeRef node) {
            if (clade.getTaxon() == null) {
                minCladeCredibility[0] += clade.getCredibility();
            }
        }

        @Override
        public boolean expectAllClades() {
            return true;
        }
    });
    return minCladeCredibility[0] / tree.getInternalNodeCount();
}

public double getMedianCladeCredibility(Tree tree) {
    final double[] cladeCredibility = new double[tree.getInternalNodeCount()];
    final int[] i = {0};
    traverseTree(tree, new CladeAction() {
        @Override
        public void actOnClade(Clade clade, Tree tree, NodeRef node) {
            if (clade.getTaxon() == null) {
                cladeCredibility[i[0]] = clade.getCredibility();
                i[0] += 1;
            }
        }

        @Override
        public boolean expectAllClades() {
            return true;
        }
    });


    return DiscreteStatistics.median(cladeCredibility);
}

/**
 * Returns the number of clades in the tree with threshold credibility or higher
 *
 * @param tree
 * @param threshold
 * @return
 */
public int getTopCladeCount(Tree tree, double threshold) {
    final int[] count = {0};
    traverseTree(tree, new CladeAction() {
        @Override
        public void actOnClade(Clade clade, Tree tree, NodeRef node) {
            if (clade.getTaxon() == null && clade.getCredibility() > threshold) {
                count[0] += 1;
            }
        }

        @Override
        public boolean expectAllClades() {
            return true;
        }
    });
    return count[0];
}

/**
 * Returns the set of clades in the tree with threshold credibility or higher
 *
 * @param tree
 * @param threshold
 * @return
 */
public Set<BiClade> getTopClades(Tree tree, double threshold) {
    Set<BiClade> clades = new HashSet<>();
    traverseTree(tree, new CladeAction() {
        @Override
        public void actOnClade(Clade clade, Tree tree, NodeRef node) {
            if (clade.getTaxon() == null && clade.getCredibility() > threshold) {
                clades.add((BiClade) clade);
            }
        }

        @Override
        public boolean expectAllClades() {
            return true;
        }
    });
    return clades;
}

/**
 * Returns the number of clades in the clade system with threshold credibility or higher
 *
 * @param threshold
 * @return
 */
public int getTopCladeCount(double threshold) {
    int count = 0;
    for (Clade clade : cladeMap.values()) {
        if (clade.getCredibility() > threshold) {
            count += 1;
        }
    }
    return count;
}

public int getCladeFrequencyCount(int cladeCount) {
    int count = 0;
    for (Clade clade : cladeMap.values()) {
        if (clade.getCount() == cladeCount) {
            count += 1;
        }
    }
    return count;
}

/**
 * Returns the set of clades in the clade system with threshold credibility or higher
 *
 * @param threshold
 * @return
 */
public Set<BiClade> getTopClades(double threshold) {
    Set<BiClade> clades = new HashSet<>();
    for (BiClade clade : cladeMap.values()) {
        if (clade.getSize() == 1 || clade.getCredibility() >= threshold) {
            clades.add(clade);
        }
    }
    return clades;
}

public List<BiClade> getTopCladeList(double threshold) {
    List<BiClade> clades = new ArrayList<>();
    for (BiClade clade : cladeMap.values()) {
        if (clade.getSize() == 1 || clade.getCredibility() >= threshold) {
            clades.add(clade);
        }
    }
    return clades;
}

public int getCladeCount() {
    return cladeMap.keySet().size();
}

public int getCommonCladeCount(CladeSystem referenceCladeSystem) {
    int count = 0;
    for (Object key : cladeMap.keySet()) {
        if (referenceCladeSystem.cladeMap.containsKey(key)) {
            count++;
        }
    }
    return count;
}

    /**
     * @return the tip clades, in the order of the taxa (none if no taxa have been set)
     */
    public Collection<BiClade> getTipClades() {
        Tips tips = this.tips;
        if (tips == null) {
            return Collections.emptyList();
        }
        return Collections.unmodifiableList(Arrays.asList(tips.tipClades));
    }

    /**
     * @return the clades (excluding the tip clades)
     */
    public Collection<BiClade> getClades() {
        return cladeMap.values();
    }

    public Map<Object, BiClade> getCladeMap() {
        return cladeMap;
    }

    //
    // Private stuff
    //

    // the taxa and their tip clades: set once, under the lock
    private final Object tipsLock = new Object();
    private volatile Tips tips = null;

    // a map of clade key to clade (excluding tip clades)
    private final Map<Object, BiClade> cladeMap = new ConcurrentHashMap<>();

    // the same clade whichever tree it is taken from
    private volatile Clade rootClade;

}