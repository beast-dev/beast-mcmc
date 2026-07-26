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

package dr.app.tools.treeannotator;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.stats.DiscreteStatistics;

import java.util.*;

/**
 * @author Andrew Rambaut
 * @version $
 */
public final class CladeSystem {
    private enum MrcaCollectionMode {
        LEGACY,
        ALTERNATIVE_BY_CLADE
    }

    // Compile-time strategy switch for MRCA clade height collection.
    private static final MrcaCollectionMode MRCA_COLLECTION_MODE = MrcaCollectionMode.ALTERNATIVE_BY_CLADE;
    private static final boolean PARALLELIZE_ALTERNATIVE_BY_CLADE = true;
    private static final int PARALLELIZE_ALTERNATIVE_MIN_CLADES = 64;
    private static final boolean ENABLE_MRCA_TIMING = true;
    private static final int MRCA_TIMING_REPORT_EVERY = 1000;

    private static long legacyMrcaCallCount = 0;
    private static long alternativeMrcaCallCount = 0;
    private static long legacyMrcaNanos = 0;
    private static long alternativeMrcaNanos = 0;

    private final boolean keepSubClades;
    private final boolean keepParents;
    private double treeCount = 0;

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
            assert sourceClade != null;
            clade.setCredibility(sourceClade.getCredibility());
            clade.setCount(sourceClade.getCount());
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
        synchronized (taxonNumberMap) {
            if (taxonList == null) {
                setTaxonList(tree);
            }
        }

        if (treeCount == 0) {
            // these will always be the same so create them once
            synchronized (tipClades) {
                addTipClades(tree);
            }
        }

        // Recurse over the tree and add all the clades (or increment their
        // frequency if already present). The root clade is added too (for
        // annotation purposes).
        Clade rootClade = addClades(tree, tree.getRoot(), storeBitsets);
        if (this.rootClade == null) {
            this.rootClade = rootClade;
        }
        assert rootClade == this.rootClade;
        assert rootClade.getSize() == tree.getExternalNodeCount();

        treeCount += 1;
    }

    public void setTaxonList(TaxonList taxonList) {
        this.taxonList = taxonList;
        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            taxonNumberMap.put(taxonList.getTaxon(i), i);
        }
    }

    public Clade getRootClade() {
        return rootClade;
    }

    /**
     * add all the tips in a tree
     */
    private void addTipClades(Tree tree) {
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef tip = tree.getExternalNode(i);
            int index = tip.getNumber();
            Taxon taxon = tree.getNodeTaxon(tip);
            if (taxonNumberMap != null) {
                index = taxonNumberMap.get(taxon);
            }
            BiClade clade = new BiClade(index, taxon);
            tipClades.put(index, clade);
            tipCladeMap.put(clade.getKey(), clade);
        }
    }

    /**
     * recursively add all the clades in a tree
     */
    private BiClade addClades(Tree tree, NodeRef node, boolean storeBitsets) {
        BiClade clade;
        if (tree.isExternal(node)) {
            // all tip clades should already be there
            int index = node.getNumber();
            if (taxonNumberMap != null) {
                index = taxonNumberMap.get(tree.getNodeTaxon(node));
            }
            clade = tipClades.get(index);
            assert clade != null && clade.getTaxon().equals(tree.getNodeTaxon(node));

            if (storeBitsets) {
                BitsetKey bitset = new BitsetKey(tree.getExternalNodeCount());
                bitset.set(index);
                bitsetKeyMap.put(clade.getKey(), bitset);
            }
        } else {
            assert tree.getChildCount(node) == 2 : "requires a strictly bifurcating tree";

            BiClade clade1 = addClades(tree, tree.getChild(node, 0), storeBitsets);
            BiClade clade2 = addClades(tree, tree.getChild(node, 1), storeBitsets);
            synchronized (cladeMap) {
                clade = getOrAddClade(clade1, clade2);
            }

            if (keepParents) {
                clade1.addParent(clade);
                clade2.addParent(clade);
            }

            if (storeBitsets) {
                getOrStoreBitset(clade, clade1, clade2);
            }
        }
        assert clade != null;

        synchronized (clade) {
            clade.setCount(clade.getCount() + 1);
        }

        return clade;
    }

    /**
     * Store the bitset for the parent of two clades
     *
     * @param cladeKey
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
        if (bitset == null) {
            BitsetKey bitset1 = bitsetKeyMap.get(key1);
            BitsetKey bitset2 = bitsetKeyMap.get(key2);
            assert bitset1 != null && bitset2 != null;

            bitset = new BitsetKey(bitset1);
            bitset.or(bitset2);

            bitsetKeyMap.put(cladeKey, bitset);
        }
        return bitset;
    }

    private final Map<Object, BitsetKey> bitsetKeyMap = new HashMap<>();

    /**
     * see if a clade exists otherwise create it
     */
    private BiClade getOrAddClade(Clade child1, Clade child2) {
        Object key = BiClade.getParentKey(child1.getKey(), child2.getKey());
        BiClade clade = cladeMap.get(key);
        if (clade == null) {
            if (keepSubClades) {
                clade = new BiClade(child1, child2);
            } else {
                clade = new BiClade(key, child1.getSize() + child2.getSize());
            }
            cladeMap.put(clade.getKey(), clade);
            requiresMrcaKeyList = null;
        } else {
            synchronized (clade) {
                if (keepSubClades) {
                    clade.addSubClades(child1, child2);
                }
            }
//            }
        }

        return clade;
    }

    public Clade getClade(Object key) {
        Clade clade = cladeMap.get(key);
        if (clade == null) {
            clade = tipCladeMap.get(key);
        }
//        assert clade != null;
        return clade;
    }

    public void traverseTree(Tree tree, CladeAction action) {
        traverseTree(tree, tree.getRoot(), action);
    }

    private Object traverseTree(Tree tree, NodeRef node, CladeAction action) {

        Object key;

        if (tree.isExternal(node)) {
            int index = node.getNumber();
            if (taxonNumberMap != null) {
                index = taxonNumberMap.get(tree.getNodeTaxon(node));
            }
            key = BiClade.getTaxonKey(index);
        } else {
            assert tree.getChildCount(node) == 2;

            Object key1 = traverseTree(tree, tree.getChild(node, 0), action);
            Object key2 = traverseTree(tree, tree.getChild(node, 1), action);

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

//    public void traverseNonBinaryTree(Tree tree, CladeAction action) {
//        traverseNonBinaryTree(tree, tree.getRoot(), action);
//    }

//    private Object traverseNonBinaryTree(Tree tree, NodeRef node, CladeAction action) {
//
//        Object key;
//
//        if (tree.isExternal(node)) {

    /// /            key = node.getNumber();
//            key = taxonNumberMap.get(tree.getNodeTaxon(node));
//        } else {
//            List<Object> keys = new ArrayList<>();
//            for (int i = 0; i < tree.getChildCount(node); i++) {
//                keys.add(traverseNonBinaryTree(tree, tree.getChild(node, i), action));
//            }
//            key = BiClade.getParentKey(keys.toArray());
//        }
//
//        Clade clade = getClade(key);
//        if (clade != null) {
//            action.actOnClade(clade, tree, node);
//        } else {
//            assert action.expectAllClades();
//        }
//
//        return key;
//    }
    public void collectCladeHeights(Tree tree, boolean mrcaCladeHeights) {
        if (!mrcaCladeHeights) {
            collectCladeHeightsLegacy(tree, tree.getRoot(), false);
            return;
        }

        ensureRequiresMrcaKeyList();
        final long mrcaStartTime = ENABLE_MRCA_TIMING ? System.nanoTime() : 0;

        if (MRCA_COLLECTION_MODE == MrcaCollectionMode.LEGACY) {
            requiresMrcaKeySet.clear();
            requiresMrcaKeySet.addAll(requiresMrcaKeyList);
            collectCladeHeightsLegacy(tree, tree.getRoot(), true);

            assert requiresMrcaKeySet.isEmpty();
            if (ENABLE_MRCA_TIMING) {
                recordMrcaTiming(MrcaCollectionMode.LEGACY, System.nanoTime() - mrcaStartTime);
            }
            return;
        }

        collectCladeHeightsAlternative(tree);
        if (ENABLE_MRCA_TIMING) {
            recordMrcaTiming(MrcaCollectionMode.ALTERNATIVE_BY_CLADE, System.nanoTime() - mrcaStartTime);
        }
    }

    private static synchronized void recordMrcaTiming(MrcaCollectionMode mode, long elapsedNanos) {
        if (mode == MrcaCollectionMode.LEGACY) {
            legacyMrcaCallCount++;
            legacyMrcaNanos += elapsedNanos;
        } else {
            alternativeMrcaCallCount++;
            alternativeMrcaNanos += elapsedNanos;
        }

        final long totalCalls = legacyMrcaCallCount + alternativeMrcaCallCount;
        if (MRCA_TIMING_REPORT_EVERY > 0 && totalCalls % MRCA_TIMING_REPORT_EVERY == 0) {
            final double legacyMeanMs = legacyMrcaCallCount == 0 ? 0.0 : (legacyMrcaNanos / 1_000_000.0) / legacyMrcaCallCount;
            final double alternativeMeanMs = alternativeMrcaCallCount == 0 ? 0.0 : (alternativeMrcaNanos / 1_000_000.0) / alternativeMrcaCallCount;

            System.err.println("MRCA timing summary [calls=" + totalCalls + "]: " +
                    "legacy_calls=" + legacyMrcaCallCount +
                    ", legacy_mean_ms=" + legacyMeanMs +
                    ", alt_calls=" + alternativeMrcaCallCount +
                    ", alt_mean_ms=" + alternativeMeanMs);
        }
    }

    public static synchronized String getMrcaTimingSummary() {
        if (!ENABLE_MRCA_TIMING) {
            return "MRCA timing is disabled.";
        }

        final long totalCalls = legacyMrcaCallCount + alternativeMrcaCallCount;
        final double legacyMeanMs = legacyMrcaCallCount == 0 ? 0.0 : (legacyMrcaNanos / 1_000_000.0) / legacyMrcaCallCount;
        final double alternativeMeanMs = alternativeMrcaCallCount == 0 ? 0.0 : (alternativeMrcaNanos / 1_000_000.0) / alternativeMrcaCallCount;

        return "MRCA timing final summary [calls=" + totalCalls + "]: " +
                "legacy_calls=" + legacyMrcaCallCount +
                ", legacy_mean_ms=" + legacyMeanMs +
                ", alt_calls=" + alternativeMrcaCallCount +
                ", alt_mean_ms=" + alternativeMeanMs;
    }

    private void ensureRequiresMrcaKeyList() {
        if (requiresMrcaKeyList == null) {
            requiresMrcaKeyList = new ArrayList<>(cladeMap.values());
            // sort in descending size
            requiresMrcaKeyList.sort((o1, o2) -> o2.getSize() - o1.getSize());
        }
    }

    private Object collectCladeHeightsLegacy(Tree tree, NodeRef node, boolean mrcaCladeHeights) {

        Object key;

        if (tree.isExternal(node)) {
            int index = node.getNumber();
            if (taxonNumberMap != null) {
                index = taxonNumberMap.get(tree.getNodeTaxon(node));
            }
            key = BiClade.getTaxonKey(index);

            BiClade tip = (BiClade) getClade(key);
            tip.addHeightValue(tree.getNodeHeight(node));

            if (mrcaCladeHeights) {
                ensureTipBitset(tree, key, index);
            } else {
                BitsetKey bitset = bitsetKeyMap.get(key);
                assert bitset != null;
            }
        } else {
            assert tree.getChildCount(node) == 2;

            Object key1 = collectCladeHeightsLegacy(tree, tree.getChild(node, 0), mrcaCladeHeights);
            Object key2 = collectCladeHeightsLegacy(tree, tree.getChild(node, 1), mrcaCladeHeights);

            key = BiClade.getParentKey(key1, key2);

            BiClade nodeClade = (BiClade) getClade(key);

            if (mrcaCladeHeights) {
                final double nodeHeight = tree.getNodeHeight(node);
                // if the clade matching the node doesn't exist then find the
                // node which is the MRCA of the tips in the clade and add that
                // as the height of the clade.
                // This is the CA (clade ancestor) method of:
                // Heled and Bouckaert, (2013) Looking for trees in the forest:
                // summary tree from posterior samples'. BMC Evolutionary Biology 13:221;

                // get the bitset for clade and store it in the cache if not there
                BitsetKey nodeBitset = getOrStoreBitset(key, key1, key2);

                if (nodeClade != null) {
                    // parent is in the Clade set...
                    nodeClade.addHeightValue(nodeHeight);
                    requiresMrcaKeySet.remove(nodeClade);
                }

                // but this nodeClade may still be the MRCA of other unmatched clades...
                int nodeCladeSize = nodeBitset.cardinality();

                // then check if this node is the MRCA of one of the
                // clades in the list.

                for (Iterator<BiClade> iterator = requiresMrcaKeySet.iterator(); iterator.hasNext(); ) {
                    BiClade clade = iterator.next();
                    if (clade.getSize() < nodeCladeSize) {
                        BitsetKey cladeBitset = bitsetKeyMap.get(clade.getKey());
                        if (nodeBitset.isSubset(cladeBitset)) {
                            clade.addHeightValue(nodeHeight);
                            iterator.remove();
                        }
                    }
                }

//                }
            } else {
                // only add the node height to the height list of the matching clade
                if (nodeClade != null) {
                    // parent is in the Clade set
                    nodeClade.addHeightValue(tree.getNodeHeight(node));
                }
            }

        }
        return key;
    }

    private void collectCladeHeightsAlternative(Tree tree) {
        final Set<BiClade> unresolved = new HashSet<>(requiresMrcaKeyList);
        final List<NodeSummary> nodeSummaries = new ArrayList<>(tree.getInternalNodeCount());

        collectCladeHeightsAlternative(tree, tree.getRoot(), unresolved, nodeSummaries);

        nodeSummaries.sort(Comparator.comparingInt(summary -> summary.cladeSize));

        final List<BiClade> unresolvedList = new ArrayList<>(unresolved);
        final List<BiClade> stillUnresolved = Collections.synchronizedList(new ArrayList<>());

        java.util.stream.Stream<BiClade> unresolvedStream = unresolvedList.stream();
        if (PARALLELIZE_ALTERNATIVE_BY_CLADE && unresolvedList.size() >= PARALLELIZE_ALTERNATIVE_MIN_CLADES) {
            unresolvedStream = unresolvedStream.parallel();
        }

        unresolvedStream.forEach(clade -> {
            final BitsetKey cladeBitset = bitsetKeyMap.get(clade.getKey());
            if (cladeBitset == null) {
                stillUnresolved.add(clade);
                return;
            }

            final Double mrcaHeight = findMrcaHeight(nodeSummaries, clade, cladeBitset);
            if (mrcaHeight == null) {
                stillUnresolved.add(clade);
                return;
            }

            clade.addHeightValue(mrcaHeight);
        });

        assert stillUnresolved.isEmpty() : "Failed to resolve MRCA for " + stillUnresolved.size() + " clades";
    }

    private Object collectCladeHeightsAlternative(Tree tree,
                                                  NodeRef node,
                                                  Set<BiClade> unresolved,
                                                  List<NodeSummary> nodeSummaries) {

        Object key;

        if (tree.isExternal(node)) {
            int index = node.getNumber();
            if (taxonNumberMap != null) {
                index = taxonNumberMap.get(tree.getNodeTaxon(node));
            }
            key = BiClade.getTaxonKey(index);

            BiClade tip = (BiClade) getClade(key);
            tip.addHeightValue(tree.getNodeHeight(node));

            ensureTipBitset(tree, key, index);
        } else {
            assert tree.getChildCount(node) == 2;

            Object key1 = collectCladeHeightsAlternative(tree, tree.getChild(node, 0), unresolved, nodeSummaries);
            Object key2 = collectCladeHeightsAlternative(tree, tree.getChild(node, 1), unresolved, nodeSummaries);

            key = BiClade.getParentKey(key1, key2);

            final BitsetKey nodeBitset = getOrStoreBitset(key, key1, key2);
            final int cladeSize = nodeBitset.cardinality();
            final double nodeHeight = tree.getNodeHeight(node);
            nodeSummaries.add(new NodeSummary(nodeBitset, cladeSize, nodeHeight));

            final BiClade nodeClade = (BiClade) getClade(key);
            if (nodeClade != null) {
                nodeClade.addHeightValue(nodeHeight);
                unresolved.remove(nodeClade);
            }
        }

        return key;
    }

    private void ensureTipBitset(Tree tree, Object key, int index) {
        if (bitsetKeyMap.containsKey(key)) {
            return;
        }

        final BitsetKey bitset = new BitsetKey(tree.getExternalNodeCount());
        bitset.set(index);
        bitsetKeyMap.put(key, bitset);
    }

    private Double findMrcaHeight(List<NodeSummary> nodeSummaries, BiClade clade, BitsetKey cladeBitset) {
        for (NodeSummary nodeSummary : nodeSummaries) {
            if (nodeSummary.cladeSize <= clade.getSize()) {
                continue;
            }
            if (nodeSummary.bitset.isSubset(cladeBitset)) {
                return nodeSummary.height;
            }
        }
        return null;
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

    private List<BiClade> requiresMrcaKeyList = null;
    private final Set<BiClade> requiresMrcaKeySet = new LinkedHashSet<>();


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

    Collection<BiClade> getTipClades() {
        return tipClades.values();
    }

    Collection<BiClade> getClades() {
        return cladeMap.values();
    }

    Map<Object, BiClade> getCladeMap() {
        return cladeMap;
    }

    //
    // Private stuff
    //
    TaxonList taxonList = null;
    private final Map<Taxon, Integer> taxonNumberMap = new HashMap<>();

    // a map of taxon index to clade
    private final Map<Integer, BiClade> tipClades = new HashMap<>();
    private final Map<Object, BiClade> tipCladeMap = new HashMap<>();

    // a map of clade key to clade (excluding tip clades)
    private final Map<Object, BiClade> cladeMap = new HashMap<>();

    Clade rootClade;

}
