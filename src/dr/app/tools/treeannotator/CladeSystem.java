/*
 * CladeSystem.java
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

package dr.app.tools.treeannotator;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.stats.DiscreteStatistics;

import java.util.*;
import java.util.concurrent.*;

/**
 * @author Andrew Rambaut
 * @version $
 */
final class CladeSystem {
    private final boolean keepSubClades;
    private double treeCount = 0;

    /**
     * Constructor starting with an empty clade system
     * @param keepSubClades whether to keep all subtrees in each clade
     */
    public CladeSystem(boolean keepSubClades) {
        this.keepSubClades = keepSubClades;
    }

    /**
     * Constructor adding a single target tree
     */
    public CladeSystem(Tree targetTree) {
        this.keepSubClades = false;
        add(targetTree);
    }
    /**
     * adds all the clades in the tree
     */
    public void add(Tree tree) {
        if (taxonList == null) {
            setTaxonList(tree);
        }

        if (treeCount == 0) {
            // these will always be the same so create them once
            addTipClades(tree);
        }

        // Recurse over the tree and add all the clades (or increment their
        // frequency if already present). The root clade is added too (for
        // annotation purposes).
        rootClade = addClades(tree, tree.getRoot());

        assert rootClade.getSize() == tree.getExternalNodeCount();

        treeCount += 1;
    }

    public void setTaxonList(TaxonList taxonList) {
        this.taxonList = taxonList;
        taxonNumberMap = new HashMap<>();
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
            Clade clade = new BiClade(index, taxon);
            tipClades.put(index, clade);
        }
    }

    /**
     * recursively add all the clades in a tree
     */
    private Clade addClades(Tree tree, NodeRef node) {
        Clade clade;
        if (tree.isExternal(node)) {
            // all tip clades should already be there
            int index = node.getNumber();
            if (taxonNumberMap != null) {
                index = taxonNumberMap.get(tree.getNodeTaxon(node));
            }
            clade = tipClades.get(index);
//            assert clade != null && clade.getTaxon().equals(tree.getNodeTaxon(node));
        } else {
            assert tree.getChildCount(node) == 2 : "requires a strictly bifurcating tree";

            Clade clade1 = addClades(tree, tree.getChild(node, 0));
            Clade clade2 = addClades(tree, tree.getChild(node, 1));
            clade = getOrAddClade(clade1, clade2);
        }
        assert clade != null;

        clade.setCount(clade.getCount() + 1);

        return clade;
    }

    /**
     * see if a clade exists otherwise create it
     */
    private Clade getOrAddClade(Clade child1, Clade child2) {
        Object key = BiClade.makeKey(child1.getKey(), child2.getKey());
        BiClade clade = (BiClade)cladeMap.get(key);
        if (clade == null) {
            if (keepSubClades) {
                clade = new BiClade(child1, child2);
            } else {
                clade = new BiClade(key, child1.getSize() + child2.getSize());
            }
            cladeMap.put(clade.getKey(), clade);
        } else {
            if (keepSubClades) {
                clade.addSubClades(child1, child2);
            }
        }
        return clade;
    }

    public Clade getClade(Object key) {
        if (key instanceof Integer) {
            return tipClades.get(key);
        }
        return cladeMap.get(key);
    }

    public void traverseTree(Tree tree, CladeAction action) {
        traverseTree(tree, tree.getRoot(), action);
    }

    private Object traverseTree(Tree tree, NodeRef node, CladeAction action) {

        Object key;

        if (tree.isExternal(node)) {
            key = node.getNumber();
            if (taxonNumberMap != null) {
                key = taxonNumberMap.get(tree.getNodeTaxon(node));
            }
        } else {
            assert tree.getChildCount(node) == 2;

            Object key1 = traverseTree(tree, tree.getChild(node, 0), action);
            Object key2 = traverseTree(tree, tree.getChild(node, 1), action);

            key = BiClade.makeKey(key1, key2);
        }

        Clade clade = getClade(key);
        if (clade != null) {
            action.actOnClade(clade, tree, node);
        } else {
            assert action.expectAllClades();
        }

        return key;
    }

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
        final int[] i = { 0 };
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
     * @param tree
     * @param threshold
     * @return
     */
    public int getTopCladeCount(Tree tree, double threshold) {
        final int[] count = {0};
        traverseTree(tree, new CladeAction() {
            @Override
            public void actOnClade(Clade clade, Tree tree, NodeRef node) {
                if (clade.getTaxon() == null && clade.getCredibility() >= threshold) {
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
     * @param tree
     * @param threshold
     * @return
     */
    public Set<Clade> getTopClades(Tree tree, double threshold) {
        Set<Clade> clades = new HashSet<>();
        traverseTree(tree, new CladeAction() {
            @Override
            public void actOnClade(Clade clade, Tree tree, NodeRef node) {
                if (clade.getTaxon() == null && clade.getCredibility() >= threshold) {
                    clades.add(clade);
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
     * @param threshold
     * @return
     */
    public int getTopCladeCount(double threshold) {
        int count = 0;
        for (Clade clade : cladeMap.values()) {
            if (clade.getCredibility() >= threshold) {
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
     * @param threshold
     * @return
     */
    public Set<Clade> getTopClades(double threshold) {
        Set<Clade> clades = new HashSet<>();
        for (Clade clade : cladeMap.values()) {
            if (clade.getCredibility() >= threshold) {
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
                count ++;
            }
        }
        return count;
    }

    final int binCount = 100;
    final Map<Object, Clade>[] cladeMapBySize = new Map[binCount];
    int binWidth;

    private void binCladesBySize() {
        int n = taxonList.getTaxonCount();
        binWidth = (n / binCount) + 1;
        for (Map.Entry<Object, Clade> entry : cladeMap.entrySet()) {
            int bin = entry.getValue().getSize() / binWidth;
            if (cladeMapBySize[bin] == null) {
                cladeMapBySize[bin] = new HashMap<>();
            }
            cladeMapBySize[bin].put(entry.getKey(), entry.getValue());
        }
    }

    private Clade getCladeBySize(Object key, int size) {
        int bin = size / binWidth;
        return cladeMapBySize[bin].get(key);
    }

    public void embiggenBiClades(final int minCladeSize, final int minCladeCount) {
        binCladesBySize();

        List<Clade> allClades = new ArrayList<>(cladeMap.values());
        allClades.addAll(tipClades.values());
        Clade[] clades = new Clade[allClades.size()];
        clades = allClades.toArray(clades);

        // sort by number of trees containing clade
        Arrays.sort(clades, (o1, o2) -> o2.getCount() - o1.getCount());
        int n = 0;
        // find the point at which the count drops below minCount
        while (n < clades.length && clades[n].getCount() > minCladeCount - 1) {
            n++;
        }

        // truncate the array at this pont
        clades = Arrays.copyOf(clades, n);

        // sort by descending size
        Arrays.sort(clades, (o1, o2) -> o2.getSize() - o1.getSize());

        int maxSize = clades[0].getSize();

        int[] sizeIndices = new int[maxSize];
        n = 0;
        int currentSize = maxSize;
        while (n < clades.length) {
            if (clades[n].getSize() < currentSize) {
                currentSize -= 1;
                sizeIndices[currentSize] = n;
            }
            n++;
        }
        sizeIndices[0] = clades.length;

        n = sizeIndices[minCladeSize - 1];

        long x = (((((long)n) - 1) * n) / 2);

        System.err.println("Expanding with " + x + " clade pairs");
        System.err.println("0              25             50             75            100");
        System.err.println("|--------------|--------------|--------------|--------------|");

        long stepSize = Math.max(x / 60, 1);

        long k = 0;

        BitSet bits = new BitSet();

        for (int i = sizeIndices[maxSize - 1]; i < n - 1; i++) {
            BiClade clade1 = (BiClade)clades[i];
            if (clade1.getSize() >= 2) {
                BitSet bits1 = ((BitSet) clade1.getKey());
                for (int j = Math.max(i + 1, sizeIndices[maxSize - clade1.getSize()]); j < n; j++) {
                    BiClade clade2 = (BiClade) clades[j];
                    BiClade clade = null;

                    bits.clear();
                    bits.or(bits1);

                    Object key2 = clade2.getKey();
                    if (key2 instanceof Integer) {
                        bits.set((Integer) key2);
                    } else {
                        bits.or((BitSet) key2);
                    }

                    int size = clade1.getSize() + clade2.getSize();
                    if (bits.cardinality() == size) {
//                        clade = (BiClade) cladeMap.get(bits);
                        clade = (BiClade) getCladeBySize(bits, size);
                        if (clade != null) {
                            clade.addSubClades(clade1, clade2);
                        }
                    }

                    if (k > 0 && k % stepSize == 0) {
                        System.err.print("*");
                        System.err.flush();
                    }
                    k++;
                }
            }
        }
        System.err.println();
        System.err.println("" + k + " additional clade pairs examined");
    }

    public void embiggenBiClades(final int minCladeSize, final int minCladeCount, final int threadCount) {
        List<Clade> allClades = new ArrayList<>(cladeMap.values());
        allClades.addAll(tipClades.values());
        Clade[] cladeArray = new Clade[allClades.size()];
        cladeArray = allClades.toArray(cladeArray);

        // sort by number of trees containing clade
        Arrays.sort(cladeArray, (o1, o2) -> o2.getCount() - o1.getCount());
        int n = 0;
        // find the point at which the count drops below minCount
        while (n < cladeArray.length && cladeArray[n].getCount() > minCladeCount - 1) {
            n++;
        }

        // truncate the array at this pont
        cladeArray = Arrays.copyOf(cladeArray, n);

        // sort by descending size
        Arrays.sort(cladeArray, (o1, o2) -> o2.getSize() - o1.getSize());

        int maxSize = cladeArray[0].getSize();

        // store the indices where the clades drop in size
        // this will be used to skip to where clades are small enough to be subclades
        int[] sizeIndices = new int[maxSize];
        n = 0;
        int currentSize = maxSize;
        while (n < cladeArray.length) {
            if (cladeArray[n].getSize() < currentSize) {
                currentSize -= 1;
                sizeIndices[currentSize] = n;
            }
            n++;
        }
        sizeIndices[0] = cladeArray.length;

        n = sizeIndices[minCladeSize - 1];

        // make the clade array final so it can be accessed in an anonymous class
        final Clade[] clades = cladeArray;

        // count the exact number of clade pairs (for reporting purposes)
        long count = 0;
        for (int u = sizeIndices[maxSize - 1]; u < n - 1; u++) {
            BiClade clade1 = (BiClade) clades[u];
            count += n - Math.max(u + 1, sizeIndices[maxSize - clade1.getSize()]);
        }

        System.err.printf("Embiggening with up to %,d clade pairs...", count);
        System.err.println();
        System.err.println("0              25             50             75            100");
        System.err.println("|--------------|--------------|--------------|--------------|");
        final int stepSize = Math.max((n - 1 - sizeIndices[maxSize - 1]) / 60, 1);

        ExecutorService pool;
        if (threadCount <= 0) {
            pool = Executors.newCachedThreadPool();
        } else {
            pool = Executors.newFixedThreadPool(threadCount);
        }
        List<Future<?>> futures = new ArrayList<>();

        final int[] k = { 0 };
        for (int i = sizeIndices[maxSize - 1]; i < n - 1; i++) {
            BiClade clade1 = (BiClade)clades[i];
            // clade1 must be more than just a tip...
            if (clade1.getSize() >= 2) {
                // get the bitset for clade1 that acts as its hash key
                final BitSet bits1 = ((BitSet) clade1.getKey());
                // get final versions of the start and end of the iteration
                final int from = Math.max(i + 1, sizeIndices[maxSize - clade1.getSize()]);
                final int to = n;
                // submit the thread to the pool and store the future in the list
                futures.add(pool.submit(() -> {
                    // create and reuse a bitset to avoid reallocating it
                    BitSet bits = new BitSet();
                    for (int j = from; j < to; j++) {
                        BiClade clade2 = (BiClade) clades[j];

                        BiClade clade = null;

                        // clear the bitset and make a copy of clade1's bits
                        bits.clear();
                        bits.or(bits1);

                        // get clade2's bits and add them to the bitset
                        Object key2 = clade2.getKey();
                        if (key2 instanceof Integer) {
                            bits.set((Integer) key2);
                        } else {
                            bits.or((BitSet) key2);
                        }

                        // if the cardinality of the bitset is not the same as the sum
                        // of the sizes of the two clades then they must have had some
                        // tips in common.
                        if (bits.cardinality() == clade1.getSize() + clade2.getSize()) {
                            clade = (BiClade) cladeMap.get(bits);
                            if (clade != null) {
                                clade.addSubClades(clade1, clade2);
                            }
                        }
                    }
                    synchronized (k) {
                        if (k[0] > 0 && k[0] % stepSize == 0) {
                            System.err.print("*");
                            System.err.flush();
                        }
                        k[0]++;
                    }

                }));
            }
        }

        try {
            // wait for all the threads to run to completion
            for (Future<?> f: futures) {
                f.get();
            }
        } catch (InterruptedException | ExecutionException e) {
            throw new RuntimeException(e);
        }

        System.err.println("  ");

    }


    //
    // Private stuff
    //
    TaxonList taxonList = null;
    Map<Taxon, Integer> taxonNumberMap = null;

    private final Map<Object, Clade> tipClades = new HashMap<>();
    private final Map<Object, Clade> cladeMap = new HashMap<>();

    Clade rootClade;

}
