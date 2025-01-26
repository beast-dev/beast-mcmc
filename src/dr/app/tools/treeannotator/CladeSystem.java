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
     *
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
        Clade rootClade = addClades(tree, tree.getRoot());
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
            synchronized (cladeMap) {
                clade = getOrAddClade(clade1, clade2);
            }
        }
        assert clade != null;

        synchronized (clade) {
            clade.setCount(clade.getCount() + 1);
        }

        return clade;
    }

    /**
     * see if a clade exists otherwise create it
     */
    private Clade getOrAddClade(Clade child1, Clade child2) {
        Object key = BiClade.makeKey(child1.getKey(), child2.getKey());
        BiClade clade = (BiClade) cladeMap.get(key);
        if (clade == null) {
            if (keepSubClades) {
                clade = new BiClade(child1, child2);
            } else {
                clade = new BiClade(key, child1.getSize() + child2.getSize());
            }
            cladeMap.put(clade.getKey(), clade);
        } else {
            synchronized (clade) {
                if (keepSubClades) {
                    clade.addSubClades(child1, child2);
                }
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
     *
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
     *
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
     *
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
                count++;
            }
        }
        return count;
    }

    Map<Object, BiClade>[] cladeMapBySize = null;
    //    Map<Integer, Set<BiClade>>[] cladeSetByTipBySize = null;
    int binWidth;

    private void binCladesBySize() {
        int maxSize = taxonList.getTaxonCount();
//        int binCount = 805;
//        binWidth = maxSize / binCount;
        binWidth = 2;
        int binCount = maxSize / binWidth;
        cladeMapBySize = new Map[binCount + 2];

//        cladeSetByTipBySize = new Map[maxSize + 1];
//        for (int i = 0; i <= maxSize; i++) {
//            cladeSetByTipBySize[i] = new HashMap<>();
//        }

        for (Map.Entry<Object, Clade> entry : cladeMap.entrySet()) {
            BiClade clade = (BiClade) entry.getValue();
            if (clade.size > 1) {
                int bin = clade.size / binWidth;
                if (cladeMapBySize[bin] == null) {
                    cladeMapBySize[bin] = new HashMap<>();
                }
                cladeMapBySize[bin].put(entry.getKey(), clade);

//                BitSet bits = (BitSet) clade.key;
//                for (int i = bits.nextSetBit(0); i >= 0; i = bits.nextSetBit(i + 1)) {
////                for (int j = clade.size; j <= maxSize; j++ ) {
//                    Set<BiClade> cladeSet = cladeSetByTipBySize[clade.size].getOrDefault(i, null);
//                    if (cladeSet == null) {
//                        cladeSet = new HashSet<>();
//                        cladeSetByTipBySize[clade.size].put(i, cladeSet);
//                    }
//                    cladeSet.add(clade);
////                }
//                    // if (i == Integer.MAX_VALUE) {          break; // or (i+1) would overflow      }
//                }
            }
        }
    }

//    private boolean isSubClade(BiClade clade) {
//        BitSet bits = (BitSet)clade.key;
//        Set<BiClade> candidates = null;
//        for (int i = bits.nextSetBit(0); i >= 0; i = bits.nextSetBit(i+1)) {
//            Set<BiClade> cladeSet = new HashSet<>();
//            for (int j = clade.size + 1; j < cladeSetByTipBySize.length; j++ ) {
//                Set<BiClade> c = cladeSetByTipBySize[j].get(i);
//                if (c != null) {
//                    cladeSet.addAll(c);
//                }
//            }
//            if (candidates == null) {
//                candidates = new HashSet<>(cladeSet);
//            } else {
//                candidates.retainAll(cladeSet);
//                if (candidates.isEmpty()) {
//                    return false;
//                }
//            }
//        }
//        return true;
//    }
//
//    /**
//     * get all the larger clades with that contain all the taxa in clade
//     * @param clade
//     * @return
//     */
//    private Set<BiClade> getSuperClades(BiClade clade) {
//        BitSet bits = (BitSet)clade.key;
//        Set<BiClade> candidates = null;
//        for (int i = bits.nextSetBit(0); i >= 0; i = bits.nextSetBit(i+1)) {
//            Set<BiClade> cladeSet = new HashSet<>();
//            for (int j = clade.size + 1; j < cladeSetByTipBySize.length; j++ ) {
//                Set<BiClade> c = cladeSetByTipBySize[j].get(i);
//                if (c != null) {
//                    cladeSet.addAll(c);
//                }
//            }
//            if (candidates == null) {
//                candidates = new HashSet<>(cladeSet);
//            } else {
//                candidates.retainAll(cladeSet);
//                if (candidates.isEmpty()) {
//                    return null;
//                }
//            }
//        }
//        return candidates;
//    }

    private BiClade getCladeBySize(Object key, int size) {
        int bin = size / binWidth;
        if (cladeMapBySize[bin] == null) {
            return null;
        }
        return cladeMapBySize[bin].get(key);
    }

    /**
     * Version of embiggening that is single threaded and used
     * for testing optimisation.   
     * @param minCladeSize
     * @param minCladeCount
     */
    public void embiggenBiClades(final int minCladeSize, final int minCladeCount) {
        binCladesBySize();

        List<Clade> allClades = new ArrayList<>(cladeMap.values());
        allClades.addAll(tipClades.values());
        BiClade[] clades = new BiClade[allClades.size()];
        clades = allClades.toArray(clades);

        // sort by number of trees containing clade
        Arrays.sort(clades, (o1, o2) -> o2.getCount() - o1.getCount());
        int n = 0;
        // find the point at which the count drops below minCount
        while (n < clades.length && clades[n].getCount() > minCladeCount - 1) {
            n++;
        }

        // truncate the array at this point
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


//        n = sizeIndices[Math.max(1, minCladeSize - 1)];
//        int y = 0;
//        boolean[] notSubset = new boolean[clades.length];
//        for (int i = n - 1; i > 0; i--) {
//            Clade clade1 = clades[i];
//            notSubset[i] = true;
//            for (int j = i - 1; j > 0; j--) {
//                Clade clade2 = clades[j];
//                if (isSubset((BitSet) clade1.getKey(), (BitSet) clade2.getKey())) {
//                    notSubset[i] = false;
//                    y ++;
//                    break;
//                }
//            }
//        }
//        System.err.println(n + " clades tested, " + y + " are subsets of larger clades");

        n = sizeIndices[minCladeSize - 1];

//        System.err.println();
//        System.err.println("Finding superclades...");
//        System.err.println("0              25             50             75            100");
//        System.err.println("|--------------|--------------|--------------|--------------|");
//
//        long stepSize = Math.max(n / 60, 1);
//        long k = 0;
//        Set<BiClade>[] superClades = new Set[n];
//        for (int i = sizeIndices[maxSize - 1]; i < n - 1; i++) {
//            BiClade clade = clades[i];
//            if (clade.getSize() >= 2) {
//                superClades[i] = getSuperClades(clade);
//            }
//            if (k > 0 && k % stepSize == 0) {
//                System.err.print("*");
//                System.err.flush();
//            }
//            k++;
//        }

        // count the exact number of clade pairs (for reporting purposes)
        long count = 0;
        for (int u = sizeIndices[maxSize - 1]; u < n - 1; u++) {
            BiClade clade1 = clades[u];
            count += n - Math.max(u + 1, sizeIndices[maxSize - clade1.getSize()]);
        }

        System.err.printf("Embiggening with up to %,d clade pairs...", count);
        System.err.println();
        System.err.println("0              25             50             75            100");
        System.err.println("|--------------|--------------|--------------|--------------|");

        long stepSize = Math.max(count / 60, 1);
        long k = 0;
        long embiggulationCount = 0;

        // create and reuse a bitset to avoid reallocating it
        final BitSet bits = new BitSet();

        int rejectCount = 0;

        for (int i = sizeIndices[maxSize - 1]; i < n - 1; i++) {
//            if (notSubset[i]) {
//                rejectCount += 1;
//                continue;
//            }
            BiClade clade1 = clades[i];
            if (clade1.getSize() >= 2) {
//                Set<BiClade> superClades1 = superClades[i];
//                if (superClades1 == null) {
//                    continue;
//                }
                BitSet bits1 = ((BitSet) clade1.getKey());

                for (int j = Math.max(i + 1, sizeIndices[maxSize - clade1.getSize()]); j < n; j++) {
                    BiClade clade2 = clades[j];
                    if (clade2.index == clade1.index) {
                        continue;
                    }

//                    if (clade2.size >= 2) {
//                        Set<BiClade> superClades2 = superClades[j];
//                        superClades2.retainAll(superClades1);
//                        if (superClades2.isEmpty()) {
//                            continue;
//                        }
//                    }

                    bits.clear();
                    bits.or(bits1);

                    if (clade2.key instanceof Integer) {
                        bits.set((Integer) clade2.key);
                    } else {
                        bits.or((BitSet) clade2.key);
                    }

                    int size = clade1.size + clade2.size;
                    if (bits.cardinality() == size) {
//                        BiClade clade = (BiClade) cladeMap.get(bits);
                        BiClade clade = getCladeBySize(bits, size);
                        if (clade != null) {
                            clade.addSubClades(clade1, clade2);
                            embiggulationCount++;
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
        System.err.println(k + " additional clade pairs examined");
        System.err.println(embiggulationCount+" additional clade pairs added");
//        System.err.println(rejectCount + " rejected");
    }

    private boolean isSubset(BitSet bits1, BitSet bits2) {
        for (int i = bits1.nextSetBit(0); i >= 0; i = bits1.nextSetBit(i+1)) {
            if (!bits2.get(i)) {
                return false;
            }
            // if (i == Integer.MAX_VALUE) {          break; // or (i+1) would overflow      }
        }
        return true;
    }

    public void embiggenBiClades(final int minCladeSize, final int minCladeCount, final int threadCount) {
        binCladesBySize();

        // pull down the clades and tip clades into an array to iterate over easily
        List<Clade> allClades = new ArrayList<>(cladeMap.values());
        allClades.addAll(tipClades.values());
        BiClade[] cladeArray = new BiClade[allClades.size()];
        cladeArray = allClades.toArray(cladeArray);

        // sort by number of trees containing clade
        Arrays.sort(cladeArray, (o1, o2) -> o2.getCount() - o1.getCount());
        int n = 0;
        // find the point at which the count drops below minCount
        while (n < cladeArray.length && cladeArray[n].getCount() > minCladeCount - 1) {
            n++;
        }

        // truncate the array at this point
        cladeArray = Arrays.copyOf(cladeArray, n);

        // sort by descending size
        Arrays.sort(cladeArray, (o1, o2) -> o2.getSize() - o1.getSize());

        // get the size of the biggest clade - which will be the root clade
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
        final BiClade[] clades = cladeArray;

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
        final int[] embiggulationCount = { 0 };
        for (int i = sizeIndices[maxSize - 1]; i < n - 1; i++) {
            BiClade clade1 = clades[i];
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
                    final BitSet bits = new BitSet();
                    for (int j = from; j < to; j++) {
                        BiClade clade2 = clades[j];

                        BiClade clade = null;

                        // clear the bitset and make a copy of clade1's bits
                        bits.clear();
                        bits.or(bits1);

                        // get clade2's bits and add them to the bitset
                        if (clade2.key instanceof Integer) {
                            bits.set((Integer) clade2.key);
                        } else {
                            bits.or((BitSet) clade2.key);
                        }

                        // if the cardinality of the bitset is not the same as the sum
                        // of the sizes of the two clades then they must have had some
                        // tips in common.
                        int size = clade1.size + clade2.size;
                        if (bits.cardinality() == size) {
//                            clade = (BiClade)cladeMap.get(bits);
                            clade = getCladeBySize(bits, size);
                            if (clade != null) {
                                clade.addSubClades(clade1, clade2);
                                synchronized (embiggulationCount) {
                                    embiggulationCount[0]++;
                                }
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
        System.err.println(embiggulationCount[0] + " additional clade pairs added");

        System.err.print("  ");

    }


    //
    // Private stuff
    //
    TaxonList taxonList = null;
    private final Map<Taxon, Integer> taxonNumberMap = new HashMap<>();

    private final Map<Object, Clade> tipClades = new HashMap<>();
    private final Map<Object, Clade> cladeMap = new HashMap<>();

    Clade rootClade;

}
