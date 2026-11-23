/*
 * CladeSystemTest.java
 *
 * Copyright © 2026, the BEAST Development Team.
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

package test.dr.evolution.clades;

import dr.evolution.clades.BiClade;
import dr.evolution.clades.Clade;
import dr.evolution.clades.CladeAction;
import dr.evolution.clades.CladeSystem;
import dr.evolution.clades.TaxonHashCladeKeys;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.util.Pair;
import junit.framework.TestCase;

import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * @author Andrew Rambaut
 */
public class CladeSystemTest extends TestCase {

    private static Tree tree(String newick) throws Exception {
        return new NewickImporter(newick).importTree(null);
    }

    private static long key(String... taxa) {
        long key = 0;
        for (String taxon : taxa) {
            key ^= TaxonHashCladeKeys.getTaxonKey(taxon);
        }
        return key;
    }

    /**
     * Random rooted bifurcating trees with node heights, as Newick strings. The tips are
     * written in a random order so trees differ in how their taxa are numbered.
     */
    private static List<Tree> randomTrees(int treeCount, int taxonCount, long seed) throws Exception {
        Random random = new Random(seed);
        List<Tree> trees = new ArrayList<Tree>();
        for (int t = 0; t < treeCount; t++) {
            List<String> subtrees = new ArrayList<String>();
            List<Double> heights = new ArrayList<Double>();
            for (int i = 0; i < taxonCount; i++) {
                subtrees.add("T" + i);
                heights.add(0.0);
            }
            double height = 0.0;
            while (subtrees.size() > 1) {
                // mostly join neighbours so that the trees share many clades
                int i = random.nextInt(subtrees.size() - 1);
                int j = (random.nextInt(4) == 0 ? random.nextInt(subtrees.size()) : i + 1);
                if (i == j) {
                    j = i + 1;
                }
                height += 0.1 + random.nextDouble();
                String left = subtrees.get(i) + ":" + (height - heights.get(i));
                String right = subtrees.get(j) + ":" + (height - heights.get(j));
                String joined = random.nextBoolean() ? "(" + left + "," + right + ")" : "(" + right + "," + left + ")";
                subtrees.set(Math.min(i, j), joined);
                heights.set(Math.min(i, j), height);
                subtrees.remove(Math.max(i, j));
                heights.remove(Math.max(i, j));
            }
            trees.add(tree(subtrees.get(0) + ";"));
        }
        return trees;
    }

    /**
     * Everything about a clade system that should not depend on how it was filled.
     */
    private static Map<Object, String> describe(CladeSystem cladeSystem) {
        Map<Object, String> description = new TreeMap<Object, String>();
        List<BiClade> clades = new ArrayList<BiClade>(cladeSystem.getClades());
        clades.addAll(cladeSystem.getTipClades());
        for (BiClade clade : clades) {
            List<String> subClades = new ArrayList<String>();
            if (clade.getSubClades() != null) {
                for (Pair<BiClade, BiClade> pair : clade.getSubClades()) {
                    long first = (Long) pair.first.getKey();
                    long second = (Long) pair.second.getKey();
                    subClades.add(Math.min(first, second) + "+" + Math.max(first, second));
                }
            }
            Collections.sort(subClades);
            List<Long> parents = new ArrayList<Long>();
            if (clade.getParentClades() != null) {
                for (BiClade parent : clade.getParentClades()) {
                    parents.add((Long) parent.getKey());
                }
            }
            Collections.sort(parents);
            List<Double> heights = new ArrayList<Double>(clade.getHeightValues());
            Collections.sort(heights);
            description.put(clade.getKey(), "size " + clade.getSize() + " count " + clade.getCount() +
                    " subclades " + subClades + " parents " + parents + " heights " + heights);
        }
        return description;
    }

    public void testCountsAndCredibilities() throws Exception {
        CladeSystem cladeSystem = new CladeSystem(false, false);
        cladeSystem.add(tree("((A:1,B:1):1,(C:1,D:1):1);"));
        cladeSystem.add(tree("((B:1,A:1):1,(D:1,C:1):1);"));
        cladeSystem.add(tree("((A:1,C:1):1,(B:1,D:1):1);"));
        cladeSystem.calculateCladeCredibilities(3);

        assertEquals(3, cladeSystem.getTreeCount());
        // AB, CD, AC, BD and the root
        assertEquals(5, cladeSystem.getCladeCount());
        assertEquals(4, cladeSystem.getTipClades().size());

        Clade ab = cladeSystem.getClade(key("A", "B"));
        assertEquals(2, ab.getCount());
        assertEquals(2.0 / 3.0, ab.getCredibility(), 1E-12);
        assertEquals(2, ab.getSize());

        assertEquals(1, cladeSystem.getClade(key("A", "C")).getCount());
        assertEquals(1, cladeSystem.getClade(key("B", "D")).getCount());
        assertNull(cladeSystem.getClade(key("A", "D")));

        Clade root = cladeSystem.getClade(key("A", "B", "C", "D"));
        assertSame(cladeSystem.getRootClade(), root);
        assertEquals(3, root.getCount());
        assertEquals(1.0, root.getCredibility(), 1E-12);

        Clade a = cladeSystem.getClade(TaxonHashCladeKeys.getTaxonKey("A"));
        assertEquals("A", a.getTaxon().getId());
        assertEquals(3, a.getCount());
    }

    /**
     * A clade's key depends only on the names of its taxa: not on the order the taxa come in,
     * how the clade is built up, or which clade system it is in.
     */
    public void testKeysDependOnlyOnNames() throws Exception {
        assertEquals(key("A") ^ key("B"), TaxonHashCladeKeys.getParentKey(key("A"), key("B")));
        assertEquals(BiClade.getParentKey(key("A"), key("B")), BiClade.getParentKey(key("B"), key("A")));
        assertTrue(key("A") != key("B"));
        // similar names give unrelated keys
        assertTrue(Long.bitCount(key("T1") ^ key("T2")) > 8);

        CladeSystem first = new CladeSystem(false, false);
        first.add(tree("(((A:1,B:1):1,C:2):1,D:3);"));
        CladeSystem second = new CladeSystem(false, false);
        second.add(tree("(D:3,(C:2,(B:1,A:1):1):1);"));

        assertEquals(first.getCladeMap().keySet(), second.getCladeMap().keySet());
        assertTrue(first.getCladeMap().containsKey(key("A", "B", "C")));
        assertEquals(describe(first), describe(second));
    }

    /**
     * Adding the trees from many threads gives the same clades, counts, sub-clades and
     * parents as adding them from one.
     */
    public void testAddFromManyThreads() throws Exception {
        final List<Tree> trees = randomTrees(400, 16, 17);

        CladeSystem serial = new CladeSystem(true, true);
        for (Tree tree : trees) {
            serial.add(tree);
        }
        Map<Object, String> expected = describe(serial);
        assertTrue("the trees should share clades and differ", expected.size() > 60);

        ExecutorService pool = Executors.newFixedThreadPool(8);
        try {
            for (int repeat = 0; repeat < 25; repeat++) {
                final CladeSystem threaded = new CladeSystem(true, true);
                List<Future<?>> futures = new ArrayList<Future<?>>();
                for (final Tree tree : trees) {
                    futures.add(pool.submit(() -> threaded.add(tree)));
                }
                for (Future<?> future : futures) {
                    future.get();
                }
                assertEquals(trees.size(), threaded.getTreeCount());
                assertEquals("repeat " + repeat, expected, describe(threaded));
                assertEquals(serial.getRootClade().getKey(), threaded.getRootClade().getKey());
            }
        } finally {
            pool.shutdown();
        }
    }

    /**
     * Collecting heights (by clade, and as the height of each clade's common ancestor) and
     * traversing the trees from many threads gives the same as from one.
     */
    public void testCollectFromManyThreads() throws Exception {
        final List<Tree> trees = randomTrees(200, 12, 5);
        final Tree target = trees.get(0);

        for (final boolean mrcaHeights : new boolean[]{false, true}) {
            CladeSystem all = new CladeSystem(false, false);
            for (Tree tree : trees) {
                all.add(tree);
            }
            all.calculateCladeCredibilities(trees.size());

            CladeSystem serial = new CladeSystem(target, all, mrcaHeights);
            final AtomicInteger serialVisits = new AtomicInteger();
            for (Tree tree : trees) {
                serial.collectCladeHeights(tree, mrcaHeights);
                serial.traverseTree(tree, counter(serialVisits));
            }
            Map<Object, String> expected = describe(serial);

            if (mrcaHeights) {
                // every clade of the target gets a height from every tree
                for (BiClade clade : serial.getClades()) {
                    assertEquals(trees.size(), clade.getHeightValues().size());
                }
            }

            ExecutorService pool = Executors.newFixedThreadPool(8);
            try {
                for (int repeat = 0; repeat < 25; repeat++) {
                    final CladeSystem threaded = new CladeSystem(target, all, mrcaHeights);
                    final AtomicInteger visits = new AtomicInteger();
                    List<Future<?>> futures = new ArrayList<Future<?>>();
                    for (final Tree tree : trees) {
                        futures.add(pool.submit(() -> {
                            threaded.collectCladeHeights(tree, mrcaHeights);
                            threaded.traverseTree(tree, counter(visits));
                        }));
                    }
                    for (Future<?> future : futures) {
                        future.get();
                    }
                    assertEquals("repeat " + repeat, expected, describe(threaded));
                    assertEquals(serialVisits.get(), visits.get());
                }
            } finally {
                pool.shutdown();
            }
        }
    }

    private static CladeAction counter(final AtomicInteger count) {
        return new CladeAction() {
            @Override
            public void actOnClade(Clade clade, Tree tree, NodeRef node) {
                count.incrementAndGet();
            }

            @Override
            public boolean expectAllClades() {
                // as the actions that collect from every tree do: most clades are not in the target
                return true;
            }
        };
    }

    public void testTaxaAreSetOnce() throws Exception {
        CladeSystem cladeSystem = new CladeSystem(false, false);
        Tree first = tree("((A:1,B:1):1,C:2);");
        cladeSystem.setTaxonList(first);
        cladeSystem.add(first);
        try {
            cladeSystem.setTaxonList(first);
            fail("the taxa cannot be set twice");
        } catch (IllegalStateException e) {
            // expected
        }
        try {
            cladeSystem.add(tree("((A:1,B:1):1,X:2);"));
            fail("a tree with another taxon should be refused");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage().contains("X"));
        }
    }
}
