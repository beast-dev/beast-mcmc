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
import dr.stats.DiscreteStatistics;

import java.util.*;

/**
 * @author Andrew Rambaut
 * @version $
 */
public final class CladeSystem {
    private final boolean keepSubClades;
    private final boolean keepParents;
    private double treeCount = 0;
    private boolean storeTipHeights = false;

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
    public CladeSystem(Tree targetTree) {
        this.keepSubClades = false;
        this.keepParents = false;
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
            BiClade clade = new BiClade(index, taxon);
            tipClades.put(index, clade);
        }
    }

    /**
     * recursively add all the clades in a tree
     */
    private BiClade addClades(Tree tree, NodeRef node) {
        BiClade clade;
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

            BiClade clade1 = addClades(tree, tree.getChild(node, 0));
            BiClade clade2 = addClades(tree, tree.getChild(node, 1));
            synchronized (cladeMap) {
                clade = getOrAddClade(clade1, clade2);
            }

            if (keepParents) {
                clade1.addParent(clade);
                clade2.addParent(clade);
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
    private BiClade getOrAddClade(Clade child1, Clade child2) {
        Object key = BiClade.makeKey(child1.getKey(), child2.getKey());
        BiClade clade = cladeMap.get(key);
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
//            }
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

    public void traverseNonBinaryTree(Tree tree, CladeAction action) {
        traverseNonBinaryTree(tree, tree.getRoot(), action);
    }

    private Object traverseNonBinaryTree(Tree tree, NodeRef node, CladeAction action) {

        Object key;

        if (tree.isExternal(node)) {
//            key = node.getNumber();
            key = taxonNumberMap.get(tree.getNodeTaxon(node));
        } else {
            List<Object> keys = new ArrayList<>();
            for (int i = 0; i < tree.getChildCount(node); i++) {
                keys.add(traverseNonBinaryTree(tree, tree.getChild(node, i), action));
            }
            key = BiClade.makeKey(keys.toArray());
        }

        Clade clade = getClade(key);
        if (clade != null) {
            action.actOnClade(clade, tree, node);
        } else {
            assert action.expectAllClades();
        }

        return key;
    }

    public void collectCladeHeights(Tree tree) {
        collectCladeHeights(tree, tree.getRoot());
    }

    private Object collectCladeHeights(Tree tree, NodeRef node) {

        Object key;

        if (tree.isExternal(node)) {
            key = node.getNumber();
            if (taxonNumberMap != null) {
                key = taxonNumberMap.get(tree.getNodeTaxon(node));
            }

            if (storeTipHeights) {
                BiClade tip = (BiClade) getClade(key);
                tip.addHeightValue(tree.getNodeHeight(node));
            }
        } else {
            assert tree.getChildCount(node) == 2;

            Object key1 = collectCladeHeights(tree, tree.getChild(node, 0));
            Object key2 = collectCladeHeights(tree, tree.getChild(node, 1));

            Clade child1 = getClade(key1);
            Clade child2 = getClade(key2);

            key = BiClade.makeKey(key1, key2);

            BiClade clade = (BiClade)getClade(key);

            if (clade.getBestLeft() == child1 && clade.getBestRight() == child2) {
                clade.addChildHeightValues(tree.getNodeHeight(tree.getChild(node, 0)), tree.getNodeHeight(tree.getChild(node, 1)));
                clade.addHeightValue(tree.getNodeHeight(node));
            } else if (clade.getBestLeft() == child2 && clade.getBestRight() == child1) {
                clade.addChildHeightValues(tree.getNodeHeight(tree.getChild(node, 1)), tree.getNodeHeight(tree.getChild(node, 0)));
                clade.addHeightValue(tree.getNodeHeight(node));
            } else {
                clade.addHeightValue(tree.getNodeHeight(node));
            }
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
                    clades.add((BiClade)clade);
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

    //
    // Private stuff
    //
    TaxonList taxonList = null;
    private final Map<Taxon, Integer> taxonNumberMap = new HashMap<>();

    private final Map<Object, BiClade> tipClades = new HashMap<>();
    private final Map<Object, BiClade> cladeMap = new HashMap<>();

    Clade rootClade;

}
