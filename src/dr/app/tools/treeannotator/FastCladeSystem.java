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
import dr.evolution.util.TaxonList;

import java.util.*;

/**
 * @author Andrew Rambaut
 * @version $
 */
class FastCladeSystem implements CladeSystem {
    //
    // Public stuff
    //


    /**
     *
     */
    public FastCladeSystem() {
    }

    /**
     *
     */
    public FastCladeSystem(Tree targetTree) {
        this.targetTree = targetTree;
        add(targetTree);
    }

    /**
     * adds all the clades in the tree
     */
    @Override
    public void add(Tree tree) {
        if (taxonList == null) {
            taxonList = tree;
        }

        // Recurse over the tree and add all the clades (or increment their
        // frequency if already present). The root clade is added too (for
        // annotation purposes).
        rootClade = addClades(tree, tree.getRoot());

        assert tipClades.size() == tree.getExternalNodeCount();
    }

    @Override
    public Clade getRootClade() {
        return rootClade;
    }

    private Clade addClades(Tree tree, NodeRef node) {
        Clade clade;
        if (tree.isExternal(node)) {

//                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
            int index = node.getNumber();
            clade = new FastClade(index);
//                clade.taxon = tree.getNodeTaxon(node);

        } else {

            assert tree.getChildCount(node) == 2 : "requires a strictly bifurcating tree";

            Clade clade1 = addClades(tree, tree.getChild(node, 0));
            Clade clade2 = addClades(tree, tree.getChild(node, 1));
            clade = new FastClade(clade1, clade2, tree.getExternalNodeCount());
        }

        return addClade(clade);
    }

    private Clade addClade(Clade clade) {
        if (clade.getSize() > 1) {
            Clade c = cladeMap.get(clade.getKey());
            if (c == null) {
                cladeMap.put(clade.getKey(), clade);
                c = clade;
            }
            assert c.getSize() == clade.getSize();
            c.setCount(c.getCount() + 1);
            return c;
        } else {
            tipClades.add(clade);
            return clade;
        }
    }

    @Override
    public void collectAttributes(Set<String> attributeNames, Tree tree) {
        collectAttributes(attributeNames, tree, tree.getRoot());
    }

    private Clade collectAttributes(Set<String> attributeNames, Tree tree, NodeRef node) {

        Clade clade;

        if (tree.isExternal(node)) {

//                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
//                if (index < 0) {
//                    throw new IllegalArgumentException("Taxon, " + tree.getNodeTaxon(node).getId() + ", not found in target tree");
//                }
            int index = node.getNumber();
            clade = new FastClade(index);

        } else {
            assert tree.getChildCount(node) == 2;

            Clade clade1 = collectAttributes(attributeNames, tree, tree.getChild(node, 0));
            Clade clade2 = collectAttributes(attributeNames, tree, tree.getChild(node, 1));
            clade = new FastClade(clade1, clade2, tree.getExternalNodeCount());
        }

        collectAttributesForClade(attributeNames, clade, tree, node);

        return clade;
    }

    private void collectAttributesForClade(Set<String> attributeNames, Clade clade, Tree tree, NodeRef node) {
        if (clade != null) {

            int i = 0;
            Object[] values = new Object[attributeNames.size()];
            for (String attributeName : attributeNames) {
                boolean processed = false;

                if (!processed) {
                    Object value;
                    if (attributeName.equals("height")) {
                        value = tree.getNodeHeight(node);
                    } else if (attributeName.equals("length")) {
                        value = tree.getBranchLength(node);
                    } else {
                        value = tree.getNodeAttribute(node, attributeName);
                        if (value instanceof String && ((String) value).startsWith("\"")) {
                            value = ((String) value).replaceAll("\"", "");
                        }
                    }

                    values[i] = value;
                }
                i++;
            }
            clade.addAttributeValues(values);

            //progressStream.println(clade + " " + clade.getValuesSize());
            clade.setCount(clade.getCount() + 1);
        }
    }

//    public Map<BitSet, Clade> getCladeMap() {
//        return cladeMap;
//    }
//
//    public Clade getClade(BitSet bitSet) {
//        return cladeMap.get(bitSet);
//    }

    @Override
    public void calculateCladeCredibilities(int totalTreesUsed) {
        for (Clade clade : cladeMap.values()) {

            if (clade.getCount() > totalTreesUsed) {

                throw new AssertionError("clade.getCount=(" + clade.getCount() +
                        ") should be <= totalTreesUsed = (" + totalTreesUsed + ")");
            }

            clade.setCredibility(((double) clade.getCount()) / (double) totalTreesUsed);
        }
    }

    @Override
    public double getLogCladeCredibility(Tree tree) {
        double[] logCladeCredibility = { 0.0 };
        logCladeCredibility(tree, tree.getRoot(), logCladeCredibility);
        return logCladeCredibility[0];
    }

    private Clade logCladeCredibility(Tree tree, NodeRef node, double[] logCladeCredibility) {

        Clade clade;

        if (tree.isExternal(node)) {

//                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
            int index = node.getNumber();
            clade = new FastClade(index);
        } else {

            assert tree.getChildCount(node) == 2;

            Clade clade1 = logCladeCredibility(tree, tree.getChild(node, 0), logCladeCredibility);
            Clade clade2 = logCladeCredibility(tree, tree.getChild(node, 1), logCladeCredibility);

            clade = new FastClade(clade1, clade2, tree.getExternalNodeCount());

            logCladeCredibility[0] += Math.log(getCladeCredibility(clade.getKey()));
        }

        return clade;
    }

    private double getCladeCredibility(Object key) {
        Clade clade = cladeMap.get(key);
        assert clade != null;
//        if (clade == null) {
//            return 0.0;
//        }
        return clade.getCredibility();
    }

//    public BitSet removeClades(Tree tree, NodeRef node, boolean includeTips) {
//
//        BitSet bits = new BitSet();
//
//        if (tree.isExternal(node)) {
//
////                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
//            int index = node.getNumber();
//            bits.set(index);
//
//                removeClade(bits);
//
//        } else {
//
//            for (int i = 0; i < tree.getChildCount(node); i++) {
//
//                NodeRef node1 = tree.getChild(node, i);
//
//                bits.or(removeClades(tree, node1, includeTips));
//            }
//
//            removeClade(bits);
//        }
//
//        return bits;
//    }

//    private void removeClade(Clade clade) {
//        clade.setCount(clade.getCount() - 1);
//        if (clade.getCount() == 0) {
//            cladeSet.remove(clade);
//        }
//
//    }

    // Get tree clades as bitSets on target taxa
    // codes is an array of existing BitSet objects, which are reused

    void getTreeCladeCodes(Tree tree, BitSet[] codes) {
        getTreeCladeCodes(tree, tree.getRoot(), codes);
    }

    int getTreeCladeCodes(Tree tree, NodeRef node, BitSet[] codes) {
        final int inode = node.getNumber();
        codes[inode].clear();
        if (tree.isExternal(node)) {
//                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
            int index = node.getNumber();
            codes[inode].set(index);
        } else {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                final NodeRef child = tree.getChild(node, i);
                final int childIndex = getTreeCladeCodes(tree, child, codes);

                codes[inode].or(codes[childIndex]);
            }
        }
        return inode;
    }

    @Override
    public int getCladeCount() {
        return cladeMap.keySet().size();
    }

    //
    // Private stuff
    //
    TaxonList taxonList = null;
    Set<Clade> tipClades = new HashSet<>();
    Map<Object, Clade> cladeMap = new HashMap<>();

    Clade rootClade;

    Tree targetTree;
}
