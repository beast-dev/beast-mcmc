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

package dr.app.tools.oldtreeannotator;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.util.Pair;

import java.util.*;

/**
 * @author Andrew Rambaut
 * @version $
 */
class CladeSystem {
    //
    // Public stuff
    //


    /**
     *
     */
    public CladeSystem() {
    }

    /**
     *
     */
    public CladeSystem(Tree targetTree) {
        this.targetTree = targetTree;
        add(targetTree, true);
    }

    /**
     * adds all the clades in the tree
     */
    public void add(Tree tree, boolean includeTips) {
        if (taxonList == null) {
            setTaxonList(tree);
        }

        // Recurse over the tree and add all the clades (or increment their
        // frequency if already present). The root clade is added too (for
        // annotation purposes).
        BitSet rootBits = addClades(tree, tree.getRoot(), includeTips);
        rootClade = cladeMap.get(rootBits);
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

    private BitSet addClades(Tree tree, NodeRef node, boolean includeTips) {

        BitSet bits = new BitSet();

        if (tree.isExternal(node)) {

//                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
            int index = node.getNumber();
            if (taxonNumberMap != null) {
                index = taxonNumberMap.get(tree.getNodeTaxon(node));
            }
            bits.set(index);

            if (includeTips) {
                Clade clade = addClade(bits);
                clade.taxon = tree.getNodeTaxon(node);
            }

        } else {

            List<BitSet> subClades = new ArrayList<BitSet>();

            for (int i = 0; i < tree.getChildCount(node); i++) {

                NodeRef node1 = tree.getChild(node, i);
                BitSet subClade = addClades(tree, node1, includeTips);
                bits.or(subClade);
                subClades.add(subClade);
            }

            Clade clade = addClade(bits);

            if (subClades.size() != 2) {
                throw new IllegalArgumentException("TreeAnnotator requires strictly bifurcating trees");
            }
            clade.addSubclades(subClades.get(0), subClades.get(1));
        }

        return bits;
    }

    private Clade addClade(BitSet bits) {
        Clade clade = cladeMap.get(bits);
        if (clade == null) {
            clade = new Clade(bits);
            cladeMap.put(bits, clade);
        }
        clade.setCount(clade.getCount() + 1);

        return clade;
    }

    public void collectAttributes(Set<String> attributeNames, Tree tree) {
        collectAttributes(attributeNames, tree, tree.getRoot());
    }

    private BitSet collectAttributes(Set<String> attributeNames, Tree tree, NodeRef node) {

        BitSet bits = new BitSet();

        if (tree.isExternal(node)) {

//                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
//                if (index < 0) {
//                    throw new IllegalArgumentException("Taxon, " + tree.getNodeTaxon(node).getId() + ", not found in target tree");
//                }
            int index = node.getNumber();
            if (taxonNumberMap != null) {
                index = taxonNumberMap.get(tree.getNodeTaxon(node));
            }
            bits.set(index);

        } else {

            for (int i = 0; i < tree.getChildCount(node); i++) {

                NodeRef node1 = tree.getChild(node, i);

                bits.or(collectAttributes(attributeNames, tree, node1));
            }
        }

        collectAttributesForClade(attributeNames, bits, tree, node);

        return bits;
    }

    private void collectAttributesForClade(Set<String> attributeNames, BitSet bits, Tree tree, NodeRef node) {
        Clade clade = cladeMap.get(bits);
        if (clade != null) {

            if (clade.attributeValues == null) {
                clade.attributeValues = new ArrayList<Object[]>();
            }

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

                    //if (value == null) {
                    //    progressStream.println("attribute " + attributeNames[i] + " is null.");
                    //}

                    values[i] = value;
                }
                i++;
            }
            clade.attributeValues.add(values);

            //progressStream.println(clade + " " + clade.getValuesSize());
            clade.setCount(clade.getCount() + 1);
        }
    }

    public Map<BitSet, Clade> getCladeMap() {
        return cladeMap;
    }

    public Clade getClade(BitSet bitSet) {
        return cladeMap.get(bitSet);
    }

    public void calculateCladeCredibilities(int totalTreesUsed) {
        for (Clade clade : cladeMap.values()) {

            if (clade.getCount() > totalTreesUsed) {

                throw new AssertionError("clade.getCount=(" + clade.getCount() +
                        ") should be <= totalTreesUsed = (" + totalTreesUsed + ")");
            }

            clade.setCredibility(((double) clade.getCount()) / (double) totalTreesUsed);
        }
    }

    public double getLogCladeCredibility(Tree tree, NodeRef node, BitSet bits) {

        double logCladeCredibility = 0.0;

        if (tree.isExternal(node)) {

//                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
            int index = node.getNumber();
            if (taxonNumberMap != null) {
                index = taxonNumberMap.get(tree.getNodeTaxon(node));
            }
            bits.set(index);
        } else {

            BitSet bits2 = new BitSet();
            for (int i = 0; i < tree.getChildCount(node); i++) {

                NodeRef node1 = tree.getChild(node, i);

                logCladeCredibility += getLogCladeCredibility(tree, node1, bits2);
            }

            logCladeCredibility += Math.log(getCladeCredibility(bits2));

            if (bits != null) {
                bits.or(bits2);
            }
        }

        return logCladeCredibility;
    }

    private double getCladeCredibility(BitSet bits) {
        Clade clade = cladeMap.get(bits);
        if (clade == null) {
            return 0.0;
        }
        return clade.getCredibility();
    }

    public BitSet removeClades(Tree tree, NodeRef node, boolean includeTips) {

        BitSet bits = new BitSet();

        if (tree.isExternal(node)) {

//                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
            int index = node.getNumber();
            if (taxonNumberMap != null) {
                index = taxonNumberMap.get(tree.getNodeTaxon(node));
            }
            bits.set(index);

            if (includeTips) {
                removeClade(bits);
            }

        } else {

            for (int i = 0; i < tree.getChildCount(node); i++) {

                NodeRef node1 = tree.getChild(node, i);

                bits.or(removeClades(tree, node1, includeTips));
            }

            removeClade(bits);
        }

        return bits;
    }

    private void removeClade(BitSet bits) {
        Clade clade = cladeMap.get(bits);
        if (clade != null) {
            clade.setCount(clade.getCount() - 1);
        }

    }

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
            if (taxonNumberMap != null) {
                index = taxonNumberMap.get(tree.getNodeTaxon(node));
            }
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

    class Clade {
        public Clade(BitSet bits) {
            this.bits = bits;
            count = 0;
            credibility = 0.0;
            size = bits.cardinality();
        }

        public int getCount() {
            return count;
        }

        public void setCount(int count) {
            this.count = count;
        }

        public double getCredibility() {
            return credibility;
        }

        public void setCredibility(double credibility) {
            this.credibility = credibility;
        }

        public void addSubclades(BitSet subClade1, BitSet subClade2) {
            if (this.subClades == null) {
                this.subClades = new HashSet<>();
            }
            // Store the subclade with lowest first set bit index as the first of the pair to make
            // sure the order is the same if the pair is the same.
            if (subClade1.nextSetBit(0) < subClade2.nextSetBit(0)) {
                this.subClades.add(new Pair<>(subClade1, subClade2));
            } else {
                this.subClades.add(new Pair<>(subClade2, subClade1));
            }
        }

        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final Clade clade = (Clade) o;

            return !(bits != null ? !bits.equals(clade.bits) : clade.bits != null);

        }

        public int hashCode() {
            return (bits != null ? bits.hashCode() : 0);
        }

        public String toString() {
            return "clade " + bits.toString();
        }

        int count;
        double credibility;
        final int size;
        final BitSet bits;
        Taxon taxon = null;
        List<Object[]> attributeValues = null;
        Set<Pair<BitSet, BitSet>> subClades = null;
        Clade bestLeft = null;
        Clade bestRight = null;
        double bestSubTreeCredibility;
    }

    //
    // Private stuff
    //
    TaxonList taxonList = null;
    Map<Taxon, Integer> taxonNumberMap = null;
    Map<BitSet, Clade> cladeMap = new HashMap<>();

    Clade rootClade;

    Tree targetTree;
}
