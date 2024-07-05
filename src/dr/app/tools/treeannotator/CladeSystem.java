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
import dr.util.Pair;

import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
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
        add(targetTree);
    }

    /**
     * adds all the clades in the tree
     */
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

    public Clade getRootClade() {
        return rootClade;
    }

    private Clade addClades(Tree tree, NodeRef node) {
        Clade clade;
        if (tree.isExternal(node)) {

//                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
            int index = node.getNumber();
            clade = new Clade(index);
//                clade.taxon = tree.getNodeTaxon(node);

        } else {

            assert tree.getChildCount(node) == 2 : "requires a strictly bifurcating tree";

            Clade clade1 = addClades(tree, tree.getChild(node, 0));
            Clade clade2 = addClades(tree, tree.getChild(node, 1));
            clade = new Clade(clade1, clade2, tree.getExternalNodeCount());
        }

        return addClade(clade);
    }

    private Clade addClade(Clade clade) {
        if (clade.bits != null) {
            Clade c = cladeMap.get(clade.hash);
            if (c == null) {
                cladeMap.put(clade.hash, clade);
                c = clade;
            }
            assert c.size == clade.size;
            c.setCount(c.getCount() + 1);
            return c;
        } else {
            tipClades.add(clade);
            return clade;
        }
    }

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
            clade = new Clade(index);

        } else {
            assert tree.getChildCount(node) == 2;

            Clade clade1 = collectAttributes(attributeNames, tree, tree.getChild(node, 0));
            Clade clade2 = collectAttributes(attributeNames, tree, tree.getChild(node, 1));
            clade = new Clade(clade1, clade2, tree.getExternalNodeCount());
        }

        collectAttributesForClade(attributeNames, clade, tree, node);

        return clade;
    }

    private void collectAttributesForClade(Set<String> attributeNames, Clade clade, Tree tree, NodeRef node) {
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
// AR - we deal with this once everything
//                        } else if (attributeName.equals(location1Attribute)) {
//                            // If this is one of the two specified bivariate location names then
//                            // merge this and the other one into a single array.
//                            Object value1 = tree.getNodeAttribute(node, attributeName);
//                            Object value2 = tree.getNodeAttribute(node, location2Attribute);
//
//                            value = new Object[]{value1, value2};
//                        } else if (attributeName.equals(location2Attribute)) {
//                            // do nothing - already dealt with this...
//                            value = null;
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

//    public Map<BitSet, Clade> getCladeMap() {
//        return cladeMap;
//    }
//
//    public Clade getClade(BitSet bitSet) {
//        return cladeMap.get(bitSet);
//    }

    public void calculateCladeCredibilities(int totalTreesUsed) {
        for (Clade clade : cladeMap.values()) {

            if (clade.getCount() > totalTreesUsed) {

                throw new AssertionError("clade.getCount=(" + clade.getCount() +
                        ") should be <= totalTreesUsed = (" + totalTreesUsed + ")");
            }

            clade.setCredibility(((double) clade.getCount()) / (double) totalTreesUsed);
        }
    }

    public Clade getLogCladeCredibility(Tree tree, NodeRef node, double[] logCladeCredibility) {

        Clade clade;

        if (tree.isExternal(node)) {

//                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
            int index = node.getNumber();
            clade = new Clade(index);
        } else {

            assert tree.getChildCount(node) == 2;

            Clade clade1 = getLogCladeCredibility(tree, tree.getChild(node, 0), logCladeCredibility);
            Clade clade2 = getLogCladeCredibility(tree, tree.getChild(node, 1), logCladeCredibility);

            clade = new Clade(clade1, clade2, tree.getExternalNodeCount());

            logCladeCredibility[0] += Math.log(getCladeCredibility(clade.hash));
        }

        return clade;
    }

    private double getCladeCredibility(int hash) {
        Clade clade = cladeMap.get(hash);
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

    public int getCladeCount() {
        return cladeMap.keySet().size();
    }

    static class Clade {
        public Clade(int index) {
            this.index = index;

            count = 0;
            credibility = 1.0;
            size = 1;
            bits = null;
            hash = index;
        }

        public Clade(Clade subClade1, Clade subClade2, int tipCount) {
            count = 0;
            credibility = 0.0;
            size = subClade1.size + subClade2.size;
            index = -1;

            bits = new byte[(tipCount / 8) + 1] ;
            if (subClade1.bits == null) {
                int byteIndex = subClade1.index / 8;
                int bitMask = 1 << (subClade1.index % 8);
                bits[byteIndex] = (byte) bitMask;
            } else {
                System.arraycopy(subClade1.bits, 0, bits, 0, subClade1.bits.length);
            }

            if (subClade2.bits == null) {
                int byteIndex = subClade2.index / 8;
                int bitMask = 1 << (subClade2.index % 8);
                bits[byteIndex] |= (byte) bitMask;
            } else {
                for (int i = 0; i < bits.length; i++) {
                    bits[i] |= subClade2.bits[i];
                }
            }

//            Clade sc1, sc2;
//            if (subClade1.index < subClade2.index) {
//                index = subClade1.index;
//                sc1 = subClade1;
//                sc2 = subClade2;
//            } else {
//                index = subClade2.index;
//                sc1 = subClade2;
//                sc2 = subClade1;
//            }
//
//            int maxIndex = Math.max(sc1.maxIndex(), sc2.maxIndex());
//            bits = new byte[maxIndex - index + 1];
//
//            if (sc1.bits == null) {
//                bits[0] = 0b10000000;
//            } else {
//                System.arraycopy(sc1.bits, 0, bits, 0, sc1.bits.length);
//            }
//
//            if (sc2.bits == null) {
//                bits[sc2.index - index] = 1;
//            } else {
//                for (int i = 0; i < sc2.bits.length; i++) {
//                    bits[i + sc2.index - index] |= sc2.bits[i];
//                }
//            }

            MessageDigest digest = null;
            try {
                digest = MessageDigest.getInstance("SHA-256");
                byte[] encodedhash = digest.digest(bits);
                hash = Arrays.hashCode(encodedhash);
            } catch (NoSuchAlgorithmException e) {
                throw new RuntimeException(e);
            }

        }

        private int maxIndex() {
            if (bits != null) {
                return (index + bits.length - 1);
            } else {
                return index;
            }
        }

        private int byteIndex(int index) {
            return index >> 8;
        }
        private int bitIndex(int index) {
            return index | 8;
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

        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            if (((Clade) o).size != size) return false;

            return !(bits != null ? !Arrays.equals(bits, ((Clade) o).bits) : ((Clade) o).bits != null);

        }

        public int hashCode() {
            return hash;
        }

        public String toString() {
            return "clade " + hashCode();
        }

        int count;
        double credibility;
        final int size;
        final byte[] bits;
        final int index;

        final int hash;
        List<Object[]> attributeValues = null;
        Set<Pair<Clade, Clade>> subClades = null;
        Clade bestLeft = null;
        Clade bestRight = null;
        double bestSubTreeCredibility;
    }

    //
    // Private stuff
    //
    TaxonList taxonList = null;
    Set<Clade> tipClades = new HashSet<>();
    Map<Integer, Clade> cladeMap = new HashMap<>();

    Clade rootClade;

    Tree targetTree;
}
