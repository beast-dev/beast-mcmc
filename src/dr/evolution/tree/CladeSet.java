/*
 * CladeSet.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evolution.tree;

import dr.evolution.util.TaxonList;
import dr.util.FrequencySet;

import java.util.*;

/**
 * Stores a set of unique clades (and their node heights) for a tree
 *
 * @version $Id: CladeSet.java,v 1.8 2005/05/24 20:25:56 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class CladeSet extends FrequencySet<BitSet> {
    //
    // Public stuff
    //

    public CladeSet() {}

    /**
     * @param tree
     */
    public CladeSet(Tree tree)
    {
        this(tree, tree);
    }

    /**
     * @param taxonList  a set of taxa used to label the tips
     */
    public CladeSet(Tree tree, TaxonList taxonList)
    {
        this.taxonList = taxonList;
        add(tree);
    }

    /** get number of unique clades */
    public int getCladeCount()
    {
        return size();
    }

    /** get clade bit set */
    public String getClade(int index)
    {
        BitSet bits = get(index);

        StringBuffer buffer = new StringBuffer("{");
        boolean first = true;
        for (String taxonId : getTaxaSet(bits)) {
            if (!first) {
                buffer.append(", ");
            } else {
                first = false;
            }
            buffer.append(taxonId);
        }
        buffer.append("}");
        return buffer.toString();
    }

    private SortedSet<String> getTaxaSet(BitSet bits) {

        SortedSet<String> taxaSet = new TreeSet<String>();

        for (int i = 0; i < bits.length(); i++) {
            if (bits.get(i)) {
                taxaSet.add(taxonList.getTaxonId(i));
            }
        }
        return taxaSet;
    }


    /** get clade frequency */
    int getCladeFrequency(int index)
    {
        return getFrequency(index);
    }

    /** adds all the clades in the tree */
    public void add(Tree tree)
    {
        if (taxonList == null) {
            taxonList = tree;
        }

        totalTrees += 1;

        // Recurse over the tree and add all the clades (or increment their
        // frequency if already present). The root clade is not added.
        addClades(tree, tree.getRoot(), null);
    }

    private void addClades(Tree tree, NodeRef node, BitSet bits) {

        if (tree.isExternal(node)) {

            if (taxonList != null) {
                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
                bits.set(index);
            } else {
                bits.set(node.getNumber());
            }
        } else {

            BitSet bits2 = new BitSet();
            for (int i = 0; i < tree.getChildCount(node); i++) {

                NodeRef node1 = tree.getChild(node, i);

                addClades(tree, node1, bits2);
            }

            add(bits2, 1);
            addNodeHeight(bits2, tree.getNodeHeight(node));

            if (bits != null) {
                bits.or(bits2);
            }
        }
    }

    public double getMeanNodeHeight(int i) {
        BitSet bits = get(i);

        return getTotalNodeHeight(bits) / getFrequency(i);
    }

    private double getTotalNodeHeight(BitSet bits) {
        Double tnh = totalNodeHeight.get(bits);
        if (tnh == null) return 0.0;
        return tnh;
    }

    private void addNodeHeight(BitSet bits, double height) {
        totalNodeHeight.put(bits, (getTotalNodeHeight(bits) + height));
    }

    // Generifying found that this code was buggy. Kuckily it is not used anymore.

//    /** adds all the clades in the CladeSet */
//    public void add(CladeSet cladeSet)
//    {
//        for (int i = 0, n = cladeSet.getCladeCount(); i < n; i++) {
//            add(cladeSet.getClade(i), cladeSet.getCladeFrequency(i));
//        }
//    }

    private BitSet annotate(MutableTree tree, NodeRef node, String freqAttrName) {
        BitSet b = null;
        if (tree.isExternal(node)) {
            int index;
            if (taxonList != null) {
                index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
            } else {
                index = node.getNumber();
            }
            b = new BitSet(tree.getExternalNodeCount());
            b.set(index);

        } else {

            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                BitSet b1 = annotate(tree, child, freqAttrName);
                if( i == 0 ) {
                    b = b1;
                } else {
                    b.or(b1);
                }
            }
            final int total = getFrequency(b);
            if( total >= 0 ) {
                tree.setNodeAttribute(node, freqAttrName, total / (double)totalTrees );
            }
        }
        return b;
    }

    /**
     * Annotate clades of tree with posterior probability
     * @param tree
     * @param freqAttrName name of attribute to set per node
     * @return sum(log(all clades probability))
     */
    public double annotate(MutableTree tree, String freqAttrName) {
        annotate(tree, tree.getRoot(), freqAttrName);

        double logClade = 0.0;
        for(int n = 0; n < tree.getInternalNodeCount(); ++n) {
            final double f = (Double)tree.getNodeAttribute(tree.getInternalNode(n), freqAttrName);
            logClade += Math.log(f);
        }
        return logClade;
    }

    public boolean hasClade(int index, Tree tree) {
        BitSet bits = get(index);

        NodeRef[] mrca = new NodeRef[1];
        findClade(bits, tree, tree.getRoot(), mrca);

        return (mrca[0] != null);
    }

    private int findClade(BitSet bitSet, Tree tree, NodeRef node, NodeRef[] cladeMRCA) {

        if (tree.isExternal(node)) {

            if (taxonList != null) {
                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
                if (bitSet.get(index)) return 1;
            } else {
                if (bitSet.get(node.getNumber())) return 1;
            }
            return -1;
        } else {
            int count = 0;
            for (int i = 0; i < tree.getChildCount(node); i++) {

                NodeRef node1 = tree.getChild(node, i);

                int childCount = findClade(bitSet, tree, node1, cladeMRCA);

                if (childCount != -1 && count != -1) {
                    count += childCount;
                } else count = -1;
            }

            if (count == bitSet.cardinality()) cladeMRCA[0] = node;

            return count;
        }
    }

    //
    // Private stuff
    //
    private TaxonList taxonList = null;
    private final Map<BitSet, Double> totalNodeHeight = new HashMap<BitSet, Double>();
    private int totalTrees = 0;
}
