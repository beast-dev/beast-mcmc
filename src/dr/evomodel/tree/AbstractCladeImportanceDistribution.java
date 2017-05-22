/*
 * AbstractCladeImportanceDistribution.java
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

/**
 *
 */
package dr.evomodel.tree;

import dr.evolution.tree.Clade;
import dr.evolution.tree.ImportanceDistribution;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Likelihood;

import java.util.BitSet;
import java.util.HashMap;
import java.util.List;

/**
 * @author shhn001
 */
public abstract class AbstractCladeImportanceDistribution implements
        ImportanceDistribution {

    /**
     *
     */
    public AbstractCladeImportanceDistribution() {
        // TODO Auto-generated constructor stub
    }

    /**
     * creates a list with all clades but just the non-complementary ones.
     * ((A,B),(C,D)) just {A,B,C,D} and {A,B} are inserted. {A,B} is
     * complementary to {C,D}
     *
     * @param tree - the tree from which the clades are extracted
     * @param node - the starting node. All clades below starting at this branch
     *             are added
     */
    protected Clade getNonComplementaryClades(Tree tree, NodeRef node,
                                              List<Clade> parentClades, List<Clade> childClade,
                                              HashMap<String, Integer> taxonMap) {

        // create a new bit set for this clade
        BitSet bits = new BitSet();
        Clade c = null;

        // check if the node is external
        if (tree.isExternal(node)) {

            // if so, the only taxon in the clade is I
            // int index = node.getNumber();
            String taxon = tree.getTaxon(node.getNumber()).getId();
            int index = taxonMap.get(taxon);
            bits.set(index);

            c = new Clade(bits, tree.getNodeHeight(node));

        } else {
            // otherwise, call all children and add its taxon together to one
            // clade

            NodeRef childNode = tree.getChild(node, 0);
            // add just my first child to the list
            // the second child is complementary to the first
            Clade leftChild = getNonComplementaryClades(tree, childNode,
                    parentClades, childClade, taxonMap);
            bits.or(leftChild.getBits());

            childNode = tree.getChild(node, 1);
            // add just my first child to the list
            // the second child is complementary to the first
            Clade rightChild = getNonComplementaryClades(tree, childNode,
                    parentClades, childClade, taxonMap);
            bits.or(rightChild.getBits());

            c = new Clade(bits, tree.getNodeHeight(node));

            if (leftChild.getSize() >= 2) {
                parentClades.add(c);
                childClade.add(leftChild);
            } else if (rightChild.getSize() >= 2) {
                parentClades.add(c);
                childClade.add(rightChild);
            }
        }

        return c;
    }

    /**
     * creates a list with all clades but just the non-complementary ones.
     * ((A,B),(C,D)) just {A,B,C,D} and {A,B} are inserted. {A,B} is
     * complementary to {C,D}
     *
     * @param tree - the tree from which the clades are extracted
     * @param node - the starting node. All clades below starting at this branch
     *             are added
     */
    protected Clade getNonComplementaryClades(Tree tree, NodeRef node,
                                              List<Clade> parentClades, List<Clade> childClade) {

        // create a new bit set for this clade
        BitSet bits = new BitSet();
        Clade c = null;

        // check if the node is external
        if (tree.isExternal(node)) {

            // if so, the only taxon in the clade is I
            int index = node.getNumber();
            bits.set(index);

            c = new Clade(bits, tree.getNodeHeight(node));

        } else {
            // otherwise, call all children and add its taxon together to one
            // clade

            NodeRef childNode = tree.getChild(node, 0);
            // add just my first child to the list
            // the second child is complementary to the first
            Clade leftChild = getNonComplementaryClades(tree, childNode,
                    parentClades, childClade);
            bits.or(leftChild.getBits());

            childNode = tree.getChild(node, 1);
            // add just my first child to the list
            // the second child is complementary to the first
            Clade rightChild = getNonComplementaryClades(tree, childNode,
                    parentClades, childClade);
            bits.or(rightChild.getBits());

            c = new Clade(bits, tree.getNodeHeight(node));

            if (leftChild.getSize() >= 2) {
                parentClades.add(c);
                childClade.add(leftChild);
            } else if (rightChild.getSize() >= 2) {
                parentClades.add(c);
                childClade.add(rightChild);
            }
        }

        return c;
    }

    /**
     * Creates a list with all clades of the tree
     *
     * @param tree - the tree from which the clades are extracted
     * @param node - the starting node. All clades below starting at this branch
     *             are added
     */
    protected Clade getClades(Tree tree, NodeRef node,
                              List<Clade> parentClades, List<Clade> childClade) {

        // create a new bit set for this clade
        BitSet bits = new BitSet();
        Clade c = null;

        // check if the node is external
        if (tree.isExternal(node)) {

            // if so, the only taxon in the clade is I
            int index = node.getNumber();
            bits.set(index);

            c = new Clade(bits, tree.getNodeHeight(node));

        } else {
            // otherwise, call all children and add its taxon together to one
            // clade

            NodeRef childNode = tree.getChild(node, 0);
            // add just my first child to the list
            // the second child is complementary to the first
            Clade leftChild = getClades(tree, childNode, parentClades,
                    childClade);
            bits.or(leftChild.getBits());

            childNode = tree.getChild(node, 1);
            // add just my first child to the list
            // the second child is complementary to the first
            Clade rightChild = getClades(tree, childNode, parentClades,
                    childClade);
            bits.or(rightChild.getBits());

            c = new Clade(bits, tree.getNodeHeight(node));

            if (leftChild.getSize() >= 2) {
                parentClades.add(c);
                childClade.add(leftChild);
            }
            if (rightChild.getSize() >= 2) {
                parentClades.add(c);
                childClade.add(rightChild);
            }
        }

        return c;
    }

    /**
     * Creates a list with all clades of the tree
     *
     * @param taxonMap - the lookup map for the taxon representing an index
     * @param tree     - the tree from which the clades are extracted
     * @param node     - the starting node. All clades below starting at this branch
     *                 are added
     */
    protected Clade getClades(Tree tree, NodeRef node,
                              List<Clade> parentClades, List<Clade> childClade,
                              HashMap<String, Integer> taxonMap) {

        // create a new bit set for this clade
        BitSet bits = new BitSet();
        Clade c = null;

        // check if the node is external
        if (tree.isExternal(node)) {

            // if so, the only taxon in the clade is I
            // int index = node.getNumber();
            String taxon = tree.getTaxon(node.getNumber()).getId();
            int index = taxonMap.get(taxon);
            bits.set(index);

            c = new Clade(bits, tree.getNodeHeight(node));

        } else {
            // otherwise, call all children and add its taxon together to one
            // clade

            NodeRef childNode = tree.getChild(node, 0);
            // add just my first child to the list
            // the second child is complementary to the first
            Clade leftChild = getClades(tree, childNode, parentClades,
                    childClade, taxonMap);
            bits.or(leftChild.getBits());

            childNode = tree.getChild(node, 1);
            // add just my first child to the list
            // the second child is complementary to the first
            Clade rightChild = getClades(tree, childNode, parentClades,
                    childClade, taxonMap);
            bits.or(rightChild.getBits());

            c = new Clade(bits, tree.getNodeHeight(node));

            if (leftChild.getSize() >= 2) {
                parentClades.add(c);
                childClade.add(leftChild);
            }
            if (rightChild.getSize() >= 2) {
                parentClades.add(c);
                childClade.add(rightChild);
            }
        }

        return c;
    }

    /**
     * Creates a list with all clades of the tree
     *
     * @param tree   - the tree from which the clades are extracted
     * @param node   - the starting node. All clades below starting at this branch
     *               are added
     * @param clades - the list in which the clades are stored
     * @param bits   - a bit set to which the current bits of the clades are added
     */
    public void getClades(Tree tree, NodeRef node, List<Clade> clades,
                          BitSet bits) {

        // create a new bit set for this clade
        BitSet bits2 = new BitSet();

        // check if the node is external
        if (tree.isExternal(node)) {

            // if so, the only taxon in the clade is I
            int index = node.getNumber();
            bits2.set(index);

        } else {

            // otherwise, call all children and add its taxon together to one
            // clade
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                getClades(tree, child, clades, bits2);
            }
            // add my bit set to the list
            clades.add(new Clade(bits2, tree.getNodeHeight(node)));
        }

        // add my bit set to the bit set I was given
        // this is needed for adding all children clades together
        if (bits != null) {
            bits.or(bits2);
        }
    }

    /**
     * Creates a list with all clades of the tree
     *
     * @param tree    - the tree from which the clades are extracted
     * @param node    - the starting node. All clades heights below starting at this
     *                branch are added
     * @param heights - the list in which the heights are stored
     */
    public void getCladesHeights(Tree tree, NodeRef node, List<Double> heights) {

        // check if the node is external
        if (tree.isExternal(node)) {

            // if so, do nothing

        } else {

            // otherwise, call all children and add its taxon together to one
            // clade
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                getCladesHeights(tree, child, heights);
            }
            // add my bit set to the list
            heights.add(tree.getNodeHeight(node));
        }
    }

    /**
     * Creates a list with all clades of the tree
     *
     * @param tree    - the tree from which the clades are extracted
     * @param node    - the starting node. All clades heights below starting at this
     *                branch are added
     * @param heights - the list in which the heights are stored
     */
    public void getRelativeCladesHeights(Tree tree, NodeRef node,
                                         List<Double> heights) {

        // check if the node is external
        if (tree.isExternal(node)) {

            // if so, do nothing

        } else {

            // otherwise, call all children and add its taxon together to one
            // clade
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                getRelativeCladesHeights(tree, child, heights);
            }
            // add my bit set to the list
            if (node != tree.getRoot()) {
                NodeRef parent = tree.getParent(node);
                heights.add(tree.getNodeHeight(node)
                        / tree.getNodeHeight(parent));
            } else {
                heights.add(1.0);
            }
        }
    }

    /**
     * Creates a list with all clades of the tree
     *
     * @param tree - the tree from which the clades are extracted
     * @param node - the starting node. All clades below starting at this branch
     *             are added
     */
    protected Clade getClade(Tree tree, NodeRef node) {

        // create a new bit set for this clade
        BitSet bits = new BitSet();

        // check if the node is external
        if (tree.isExternal(node)) {

            // if so, the only taxon in the clade am I
            int index = node.getNumber();
            bits.set(index);

        } else {

            // otherwise, call all children and add its taxon together to one
            // clade
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                Clade c = getClade(tree, child);
                bits.or(c.getBits());
            }
        }

        Clade c = new Clade(bits, tree.getNodeHeight(node));

        return c;
    }

    /**
     * Finds the parent of a given clade in a list of clades. The parent is the
     * direct parent and not the grandparent or so.
     *
     * @param clades - list of clades in which we are searching the parent
     * @param child  - the child of whom we are searching the parent
     * @return the parent clade if found, otherwise itself
     */
    protected Clade getParentClade(List<Clade> clades, Clade child) {
        Clade parent = null;
        BitSet childBits = child.getBits();
        int parentSize = Integer.MAX_VALUE;

        // look in all clades of the list which contains the child and has the
        // minimum cardinality (least taxa) -> that's the parent :-)
        for (int i = 0; i < clades.size(); i++) {
            Clade tmp = clades.get(i);
            if (!child.equals(tmp) && containsClade(tmp.getBits(), childBits)) {
                if (parent == null || parentSize > tmp.getSize()) {
                    parent = tmp;
                    parentSize = parent.getSize();
                }
            }
        }
        // if there isn't a parent, then you probably asked for the whole tree
        if (parent == null) {
            parent = child;
        }

        return parent;
    }

    /**
     * Checks if clade i contains clade j.
     *
     * @param i - the parent clade
     * @param j - the child clade
     * @return true, if i contains j
     */
    protected boolean containsClade(Clade i, Clade j) {
        return containsClade(i.getBits(), j.getBits());
    }

    /**
     * Checks if clade i contains clade j.
     *
     * @param i - the parent clade
     * @param j - the child clade
     * @return true, if i contains j
     */
    protected boolean containsClade(BitSet i, BitSet j) {
        BitSet tmpI = (BitSet) i.clone();

        // just set the bits which are either in j but not in i or in i but not
        // in j
        tmpI.xor(j);
        int numberOfBitsInEither = tmpI.cardinality();
        // which bits are just in i
        tmpI.and(i);
        int numberIfBitJustInContaining = tmpI.cardinality();

        // if the number of bits just in i is equal to the number of bits just
        // in one of i or j
        // then i contains j
        return numberOfBitsInEither == numberIfBitJustInContaining;
    }

    /*
      * (non-Javadoc)
      *
      * @see
      * dr.evolution.tree.ImportanceDistribution#addTree(dr.evolution.tree.Tree)
      */
    public abstract void addTree(Tree tree);

    /*
      * (non-Javadoc)
      *
      * @see
      * dr.evolution.tree.ImportanceDistribution#getTreeProbability(dr.evolution
      * .tree.Tree)
      */
    public abstract double getTreeProbability(Tree tree);

    /*
      * (non-Javadoc)
      *
      * @see
      * dr.evolution.tree.ImportanceDistribution#splitClade(dr.evolution.tree
      * .Clade, dr.evolution.tree.Clade[])
      */
    public abstract double splitClade(Clade parent, Clade[] children);

    public abstract double setNodeHeights(TreeModel tree, Likelihood likelihood);

    public abstract double getChanceForNodeHeights(TreeModel tree, Likelihood likelihood);

}
