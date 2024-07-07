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
final class CladeSystem {
    /**
     * Constructor starting with an empty clade system
     */
    public CladeSystem() {
    }

    /**
     * Constructor adding a single target tree
     */
    public CladeSystem(Tree targetTree) {
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

        assert rootClade.getSize() == tree.getExternalNodeCount();
    }

    public Clade getRootClade() {
        return rootClade;
    }

    /**
     * recursively add all the clades in a tree
     */
    private Clade addClades(Tree tree, NodeRef node) {
        Clade clade;
        if (tree.isExternal(node)) {
            int index = node.getNumber();
            clade = getOrAddClade(index);
        } else {
            assert tree.getChildCount(node) == 2 : "requires a strictly bifurcating tree";

            Clade clade1 = addClades(tree, tree.getChild(node, 0));
            Clade clade2 = addClades(tree, tree.getChild(node, 1));
            clade = getOrAddClade(clade1, clade2);
        }

        clade.setCount(clade.getCount() + 1);

        return clade;
    }

    /**
     * see if a tip clade exists otherwise create it
     */
    private Clade getOrAddClade(int tipIndex) {
        BiClade clade = (BiClade)tipClades.get(tipIndex);
        if (clade == null) {
            clade = new BiClade(tipIndex);
            tipClades.put(tipIndex, clade);
        }
        return clade;
    }

    /**
     * see if a clade exists otherwise create it
     */
    private Clade getOrAddClade(Clade child1, Clade child2) {
        BiClade clade = (BiClade)cladeMap.get(BiClade.getKey((BiClade) child1, (BiClade) child2));
        if (clade == null) {
            clade = new BiClade(child1, child2);
            cladeMap.put(clade.getKey(), clade);
        } else {
            clade.addSubClades(child1, child2);
        }
        return clade;
    }

    public Clade getClade(Object key) {
        if (key instanceof Integer) {
            return tipClades.get(key);
        }
        return cladeMap.get(key);
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
            clade = new BiClade(index);

        } else {
            assert tree.getChildCount(node) == 2;

            Clade clade1 = collectAttributes(attributeNames, tree, tree.getChild(node, 0));
            Clade clade2 = collectAttributes(attributeNames, tree, tree.getChild(node, 1));
            clade = new BiClade(clade1, clade2);
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

    public void calculateCladeCredibilities(int totalTreesUsed) {
        for (Clade clade : cladeMap.values()) {
            assert clade.getCount() <= totalTreesUsed : "clade.getCount=(" + clade.getCount() +
                        ") should be <= totalTreesUsed = (" + totalTreesUsed + ")";

            clade.setCredibility(((double) clade.getCount()) / (double) totalTreesUsed);
        }
    }

    public double getLogCladeCredibility(Tree tree) {
        double[] logCladeCredibility = { 0.0 };
        logCladeCredibility(tree, tree.getRoot(), logCladeCredibility);
        return logCladeCredibility[0];
    }

    private Clade logCladeCredibility(Tree tree, NodeRef node, double[] logCladeCredibility) {

        Clade clade;

        if (tree.isExternal(node)) {
            int index = node.getNumber();
            clade = new BiClade(index);
        } else {

            assert tree.getChildCount(node) == 2;

            Clade clade1 = logCladeCredibility(tree, tree.getChild(node, 0), logCladeCredibility);
            Clade clade2 = logCladeCredibility(tree, tree.getChild(node, 1), logCladeCredibility);

            clade = new BiClade(clade1, clade2);

            logCladeCredibility[0] += Math.log(getCladeCredibility(clade));
        }

        return clade;
    }

    private double getCladeCredibility(Clade keyClade) {
        Clade clade = cladeMap.get(keyClade.getKey());
        assert clade != null;
        return clade.getCredibility();
    }


    public int getCladeCount() {
        return cladeMap.keySet().size();
    }

    //
    // Private stuff
    //
    TaxonList taxonList = null;
    private final Map<Object, Clade> tipClades = new HashMap<>();
   private final Map<Object, Clade> cladeMap = new HashMap<>();

    Clade rootClade;
}
