/*
 * BranchScoreMetric.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
package dr.evolution.tree;

import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import jebl.evolution.taxa.Taxon;
import jebl.evolution.graphs.Node;
import jebl.evolution.treemetrics.RootedTreeMetric;
import jebl.evolution.trees.RootedTree;

/**
 * @author Sebastian Hoehna
 * @version 1.0
 *
 */
public class BranchScoreMetric implements RootedTreeMetric {

    public BranchScoreMetric() {
        taxonMap = null;
    }

    public BranchScoreMetric(List<Taxon> taxa) {
        taxonMap = new HashMap<Taxon, Integer>();
        for (int i = 0; i < taxa.size(); i++) {
            taxonMap.put(taxa.get(i), i);
        }
    }

    public double getMetric(RootedTree tree1, RootedTree tree2) {

        if (!tree1.getTaxa().equals(tree2.getTaxa())) {
            throw new IllegalArgumentException("Trees contain different taxa");
        }

        Map<Taxon, Integer> tm = taxonMap;

        if (tm == null) {
            List<Taxon> taxa = new ArrayList<Taxon>(tree1.getTaxa());

            if (!tree2.getTaxa().equals(taxa))
                tm = new HashMap<Taxon, Integer>();
            for (int i = 0; i < taxa.size(); i++) {
                tm.put(taxa.get(i), i);
            }
        }

        List<Clade> clades1 = new ArrayList<Clade>();
        getClades(tm, tree1, tree1.getRootNode(), clades1, null);

        List<Clade> clades2 = new ArrayList<Clade>();
        getClades(tm, tree2, tree2.getRootNode(), clades2, null);

        return getDistance(clades1, clades2);
    }

    protected void getClades(Map<Taxon, Integer> taxonMap, RootedTree tree,
                           Node node, List<Clade> clades, BitSet bits) {

        BitSet bits2 = new BitSet();

        if (tree.isExternal(node)) {

            int index = taxonMap.get(tree.getTaxon(node));
            bits2.set(index);

        } else {

            for (Node child : tree.getChildren(node)) {
                getClades(taxonMap, tree, child, clades, bits2);
            }
            clades.add(new Clade(bits2, tree.getHeight(node)));
        }

        if (bits != null) {
            bits.or(bits2);
        }
    }

    protected double getDistance(List<Clade> clades1, List<Clade> clades2) {

        Collections.sort(clades1);
        Collections.sort(clades2);
        double distance = 0.0;
        int indexClade2 = 0;
        Clade clade2 = null;
        Clade parent1, parent2 = null;
        double height1, height2;

        for (Clade clade1 : clades1) {

            parent1 = findParent(clade1, clades1);
            height1 = parent1.getHeight() - clade1.getHeight();

            if (indexClade2 < clades2.size()) {
                clade2 = clades2.get(indexClade2);
                parent2 = findParent(clade2, clades2);
            }
            while (clade1.compareTo(clade2) > 0 && indexClade2 < clades2.size()) {
                height2 = parent2.getHeight() - clade2.getHeight();
                distance += height2 * height2;
                indexClade2++;
                if (indexClade2 < clades2.size()) {
                    clade2 = clades2.get(indexClade2);
                    parent2 = findParent(clade2, clades2);
                }
            }
            if (clade1.compareTo(clade2) == 0) {
                height2 = parent2.getHeight() - clade2.getHeight();
                distance += (height1 - height2) * (height1 - height2);
                indexClade2++;
            } else {
                distance += height1 * height1;
            }
        }

        return Math.sqrt(distance);
    }

    private Clade findParent(Clade clade1, List<Clade> clades) {
        Clade parent = null;
        for (Clade clade2 : clades) {
            if (isParent(clade2, clade1)) {
                if (parent == null || parent.getSize() > clade2.getSize())
                    parent = clade2;
            }
        }

        if (parent == null){
            //the case that this clade is the whole tree
            return clade1;
        }

        return parent;
    }

    private boolean isParent(Clade parent, Clade child) {
        if (parent.getSize() <= child.getSize()) {
            return false;
        }

        tmpBits.clear();
        tmpBits.or(parent.getBits());
        tmpBits.xor(child.getBits());

        return tmpBits.cardinality() < parent.getSize();
    }

    BitSet tmpBits = new BitSet();

    private final Map<Taxon, Integer> taxonMap;

    public static void main(String[] args) {

        try {

            NewickImporter importer = new NewickImporter("((('C':0.03365591238,'A':0.7225157402):0.306488578,'B':0.4572411443):0.4673149632,('D':0.7966438427,'E':0.8063645191):0.7478901469)");
            Tree treeOne = importer.importNextTree();
            System.out.println("tree 1: " + treeOne);

            importer = new NewickImporter("(('A':0.2333369483,'B':0.3468381313):0.5562255983,('C':0.8732210915,('D':0.9124725792,'E':0.1983703848):0.5252404297):0.2000638912)");
            Tree treeTwo = importer.importNextTree();
            System.out.println("tree 2: " + treeTwo + "\n");

            double metric = (new BranchScoreMetric().getMetric(TreeUtils.asJeblTree(treeOne), TreeUtils.asJeblTree(treeTwo)));

            System.out.println("rooted branch score metric = " + metric);

        } catch(Importer.ImportException ie) {
            System.err.println(ie);
        } catch(IOException ioe) {
            System.err.println(ioe);
        }

    }

}
