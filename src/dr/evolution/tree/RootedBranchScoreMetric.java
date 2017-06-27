/*
 * RootedBranchScoreMetric.java
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

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import jebl.evolution.graphs.Node;
import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.RootedTree;

import java.io.IOException;
import java.util.*;

/**
 * @author Guy Baele
 *
 */
public class RootedBranchScoreMetric extends BranchScoreMetric {

    public RootedBranchScoreMetric() {
        taxonMap = null;
    }

    public RootedBranchScoreMetric(List<Taxon> taxa) {
        taxonMap = new HashMap<Taxon, Integer>();
        for (int i = 0; i < taxa.size(); i++) {
            //System.err.println(i + " -> " + taxa.get(i));
            taxonMap.put(taxa.get(i), i);
        }
    }

    @Override
    public double getMetric(RootedTree tree1, RootedTree tree2) {

        if (!tree1.getTaxa().equals(tree2.getTaxa())) {
            throw new IllegalArgumentException("Trees contain different taxa");
        }

        Map<Taxon, Integer> tm = taxonMap;

        if (tm == null) {
            List<Taxon> taxa = new ArrayList<Taxon>(tree1.getTaxa());

            if (!tree2.getTaxa().equals(taxa)) {
                tm = new HashMap<Taxon, Integer>();
                for (int i = 0; i < taxa.size(); i++) {
                    //System.err.println(i + " -> " + taxa.get(i));
                    tm.put(taxa.get(i), i);
                }
            }
        }

        List<Clade> clades1 = new ArrayList<Clade>();
        getClades(tm, tree1, tree1.getRootNode(), clades1, null);

        List<Clade> clades2 = new ArrayList<Clade>();
        getClades(tm, tree2, tree2.getRootNode(), clades2, null);

        return Math.sqrt(Math.pow(getDistance(clades1, clades2),2) + getExternalDistance(tree1, tree2));
    }

    private double getExternalDistance(RootedTree tree1, RootedTree tree2) {

        double distance = 0.0;

        List<Taxon> taxaOne = new ArrayList<Taxon>(tree1.getTaxa());
        List<Taxon> taxaTwo = new ArrayList<Taxon>(tree2.getTaxa());

        for (Taxon tOne : taxaOne) {
            for (Taxon tTwo : taxaTwo) {
                if (tOne.equals(tTwo)) {
                    Node parentOne = tree1.getParent(tree1.getNode(tOne));
                    Node parentTwo = tree2.getParent(tree2.getNode(tTwo));
                    distance += Math.pow((tree1.getHeight(parentOne)-tree1.getHeight(tree1.getNode(tOne)))-(tree2.getHeight(parentTwo)-tree2.getHeight(tree2.getNode(tTwo))),2);
                }
            }
        }

        return distance;

    }

    /*private void getClades(Map<Taxon, Integer> taxonMap, RootedTree tree,
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

    private double getDistance(List<Clade> clades1, List<Clade> clades2) {

        Collections.sort(clades1);
        Collections.sort(clades2);

        double distance = 0.0;
        ArrayList<Clade> checkClades = new ArrayList<Clade>();

        for (Clade clade1 : clades1) {
            boolean found = false;
            for (Clade clade2 : clades2) {
                if (clade1.compareTo(clade2) == 0) {
                    Clade parent1 = findParent(clade1, clades1);
                    Clade parent2 = findParent(clade2, clades2);
                    distance += Math.pow((parent1.getHeight()-clade1.getHeight()) - (parent2.getHeight()-clade2.getHeight()),2);
                    System.err.println(clade1 + " vs. " + clade2 + " = " + Math.pow((parent1.getHeight()-clade1.getHeight()) - (parent2.getHeight()-clade2.getHeight()),2));
                    found = true;
                    checkClades.add(clade2);
                }
            }
            if (!found) {
                Clade parent1 = findParent(clade1, clades1);
                distance += Math.pow(parent1.getHeight()-clade1.getHeight(),2);
                System.err.println(clade1 + " = " + Math.pow(parent1.getHeight()-clade1.getHeight(),2));
            }
        }

        for (Clade clade2 : clades2) {
            if (!checkClades.contains(clade2)) {
                boolean found = false;
                for (Clade clade1 : clades1) {
                    if (clade1.compareTo(clade2) == 0) {
                        Clade parent1 = findParent(clade1, clades1);
                        Clade parent2 = findParent(clade2, clades2);
                        distance += Math.pow((parent1.getHeight()-clade1.getHeight()) - (parent2.getHeight()-clade2.getHeight()),2);
                        System.err.println(clade1 + " vs. " + clade2 + " = " + Math.pow((parent1.getHeight()-clade1.getHeight()) - (parent2.getHeight()-clade2.getHeight()),2));
                        found = true;
                        checkClades.add(clade2);
                    }
                }
                if (!found) {
                    Clade parent2 = findParent(clade2, clades2);
                    distance += Math.pow(parent2.getHeight()-clade2.getHeight(),2);
                    System.err.println(clade2 + " = " + Math.pow(parent2.getHeight()-clade2.getHeight(),2));
                }
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
    */

    private final Map<Taxon, Integer> taxonMap;

    public static void main(String[] args) {

        try {

            NewickImporter importer = new NewickImporter("((('C':0.03365591238,'A':0.7225157402):0.306488578,'B':0.4572411443):0.4673149632,('D':0.7966438427,'E':0.8063645191):0.7478901469)");
            Tree treeOne = importer.importNextTree();
            System.out.println("tree 1: " + treeOne);

            importer = new NewickImporter("(('A':0.2333369483,'B':0.3468381313):0.5562255983,('C':0.8732210915,('D':0.9124725792,'E':0.1983703848):0.5252404297):0.2000638912)");
            Tree treeTwo = importer.importNextTree();
            System.out.println("tree 2: " + treeTwo + "\n");

            double metric = (new RootedBranchScoreMetric().getMetric(TreeUtils.asJeblTree(treeOne), TreeUtils.asJeblTree(treeTwo)));

            System.out.println("rooted branch score metric = " + metric);

        } catch(Importer.ImportException ie) {
            System.err.println(ie);
        } catch(IOException ioe) {
            System.err.println(ioe);
        }

    }

}
