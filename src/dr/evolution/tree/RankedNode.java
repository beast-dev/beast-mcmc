/*
 * RankedNode.java
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


import jebl.evolution.io.ImportException;
import jebl.evolution.io.NewickImporter;
import jebl.evolution.trees.RootedTree;

import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

/**
 * @author Alexei Drummond
 */
public class RankedNode {

    int rank;
    int n;
    RankedNode child1, child2;
    //Clade clade;

    BitSet cladeBits;


    public RankedNode(int i, int n) {
        rank = i - n;
        this.n = n;
        //clade = new TerminalClade(rank);
        child1 = null;
        child2 = null;

        cladeBits = new BitSet(n);
        cladeBits.set(i);
    }

    public RankedNode(int rank, RankedNode child1, RankedNode child2) {

        this.rank = rank;
        n = child1.n;


        if (!child1.cladeBits.intersects(child2.cladeBits)) {

            cladeBits = new BitSet();
            cladeBits.or(child1.cladeBits);
            cladeBits.or(child2.cladeBits);
            this.child1 = child1;
            this.child2 = child2;
        } else throw new IllegalArgumentException();
    }

    static BitSet inter = new BitSet();

    public boolean compatible(BitSet constraint) {

        if (cladeBits.intersects(constraint)) {
            inter.clear();
            inter.or(cladeBits);
            inter.and(constraint);

            //System.out.println(cladeBits + " and " + constraint + " : compatible=" + Math.min(cladeBits.cardinality(), constraint.cardinality()));

            return (inter.cardinality() == Math.min(cladeBits.cardinality(), constraint.cardinality()));
        }
        return true;
    }

    public boolean isExternal() {
        return (child1 == null) && (child2 == null);
    }

    public RankedNode[] getChildren() {
        if (isExternal()) return new RankedNode[]{};
        return new RankedNode[]{child1, child2};
    }

    public RootedTree getTree() {

        try {
            NewickImporter importer = new NewickImporter(new StringReader(toNewick()), true);
            return (RootedTree) importer.importNextTree();
        } catch (Exception e) {
            return null;
        }

    }

    public String toNewick() {
        return toNewick(this, null);
    }

    private String toNewick(RankedNode node, RankedNode parent) {

        if (node.isExternal()) {
            return ((char) (node.rank + node.n + 65)) + ":" + (parent.rank + 1);
        } else {

            if (node.child1.compare(node.child2) > 0) {
                RankedNode temp = node.child1;
                node.child1 = node.child2;
                node.child2 = temp;
            }

            return "(" + toNewick(node.child1, node) + ", " + toNewick(node.child2, node) + "):" + ((parent != null) ? parent.rank - node.rank : 0);
        }
    }

    private int compare(RankedNode node2) {
        if (isExternal()) {
            if (node2.isExternal()) {
                return rank - node2.rank;
            } else {
                return 1 - node2.cladeBits.cardinality();
            }
        } else {
            return cladeBits.cardinality() - node2.cladeBits.cardinality();
        }
    }

//    class Clade {
//
//        boolean[] in;
//
//        public Clade() {
//            in = new boolean[n];
//        }
//
//        public Clade(Clade c1, Clade c2) {
//            in = new boolean[n];
//
//            for (int i = 0; i < n; i++) {
//                if (c1.contains(i - n) && c2.contains(i - n)) {
//                    throw new IllegalArgumentException();
//                } else {
//                    in[i] = c1.contains(i - n) || c2.contains(i - n);
//                }
//            }
//        }
//
//        public void add(int node) {
//            in[node + n] = true;
//        }
//
//        public boolean contains(int node) {
//            return in[node + n];
//        }
//
//        public void set(Clade clade) {
//            System.arraycopy(clade.in, 0, in, 0, in.length);
//        }
//
//        public boolean equals(Clade c) {
//            for (int i = 0; i < in.length; i++) {
//                if (c.in[i] != in[i]) return false;
//            }
//            return true;
//        }
//
//        public boolean compatible(Clade clade) {
//
//            return outsize(clade) == 0 || insize(clade) == 0 || insize(clade) == clade.size();
//        }
//
//        public int size() {
//            int size = 0;
//            for (int i = 0; i < n; i++) {
//                size += in[i] ? 1 : 0;
//            }
//            return size;
//        }
//
//        private int insize(Clade clade) {
//            int insize = 0;
//            for (int i = 0; i < n; i++) {
//                if (contains(i - n)) {
//                    if (clade.contains(i - n)) {
//                        insize += 1;
//                    }
//                }
//            }
//            return insize;
//        }
//
//        private int outsize(Clade clade) {
//            int outsize = 0;
//            for (int i = 0; i < n; i++) {
//                if (contains(i - n)) {
//                    if (!clade.contains(i - n)) {
//                        outsize += 1;
//                    }
//                }
//            }
//            return outsize;
//        }
//
//
//        public boolean disjoint(Clade clade) {
//            for (int i = 0; i < n; i++) {
//                if (contains(i - n) && clade.contains(i - n)) {
//                    return false;
//                }
//            }
//            return true;
//        }
//
//        public String toString() {
//            StringBuilder builder = new StringBuilder();
//            for (boolean b : in) {
//                if (b) builder.append('1');
//                else builder.append('0');
//            }
//            return builder.toString();
//        }
//
//    }
//
//    class TerminalClade extends Clade {
//        int terminal;
//
//        public TerminalClade(int terminal) {
//            this.terminal = terminal;
//        }
//
//        public void add(int node) {
//            throw new IllegalArgumentException();
//        }
//
//        public boolean contains(int node) {
//            return terminal == node;
//        }
//
//        public void set(Clade clade) {
//            throw new IllegalArgumentException();
//        }
//
//        public String toString() {
//            StringBuilder builder = new StringBuilder();
//            for (int i = 0; i < n; i++) {
//                if (contains(i - n)) builder.append('1');
//                else builder.append('0');
//            }
//            return builder.toString();
//        }
//    }
    static long[] Rn;

    public static void main(String[] args) throws IOException, ImportException {

//        String[] constraintStrings = {
//                  "111110000000000",
//                  "111111110000000",
//                  "111111111000000",
//                  "000000000110000",
//                  "000000000001100",
//                  "000000000111100",
//                  "111111111111111"};

        String[] constraintStrings = {
                  "111110000",
                  "111111110"};

//        //String[] constraintStrings = {"00000000011", "00000001100"};

//        String[] constraintStrings = { "0000011", "1110000", "1111111"};

        int n = constraintStrings[0].length();

        Rn = new long[12];
        Rn[0] = 1;
        Rn[1] = 1;
        for (int i = 2; i < Rn.length; i++) {
            Rn[i] = Rn[i - 1] * ((i + 1) * i / 2);
            System.out.println((i + 1) + "\t" + Rn[i]);
        }

        RankedForest start = new RankedForest.Default(n, false);

        List<RankedNode> complete = new ArrayList<RankedNode>();

        List<BitSet> constraints = new ArrayList<BitSet>();

        for (String constraintString : constraintStrings) {

            constraints.add(start.getNodes().get(0).createClade(constraintString));
        }
        int[] count = new int[]{0, 0};

        long startTime = System.currentTimeMillis();

        processHistory(start, complete, constraints, count);

        long finishTime = System.currentTimeMillis();

        System.out.println("n = " + n);

        System.out.println("Constraints:");
        for (BitSet constraint : constraints) {
            System.out.println("  " + constraint);
        }

        int total = complete.size() + count[1];

        System.out.println(total + " histories found in " + count[0] + " calls");

        System.out.println("Took " + Math.round((finishTime - startTime) / 10.0) / 100.0 + " seconds");

        System.out.println("Max size of clear forest: " + max);
    }

    private BitSet createClade(String s) {
        BitSet bitSet = new BitSet();
        for (int i = 0; i < s.length(); i++) {
            if (s.charAt(i) == '1') {
                bitSet.set(i);
            }
        }
        return bitSet;
    }

    static int max = 1;

    private static void processHistory(RankedForest history, List<RankedNode> completeHistories, List<BitSet> constraints, int[] count) {

        count[0] += 1;

        List<RankedNode> nodes = history.getNodes();
        int k = nodes.size();
        //System.out.println("k = " + k);

        if (k == 1) {
            count[1] += 1;
            if (count[0] % 10000000 == 0) System.out.println(count[1]);
            //completeHistories.add(nodes.get(0));
        } else {

            if (history.isClear()) {
                int treeCount = history.getNodes().size();
                if (treeCount > max) max = treeCount;

                count[1] += Rn[treeCount - 1];
            } else {

                for (int i = 0; i < k; i++) {
                    for (int j = i + 1; j < k; j++) {

                        RankedNode parent = new RankedNode(history.rank() + 1, nodes.get(i), nodes.get(j));

                        boolean compatible = true;
                        if (constraints != null) {
                            for (BitSet constraint : constraints) {
                                if (!parent.compatible(constraint)) {
                                    //System.out.println(parent.cladeBits + " not compatible with " + constraint);
                                    compatible = false;
                                    break;
                                }
                            }
                        }
                        if (compatible) {

                            //System.out.println(parent.toNewick() + " is compatible!");

                            RankedForest newHis = new RankedForest.Parent(history, parent, constraints);

                            // check if order of constraints okay
                            if (newHis.compatibleRank(constraints)) {
                                processHistory(newHis, completeHistories, constraints, count);
                            }
                        }
                        //System.out.println(i + " " + j);
                    }
                }
            }
        }
    }
}