package dr.app.treespace;

import dr.evolution.tree.*;
import dr.evolution.util.TaxonList;
import java.util.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class CladeSystem {
    /**
     */
    public CladeSystem() {
    }

    /**
     * adds all the clades in the tree
     */
    public void add(Tree tree, boolean includeTips) {
        if (taxonList == null) {
            taxonList = tree;
        }

        // Recurse over the tree and add all the clades (or increment their
        // frequency if already present). The root clade is added too (for
        // annotation purposes).
        addClades(tree, tree.getRoot(), includeTips);
    }
//
//        public Clade getClade(NodeRef node) {
//            return null;
//        }

    private BitSet addClades(Tree tree, NodeRef node, boolean includeTips) {

        BitSet bits = new BitSet();

        if (tree.isExternal(node)) {

            int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
            bits.set(index);

            if (includeTips) {
                addClade(bits);
            }

        } else {

            for (int i = 0; i < tree.getChildCount(node); i++) {

                NodeRef node1 = tree.getChild(node, i);

                bits.or(addClades(tree, node1, includeTips));
            }

            addClade(bits);
        }

        return bits;
    }

    private void addClade(BitSet bits) {
        Clade clade = cladeMap.get(bits);
        if (clade == null) {
            clade = new Clade(bits);
            cladeMap.put(bits, clade);
        }
        clade.setCount(clade.getCount() + 1);
    }

    public Map getCladeMap() {
        return cladeMap;
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

            int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
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

            int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
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

    public List<Clade> getClades() {
        return new ArrayList<Clade>(cladeMap.values());
    }


    class Clade {
        public Clade(BitSet bits) {
            this.bits = bits;
            count = 0;
            credibility = 0.0;
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
        BitSet bits;
    }

    //
    // Private stuff
    //
    TaxonList taxonList = null;
    Map<BitSet, Clade> cladeMap = new HashMap<BitSet, Clade>();

    Tree targetTree;
}
