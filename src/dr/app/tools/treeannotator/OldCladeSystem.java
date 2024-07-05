package dr.app.tools.treeannotator;

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
class OldCladeSystem {
    //
    // Public stuff
    //

    private final TreeAnnotator treeAnnotator;

    /**
     *
     */
    public OldCladeSystem(TreeAnnotator treeAnnotator) {
        this.treeAnnotator = treeAnnotator;
    }

    /**
     *
     */
    public OldCladeSystem(TreeAnnotator treeAnnotator, Tree targetTree) {
        this.treeAnnotator = treeAnnotator;
        this.targetTree = targetTree;
        add(targetTree, true);
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
        BitSet rootBits = addClades(tree, tree.getRoot(), includeTips);
        rootClade = cladeMap.get(rootBits);
    }

    public Clade getRootClade() {
        return rootClade;
    }

    private BitSet addClades(Tree tree, NodeRef node, boolean includeTips) {

        BitSet bits = new BitSet();

        if (tree.isExternal(node)) {

//                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
            int index = node.getNumber();
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

    public void collectAttributes(Tree tree) {
        collectAttributes(tree, tree.getRoot());
    }

    private BitSet collectAttributes(Tree tree, NodeRef node) {

        BitSet bits = new BitSet();

        if (tree.isExternal(node)) {

//                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
//                if (index < 0) {
//                    throw new IllegalArgumentException("Taxon, " + tree.getNodeTaxon(node).getId() + ", not found in target tree");
//                }
            int index = node.getNumber();
            bits.set(index);

        } else {

            for (int i = 0; i < tree.getChildCount(node); i++) {

                NodeRef node1 = tree.getChild(node, i);

                bits.or(collectAttributes(tree, node1));
            }
        }

        collectAttributesForClade(bits, tree, node);

        return bits;
    }

    private void collectAttributesForClade(BitSet bits, Tree tree, NodeRef node) {
        Clade clade = cladeMap.get(bits);
        if (clade != null) {

            if (clade.attributeValues == null) {
                clade.attributeValues = new ArrayList<Object[]>();
            }

            int i = 0;
            Object[] values = new Object[treeAnnotator.attributeNames.size()];
            for (String attributeName : treeAnnotator.attributeNames) {
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
    Map<BitSet, Clade> cladeMap = new HashMap<>();

    Clade rootClade;

    Tree targetTree;
}
