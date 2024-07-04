package dr.app.tools.treeannotator;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.util.Attributable;
import dr.util.Pair;

import java.util.*;

/**
 * @author Andrew Rambaut
 * @version $
 */
class CladeSystem2 {
    //
    // Public stuff
    //

    private final TreeAnnotator treeAnnotator;

    /**
     *
     */
    public CladeSystem2(TreeAnnotator treeAnnotator) {
        this.treeAnnotator = treeAnnotator;
    }

    /**
     *
     */
    public CladeSystem2(TreeAnnotator treeAnnotator, Tree targetTree) {
        this.treeAnnotator = treeAnnotator;
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

            assert tree.getChildCount(node) == 2;

            Clade clade1 = collectAttributes(tree, tree.getChild(node, 0));
            Clade clade2 = collectAttributes(tree, tree.getChild(node, 1));
            clade = new Clade(clade1, clade2);

//            for (int i = 0; i < tree.getChildCount(node); i++) {
//
//                NodeRef node1 = tree.getChild(node, i);
//                BitSet subClade = addClades(tree, node1, includeTips);
//                bits.or(subClade);
//                subClades.add(subClade);
//            }
//
//            Clade clade = addClade(bits);
//
//            if (subClades.size() != 2) {
//                throw new IllegalArgumentException("TreeAnnotator requires strictly bifurcating trees");
//            }
//            clade.addSubclades(subClades.get(0), subClades.get(1));
        }

        return addClade(clade);
    }

    private Clade addClade(Clade clade) {
        Clade c = cladeMap.get(clade.hashCode());
        if (c == null) {
            cladeMap.put(clade.hashCode(), clade);
            c = clade;
        }
        c.setCount(c.getCount() + 1);
        return c;
    }

    public void collectAttributes(Tree tree) {
        collectAttributes(tree, tree.getRoot());
    }

    private Clade collectAttributes(Tree tree, NodeRef node) {

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

            Clade clade1 = collectAttributes(tree, tree.getChild(node, 0));
            Clade clade2 = collectAttributes(tree, tree.getChild(node, 1));
            clade = new Clade(clade1, clade2);
        }

        collectAttributesForClade(clade, tree, node);

        return clade;
    }

    private void collectAttributesForClade(Clade clade, Tree tree, NodeRef node) {
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

//    public Map<BitSet, Clade> getCladeMap() {
//        return cladeMap;
//    }
//
//    public Clade getClade(BitSet bitSet) {
//        return cladeMap.get(bitSet);
//    }

//    public void calculateCladeCredibilities(int totalTreesUsed) {
//        for (Clade clade : cladeMap.values()) {
//
//            if (clade.getCount() > totalTreesUsed) {
//
//                throw new AssertionError("clade.getCount=(" + clade.getCount() +
//                        ") should be <= totalTreesUsed = (" + totalTreesUsed + ")");
//            }
//
//            clade.setCredibility(((double) clade.getCount()) / (double) totalTreesUsed);
//        }
//    }

//    public double getLogCladeCredibility(Tree tree, NodeRef node, BitSet bits) {
//
//        double logCladeCredibility = 0.0;
//
//        if (tree.isExternal(node)) {
//
////                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
//            int index = node.getNumber();
//            bits.set(index);
//        } else {
//
//            BitSet bits2 = new BitSet();
//            for (int i = 0; i < tree.getChildCount(node); i++) {
//
//                NodeRef node1 = tree.getChild(node, i);
//
//                logCladeCredibility += getLogCladeCredibility(tree, node1, bits2);
//            }
//
//            logCladeCredibility += Math.log(getCladeCredibility(bits2));
//
//            if (bits != null) {
//                bits.or(bits2);
//            }
//        }
//
//        return logCladeCredibility;
//    }

    private double getCladeCredibility(int hash) {
        Clade clade = cladeMap.get(hash);
        if (clade == null) {
            return 0.0;
        }
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

    class Clade extends Attributable {
        public Clade(int index) {
            this.index = index;
            count = 0;
            credibility = 1.0;
            size = 1;
            bits = null;
        }

        public Clade(Clade subClade1, Clade subClade2) {
            count = 0;
            credibility = 0.0;
            size = subClade1.size + subClade2.size;

            index = Math.min(subClade1.minIndex(), subClade2.minIndex());
            int maxIndex = Math.max(subClade1.maxIndex(), subClade2.maxIndex());
            bits = new byte[maxIndex - index + 1];
            if (subClade1.minIndex() < subClade2.minIndex()) {
                if (subClade1.bits == null) {
                    bits[subClade1.index] = 1;
                } else {
                    System.arraycopy(subClade1.bits, 0, bits, 0, subClade1.bits.length);
                    if (subClade2.bits == null) {
                        bits[subClade2.index] = 1;
                    } else {
                        for (int i = 0; i < subClade2.bits.length; i++) {
                            bits[i + subClade2.minIndex()] |= subClade2.bits[0];
                        }
                    }
                }
            } else {
                if (subClade2.bits == null) {
                    bits[subClade2.index] = 1;
                } else {
                    System.arraycopy(subClade2.bits, 0, bits, 0, subClade2.bits.length);
                    if (subClade1.bits == null) {
                        bits[subClade1.index] = 1;
                    } else {
                        for (int i = 0; i < subClade1.bits.length; i++) {
                            bits[i + subClade1.minIndex()] |= subClade1.bits[0];
                        }
                    }
                }
            }

        }

        private int minIndex() {
            return index;
        }

        private int maxIndex() {
            if (index == -1) {
                return index + bits.length - 1;
            } else {
                return index;
            }
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
            return (bits != null ? Arrays.hashCode(bits) : 0);
        }

        public String toString() {
            return "clade " + hashCode();
        }

        int count;
        double credibility;
        final int size;
        final byte[] bits;
        final int index;
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
    Map<Integer, Clade> cladeMap = new HashMap<>();

    Clade rootClade;

    Tree targetTree;
}
