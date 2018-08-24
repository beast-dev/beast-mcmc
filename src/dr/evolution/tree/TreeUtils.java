package dr.evolution.tree;

import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import jebl.evolution.graphs.Node;
import jebl.evolution.trees.SimpleRootedTree;

import java.text.NumberFormat;
import java.util.*;

/**
 * Static utility functions for trees.
 */
public class TreeUtils {

    /**
     * Count number of leaves in subtree whose root is node.
     *
     * @param tree the tree
     * @param node the node to get leaf count below
     * @return the number of leaves under this node.
     */
    public static int getLeafCount(Tree tree, NodeRef node) {

        int childCount = tree.getChildCount(node);
        if (childCount == 0) return 1;

        int leafCount = 0;
        for (int i = 0; i < childCount; i++) {
            leafCount += getLeafCount(tree, tree.getChild(node, i));
        }
        return leafCount;
    }

    public static double getTreeLength(Tree tree, NodeRef node) {

        int childCount = tree.getChildCount(node);
        if (childCount == 0) return tree.getBranchLength(node);

        double length = 0;
        for (int i = 0; i < childCount; i++) {
            length += getTreeLength(tree, tree.getChild(node, i));
        }
        if (node != tree.getRoot())
            length += tree.getBranchLength(node);
        return length;

    }

    public static double getSubTreeLength(Tree tree, Set<Taxon> taxa) {

        double[] treeLength = new double[] { 0.0 };

        getSubTreeLength(tree, tree.getRoot(), taxa, treeLength);

        return treeLength[0];

    }

    private static int getSubTreeLength(Tree tree, NodeRef node, Set<Taxon> taxa, double[] treeLength) {
        int childCount = tree.getChildCount(node);
        if (childCount == 0) {
            if (taxa.contains(tree.getNodeTaxon(node))) {
                treeLength[0] += tree.getBranchLength(node);
                return 1;
            } else {
                return 0;
            }
        }

        int count = 0;
        for (int i = 0; i < childCount; i++) {
            count += getSubTreeLength(tree, tree.getChild(node, i), taxa, treeLength);
        }

        if (count > 0 && count < taxa.size()) {
            if (node != tree.getRoot()) {
                treeLength[0] += tree.getBranchLength(node);
            }
            return count;
        }

        return 0;
    }

    public static double getMinNodeHeight(Tree tree, NodeRef node) {

        int childCount = tree.getChildCount(node);
        if (childCount == 0) return tree.getNodeHeight(node);

        double minNodeHeight = Double.MAX_VALUE;
        for (int i = 0; i < childCount; i++) {
            double height = getMinNodeHeight(tree, tree.getChild(node, i));
            if (height < minNodeHeight) {
                minNodeHeight = height;
            }
        }
        return minNodeHeight;
    }

    /**
     * @param tree the tree to test fo ultrametricity
     * @return true only if all tips have height 0.0
     */
    public static boolean isUltrametric(Tree tree) {
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            if (tree.getNodeHeight(tree.getExternalNode(i)) != 0.0)
                return false;
        }
        return true;
    }

    /**
     * @param tree the tree to test if binary
     * @return true only if internal nodes have 2 children
     */
    public static boolean isBinary(Tree tree) {
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            if (tree.getChildCount(tree.getInternalNode(i)) > 2)
                return false;
        }
        return true;
    }

    /**
     * @param tree the tree to retrieve leaf set of
     * @return a set of strings which are the taxa of the tree.
     */
    public static Set<String> getLeafSet(Tree tree) {

        HashSet<String> leafSet = new HashSet<String>();
        int m = tree.getTaxonCount();

        for (int i = 0; i < m; i++) {

            Taxon taxon = tree.getTaxon(i);
            leafSet.add(taxon.getId());
        }

        return leafSet;
    }

    /**
     * @param tree the tree
     * @param taxa the taxa
     * @return Set of taxon names (id's) associated with the taxa in taxa.
     * @throws MissingTaxonException
     *          if a taxon in taxa is not contained in the tree
     */
    public static Set<String> getLeavesForTaxa(Tree tree, TaxonList taxa) throws MissingTaxonException {

        HashSet<String> leafNodes = new HashSet<String>();
        int m = taxa.getTaxonCount();
        int n = tree.getExternalNodeCount();

        for (int i = 0; i < m; i++) {

            Taxon taxon = taxa.getTaxon(i);
            boolean found = false;
            for (int j = 0; j < n; j++) {

                NodeRef node = tree.getExternalNode(j);
                if (tree.getNodeTaxon(node).getId().equals(taxon.getId())) {

                    found = true;
                    break;
                }
            }

            if (!found) {
                throw new MissingTaxonException(taxon);
            }

            leafNodes.add(taxon.getId());
        }

        return leafNodes;
    }

    /**
     * @param tree the tree
     * @param node the node to get names of leaves below
     * @return a set of taxa names (as strings) of the leaf nodes descended from the given node.
     */
    public static Set<String> getDescendantLeaves(Tree tree, NodeRef node) {

        HashSet<String> set = new HashSet<String>();
        getDescendantLeaves(tree, node, set);
        return set;
    }

    /**
     * @param tree the tree
     * @param node the node to get name of leaves below
     * @param set  will be populated with taxa names (as strings) of the leaf nodes descended from the given node.
     */
    private static void getDescendantLeaves(Tree tree, NodeRef node, Set<String> set) {

        if (tree.isExternal(node)) {
            set.add(tree.getTaxonId(node.getNumber()));
        } else {

            for (int i = 0; i < tree.getChildCount(node); i++) {

                NodeRef node1 = tree.getChild(node, i);

                getDescendantLeaves(tree, node1, set);
            }
        }
    }


    /**
     * @param tree the tree
     * @param node the node to get external nodes below
     * @return a set of noderefs of the leaf nodes descended from the given node.
     */
    public static Set<NodeRef> getExternalNodes(Tree tree, NodeRef node) {

        HashSet<NodeRef> set = new HashSet<NodeRef>();
        getExternalNodes(tree, node, set);
        return set;
    }

    /**
     * @param tree the tree
     * @param node the node to get external nodes below
     * @param set  is populated with noderefs of the leaf nodes descended from the given node.
     */
    private static void getExternalNodes(Tree tree, NodeRef node, Set<NodeRef> set) {

        if (tree.isExternal(node)) {
            set.add(node);
        } else {

            for (int i = 0; i < tree.getChildCount(node); i++) {

                NodeRef node1 = tree.getChild(node, i);

                getExternalNodes(tree, node1, set);
            }
        }
    }

    /**
     * Gets the most recent common ancestor (MRCA) node of a set of leaf nodes.
     *
     * @param tree      the Tree
     * @param leafNodes a set of names
     * @return the NodeRef of the MRCA
     */
    public static NodeRef getCommonAncestorNode(Tree tree, Set<String> leafNodes) {

        int cardinality = leafNodes.size();

        if (cardinality == 0) {
            throw new IllegalArgumentException("No leaf nodes selected");
        }

        NodeRef[] mrca = {null};
        getCommonAncestorNode(tree, tree.getRoot(), leafNodes, cardinality, mrca);

        return mrca[0];
    }

    /*
     * Private recursive function used by getCommonAncestorNode.
     */
    private static int getCommonAncestorNode(Tree tree, NodeRef node,
                                             Set<String> leafNodes, int cardinality,
                                             NodeRef[] mrca) {

        if (tree.isExternal(node)) {

            if (leafNodes.contains(tree.getTaxonId(node.getNumber()))) {
                if (cardinality == 1) {
                    mrca[0] = node;
                }
                return 1;
            } else {
                return 0;
            }
        }

        int matches = 0;

        for (int i = 0; i < tree.getChildCount(node); i++) {

            NodeRef node1 = tree.getChild(node, i);

            matches += getCommonAncestorNode(tree, node1, leafNodes, cardinality, mrca);

            if (mrca[0] != null) {
                break;
            }
        }

        if (mrca[0] == null) {
            // If we haven't already found the MRCA, test this node
            if (matches == cardinality) {
                mrca[0] = node;
            }
        }

        return matches;
    }

    /**
     * @param tree the tree
     * @param taxa the taxa
     * @return A bitset with the node numbers set.
     * @throws MissingTaxonException
     *          if a taxon in taxa is not contained in the tree
     */
    public static BitSet getTipsBitSetForTaxa(Tree tree, TaxonList taxa) throws MissingTaxonException {

        BitSet tips = new BitSet();
        for (int n: getTipsForTaxa(tree, taxa)) {
            tips.set(n);
        }
        return tips;
    }

    /**
     * @param tree the tree
     * @param taxa the taxa
     * @return A HashSet of node numbers.
     * @throws MissingTaxonException
     *          if a taxon in taxa is not contained in the tree
     */
    public static Set<Integer> getTipsForTaxa(Tree tree, TaxonList taxa) throws MissingTaxonException {

        Set<Integer> tips = new LinkedHashSet<Integer>();

        for (int i = 0; i < taxa.getTaxonCount(); i++) {

            Taxon taxon = taxa.getTaxon(i);
            boolean found = false;
            for (int j = 0; j < tree.getExternalNodeCount(); j++) {

                NodeRef node = tree.getExternalNode(j);
                if (tree.getNodeTaxon(node).getId().equals(taxon.getId())) {
                    tips.add(node.getNumber());
                    found = true;
                    break;
                }
            }

            if (!found) {
                throw new MissingTaxonException(taxon);
            }
        }

        return tips;
    }

    public static boolean isMonophyletic(Tree tree, Set<String> leafNodes) {
        return isMonophyletic(tree, leafNodes, Collections.<String>emptySet());
    }

    /**
     * Performs a monophyly test on a set of leaf nodes. The nodes are monophyletic
     * if there is a node in the tree which subtends all the taxa in the set (and
     * only those taxa).
     *
     * @param tree      a tree object to perform test on
     * @param leafNodes a set of leaf node ids
     * @param ignore    a set of ids to ignore in monophyly assessment
     * @return boolean is monophyletic?
     */
    public static boolean isMonophyletic(Tree tree, Set<String> leafNodes, Set<String> ignore) {

        int cardinality = leafNodes.size();

        if (cardinality == 1) {
            // A single selected leaf is always monophyletic
            return true;
        }

        if (cardinality == tree.getExternalNodeCount()) {
            // All leaf nodes are selected
            return true;
        }

        if (cardinality == 0) {
            throw new IllegalArgumentException("No leaf nodes selected");
        }

        int[] matchCount = {0};
        int[] leafCount = {0};
        boolean[] isMono = {false};
        isMonophyletic(tree, tree.getRoot(), leafNodes, ignore, cardinality, matchCount, leafCount, isMono);

        return isMono[0];
    }

    /**
     * Private recursive function used by isMonophyletic.
     *
     * @param tree        a tree object to perform test on
     * @param node        the node that is currently being assessed in recursive procedure
     * @param leafNodes   a set of leaf node ids
     * @param ignore      a set of ids to ignore in monophyly assessment
     * @param cardinality the size of leafNodes set
     * @return boolean is monophyletic?
     */
    private static boolean isMonophyletic(Tree tree, NodeRef node,
                                          Set<String> leafNodes, Set<String> ignore,
                                          int cardinality,
                                          int[] matchCount, int[] leafCount,
                                          boolean[] isMono) {

        if (tree.isExternal(node)) {

            String id = tree.getNodeTaxon(node).getId();
            if (leafNodes.contains(id)) {
                matchCount[0] = 1;
            } else {
                matchCount[0] = 0;
            }
            if (!ignore.contains(id)) {
                leafCount[0] = 1;
            } else {
                leafCount[0] = 0;
            }
            return false;
        }

        int mc = 0;
        int lc = 0;

        for (int i = 0; i < tree.getChildCount(node); i++) {

            NodeRef node1 = tree.getChild(node, i);

            boolean done = isMonophyletic(tree, node1, leafNodes, ignore, cardinality, matchCount, leafCount, isMono);
            mc += matchCount[0];
            lc += leafCount[0];

            if (done) {
                return true;
            }
        }

        matchCount[0] = mc;
        leafCount[0] = lc;

        // If we haven't already found the MRCA, test this node
        if (mc == lc && lc == cardinality) {
            isMono[0] = true;
            return true;
        }

        return false;
    }

    public static NodeRef getCommonAncestor(Tree tree, NodeRef n1, NodeRef n2) {
        while( n1 != n2 ) {
            if( tree.getNodeHeight(n1) < tree.getNodeHeight(n2) ) {
                n1 = tree.getParent(n1);
            } else {
                n2 = tree.getParent(n2);
            }
        }
        return n1;
    }

    public static NodeRef getCommonAncestorSafely(Tree tree, NodeRef n1, NodeRef n2) {
        while( n1 != n2 ) {

            if (tree.isRoot(n1)) return n1;

            if (tree.isRoot(n2)) return n2;

            if( tree.getNodeHeight(n1) < tree.getNodeHeight(n2) ) {
                n1 = tree.getParent(n1);
            } else {
                n2 = tree.getParent(n2);
            }
        }
        return n1;
    }

    // A lightweight version for finding the most recent common ancestor of a group of taxa.
    // return the node-ref of the MRCA.

    // would be nice to use nodeRef's, but they are not preserved :(
    public static NodeRef getCommonAncestor(Tree tree, int[] nodes) {
        NodeRef cur = tree.getNode(nodes[0]);

        for(int k = 1; k < nodes.length; ++k) {
            cur = getCommonAncestor(tree, cur, tree.getNode(nodes[k]));
        }
        return cur;
    }


    /**
     * @param tree
     * @param range
     * @return the size of the largest clade with tips in the given range of times.
     */
    public static int largestClade(Tree tree, double range) {

        return largestClade(tree, tree.getRoot(), range, new double[]{0.0, 0.0});

    }

    /**
     * @return the size of the largest clade with tips in the given range of times.
     */
    private static int largestClade(Tree tree, NodeRef node, double range, double[] currentBounds) {

        if (tree.isExternal(node)) {
            currentBounds[0] = tree.getNodeHeight(node);
            currentBounds[1] = tree.getNodeHeight(node);
            return 1;
        } else {
            // get the bounds and max clade size of the left clade
            int cladeSize1 = largestClade(tree, tree.getChild(node, 0), range, currentBounds);
            double min = currentBounds[0];
            double max = currentBounds[1];

            // get the bounds and max clade size of the right clade
            int cladeSize2 = largestClade(tree, tree.getChild(node, 1), range, currentBounds);
            min = Math.min(min, currentBounds[0]);
            max = Math.max(max, currentBounds[1]);

            // update the joint bounds
            currentBounds[0] = min;
            currentBounds[1] = max;

            // if the joint clade is valid return the joint size
            if (max - min < range) {
                return cladeSize1 + cladeSize2;
            }
            // if the joint clade is not valid return the max of the two
            return Math.max(cladeSize1, cladeSize2);
        }
    }

    /**
     * Calculates the minimum number of steps for the parsimony reconstruction of a
     * binary character defined by leafStates.
     *
     * @param tree       a tree object to perform test on
     * @param leafStates a set of booleans, one for each leaf node
     * @return number of parsimony steps
     */
    public static int getParsimonySteps(Tree tree, Set leafStates) {

        int[] score = {0};
        getParsimonySteps(tree, tree.getRoot(), leafStates, score);
        return score[0];
    }

    private static int getParsimonySteps(Tree tree, NodeRef node, Set leafStates, int[] score) {

        if (tree.isExternal(node)) {
            return (leafStates.contains(tree.getTaxonId(node.getNumber())) ? 1 : 2);

        } else {

            int uState = getParsimonySteps(tree, tree.getChild(node, 0), leafStates, score);
            int iState = uState;

            for (int i = 1; i < tree.getChildCount(node); i++) {

                int state = getParsimonySteps(tree, tree.getChild(node, i), leafStates, score);
                uState = state | uState;

                iState = state & iState;

            }

            if (iState == 0) {
                score[0] += 1;
            }

            return uState;
        }

    }

    /**
     * Calculates the parsimony reconstruction of a binary character defined
     * by leafStates at a given node.
     *
     * @param tree       a tree object to perform test on
     * @param node       a NodeRef object from tree
     * @param leafStates a set of booleans, one for each leaf node
     * @return number of parsimony steps
     */
    public static double getParsimonyState(Tree tree, NodeRef node, Set leafStates) {

        int state = getParsimonyStateAtNode(tree, node, leafStates);
        switch (state) {
            case 1:
                return 0.0;
            case 2:
                return 1.0;
            default:
                return 0.5;
        }
    }

    private static int getParsimonyStateAtNode(Tree tree, NodeRef node, Set leafStates) {

        if (tree.isExternal(node)) {
            return (leafStates.contains(tree.getTaxonId(node.getNumber())) ? 1 : 2);

        } else {

            int uState = getParsimonyStateAtNode(tree, tree.getChild(node, 0), leafStates);
            int iState = uState;

            for (int i = 1; i < tree.getChildCount(node); i++) {

                int state = getParsimonyStateAtNode(tree, tree.getChild(node, i), leafStates);
                uState = state | uState;

                iState = state & iState;

            }

            return uState;
        }

    }

    /**
     * determine preorder successor of this node
     *
     * @return next node
     */
    public static NodeRef preorderSuccessor(Tree tree, NodeRef node) {

        NodeRef next = null;

        if (tree.isExternal(node)) {
            NodeRef cn = node, ln = null; // Current and last node

            // Go up
            do {
                if (tree.isRoot(cn)) {
                    next = cn;
                    break;
                }
                ln = cn;
                cn = tree.getParent(cn);
            }
            while (tree.getChild(cn, tree.getChildCount(cn) - 1) == ln);

            // Determine next node
            if (next == null) {
                // Go down one node
                for (int i = 0; i < tree.getChildCount(cn) - 1; i++) {

                    if (tree.getChild(cn, i) == ln) {
                        next = tree.getChild(cn, i + 1);
                        break;
                    }
                }
            }
        } else {
            next = tree.getChild(node, 0);
        }

        return next;
    }

    /**
     * determine a postorder traversal list of nodes in a tree
     *
     */
    public static void postOrderTraversalList(Tree tree, int[] postOrderList) {

        final int nodeCount = tree.getNodeCount();
        if (postOrderList.length != nodeCount) {
            throw new IllegalArgumentException("Illegal list length");
        }

        int idx = nodeCount - 1;
        int cidx = nodeCount - 1;

        postOrderList[idx] = tree.getRoot().getNumber();

        while (cidx > 0) {
            NodeRef cNode = tree.getNode(postOrderList[idx]);
            for(int i = 0; i < tree.getChildCount(cNode); ++i) {
                cidx -= 1;
                postOrderList[cidx] = tree.getChild(cNode, i).getNumber();
            }
            idx -= 1;
        }
    }

    // populate  postOrderList with node numbers in posorder: parent before children.
    //  postOrderList is pre-allocated with the right size
    public static void preOrderTraversalList(Tree tree, int[] postOrderList) {

        final int nodeCount = tree.getNodeCount();
        if (postOrderList.length != nodeCount) {
            throw new IllegalArgumentException("Illegal list length");
        }
        postOrderList[0] = tree.getRoot().getNumber();
        preOrderTraversalList(tree, 0, postOrderList);
    }

    static int preOrderTraversalList(Tree tree, int idx, int[] postOrderList) {
        final NodeRef node = tree.getNode(postOrderList[idx]);
        for(int i = 0; i < tree.getChildCount(node); ++i) {
            final NodeRef child = tree.getChild(node, i);
            idx += 1;
            postOrderList[idx] = child.getNumber();
            if( ! tree.isExternal(child) ) {
              idx = preOrderTraversalList(tree, idx, postOrderList);
            }
        }
        return idx;
    }


    /**
     * determine postorder successor of a node
     *
     * @return next node
     */
    public static NodeRef postorderSuccessor(Tree tree, NodeRef node) {

        NodeRef cn = null;
        NodeRef parent = tree.getParent(node);

        if (tree.getRoot() == node) {
            cn = node;
        } else {

            // Go up one node
            if (tree.getChild(parent, tree.getChildCount(parent) - 1) == node) {
                return parent;
            }

            // Go down one node
            for (int i = 0; i < tree.getChildCount(parent) - 1; i++) {
                if (tree.getChild(parent, i) == node) {
                    cn = tree.getChild(parent, i + 1);
                    break;
                }
            }
        }

        // Go down until leaf
        while (tree.getChildCount(cn) > 0) {
            cn = tree.getChild(cn, 0);
        }

        return cn;
    }

    /**
     * Gets finds the most ancestral node with attribute set.
     */
    public static NodeRef findNodeWithAttribute(Tree tree, String attribute) {

        NodeRef root = tree.getRoot();
        NodeRef node = root;

        do {

            if (tree.getNodeAttribute(node, attribute) != null) {
                return node;
            }

            node = TreeUtils.preorderSuccessor(tree, node);

        } while (node != root);

        return null;
    }

    /**
     * Gets finds the most recent date amongst the external nodes.
     */
    public static dr.evolution.util.Date findMostRecentDate(Tree tree) {

        dr.evolution.util.Date mostRecent = null;

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            Taxon taxon = tree.getNodeTaxon(tree.getExternalNode(i));

            dr.evolution.util.Date date = (dr.evolution.util.Date) taxon.getAttribute(dr.evolution.util.Date.DATE);
            if ((date != null) && (mostRecent == null || date.after(mostRecent))) {
                mostRecent = date;
            }
        }

        return mostRecent;
    }

    /**
     * Recursive function for constructing a newick tree representation in the given buffer.
     */
    public static String newick(Tree tree) {
        StringBuffer buffer = new StringBuffer();
        newick(tree, tree.getRoot(), true, BranchLengthType.LENGTHS_AS_TIME, null, null, null, null, buffer);
        buffer.append(";");
        return buffer.toString();
    }

    /**
     * @param tree tree to return in newick format
     * @param dp   the decimal places for branch lengths
     * @return a string representation of the tree in newick format with branch lengths expressed with the given
     *         number of decimal places
     */
    public static String newick(Tree tree, int dp) {
        StringBuffer buffer = new StringBuffer();

        // use the English locale to ensure there are no commas in the number!
        NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);
        format.setMaximumFractionDigits(dp);

        newick(tree, tree.getRoot(), true, BranchLengthType.LENGTHS_AS_TIME, format, null, null, null, buffer);
        buffer.append(";");
        return buffer.toString();
    }

    public static String newick(Tree tree, BranchRates branchRates) {
        StringBuffer buffer = new StringBuffer();
        newick(tree, tree.getRoot(), true, BranchLengthType.LENGTHS_AS_SUBSTITUTIONS, null, branchRates, null, null, buffer);
        buffer.append(";");
        return buffer.toString();
    }

    public static String newick(Tree tree,
                                TreeTraitProvider[] treeTraitProviders
    ) {
        StringBuffer buffer = new StringBuffer();
        newick(tree, tree.getRoot(), true, BranchLengthType.LENGTHS_AS_TIME, null, null, treeTraitProviders, null, buffer);
        buffer.append(";");
        return buffer.toString();
    }

    /**
     * Recursive function for constructing a newick tree representation in the given buffer.
     */
    public static String newickNoLengths(Tree tree) {
        StringBuffer buffer = new StringBuffer();
        newick(tree, tree.getRoot(), true, BranchLengthType.NO_BRANCH_LENGTHS, null, null, null, null, buffer);
        buffer.append(";");
        return buffer.toString();
    }

    /**
     * Recursive function for constructing a newick tree representation in the given buffer.
     *
     * @param tree                     The tree
     * @param node                     The node [tree.getRoot()]
     * @param labels                   whether labels or numbers should be used
     * @param lengths                  What type of branch lengths: NO_BRANCH_LENGTHS, LENGTHS_AS_TIME, LENGTHS_AS_SUBSTITUTIONS
     * @param branchRates              An optional BranchRates (or null) used to scale branch times into substitutions
     * @param treeTraitProviders       An array of TreeTraitProvider
     * @param format                   formatter for branch lengths
     * @param idMap                    A map if id names to integers that is used to overide node labels when present
     * @param buffer                   The StringBuffer
     */
    public static void newick(Tree tree, NodeRef node, boolean labels, BranchLengthType lengths, NumberFormat format,
                              BranchRates branchRates,
                              TreeTraitProvider[] treeTraitProviders,
                              Map<String, Integer> idMap, StringBuffer buffer) {

        NodeRef parent = tree.getParent(node);

        if (tree.isExternal(node)) {
            if (!labels) {
                int k = node.getNumber();
                if (idMap != null) {
                    buffer.append(idMap.get(tree.getTaxonId(k)));
                } else {
                    buffer.append((k + 1));
                }
            } else {
                String label = tree.getTaxonId(node.getNumber());
                if (label.contains(" ") || label.contains(":") || label.contains(";") || label.contains(",")) {
                    buffer.append("\"");
                    buffer.append(label);
                    buffer.append("\"");
                } else {
                    buffer.append(label);
                }
            }
        } else {
            buffer.append("(");
            newick(tree, tree.getChild(node, 0), labels, lengths, format,
                    branchRates,
                    treeTraitProviders, idMap,
                    buffer);
            for (int i = 1; i < tree.getChildCount(node); i++) {
                buffer.append(",");
                newick(tree, tree.getChild(node, i), labels, lengths, format,
                        branchRates,
                        treeTraitProviders, idMap,
                        buffer);
            }
            buffer.append(")");
        }

        writeTreeTraits(buffer, tree, node, treeTraitProviders, TreeTrait.Intent.NODE);

        if (parent != null && lengths != BranchLengthType.NO_BRANCH_LENGTHS) {
            buffer.append(":");
            writeTreeTraits(buffer, tree, node, treeTraitProviders, TreeTrait.Intent.BRANCH);

            if (lengths != BranchLengthType.NO_BRANCH_LENGTHS) {
                double length = tree.getNodeHeight(parent) - tree.getNodeHeight(node);
                if (lengths == BranchLengthType.LENGTHS_AS_SUBSTITUTIONS) {
                    if (branchRates == null) {
                        throw new IllegalArgumentException("No BranchRates provided");
                    }
                    length *= branchRates.getBranchRate(tree, node);
                }
                String lengthString;
                if (format != null) {
                    lengthString = format.format(length);
                } else {
                    lengthString = String.valueOf(length);
                }

                buffer.append(lengthString);
            }
        }
    }

    private static void writeTreeTraits(StringBuffer buffer, Tree tree, NodeRef node, TreeTraitProvider[] treeTraitProviders, TreeTrait.Intent intent) {
        if (treeTraitProviders != null) {
            boolean hasAttribute = false;
            for (TreeTraitProvider ttp : treeTraitProviders) {
                TreeTrait[] tts = ttp.getTreeTraits();
                for (TreeTrait treeTrait: tts) {
                    if (treeTrait.getLoggable() && treeTrait.getIntent() == intent) {
                        String value = treeTrait.getTraitString(tree, node);

                        if (value != null) {
                            if (!hasAttribute) {
                                buffer.append("[&");
                                hasAttribute = true;
                            } else {
                                buffer.append(",");
                            }
                            buffer.append(treeTrait.getTraitName());
                            buffer.append("=");
                            buffer.append(value);

//                                if (values.length > 1) {
//                                    buffer.append("{");
//                                    buffer.append(values[0]);
//                                    for (int i = 1; i < values.length; i++) {
//                                        buffer.append(",");
//                                        buffer.append(values[i]);
//                                    }
//                                    buffer.append("}");
//                                } else {
//                                    buffer.append(values[0]);
//                                }
                        }
                    }

                }
            }
            if (hasAttribute) {
                buffer.append("]");
            }
        }
    }
    /**
     * Recursive function for constructing a newick tree representation in the given buffer.
     */
    public static String uniqueNewick(Tree tree, NodeRef node) {
        if (tree.isExternal(node)) {
            //buffer.append(tree.getNodeTaxon(node).getId());
            return tree.getNodeTaxon(node).getId();
        } else {
            StringBuffer buffer = new StringBuffer("(");

            ArrayList<String> subtrees = new ArrayList<String>();
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                subtrees.add(uniqueNewick(tree, child));
            }
            Collections.sort(subtrees);
            for (int i = 0; i < subtrees.size(); i++) {
                buffer.append(subtrees.get(i));
                if (i < subtrees.size() - 1) {
                    buffer.append(",");
                }
            }
            buffer.append(")");

            return buffer.toString();
        }
    }

    /**
     * Recursive function for constructing a newick tree representation in the given buffer.
     */
    public static Tree rotateByName(Tree tree) {

        return new SimpleTree(rotateNodeByName(tree, tree.getRoot()));
    }

    /**
     * Recursive function for constructing a newick tree representation in the given buffer.
     */
    private static SimpleNode rotateNodeByName(Tree tree, NodeRef node) {

        if (tree.isExternal(node)) {
            return new SimpleNode(tree, node);
        } else {

            SimpleNode parent = new SimpleNode(tree, node);

            NodeRef child1 = tree.getChild(node, 0);
            NodeRef child2 = tree.getChild(node, 1);

            String subtree1 = uniqueNewick(tree, child1);
            String subtree2 = uniqueNewick(tree, child2);

            if (subtree1.compareTo(subtree2) > 0) {
                parent.addChild(rotateNodeByName(tree, child2));
                parent.addChild(rotateNodeByName(tree, child1));
            } else {
                parent.addChild(rotateNodeByName(tree, child1));
                parent.addChild(rotateNodeByName(tree, child2));
            }
            return parent;
        }
    }

    public static MutableTree rotateTreeByComparator(Tree tree, Comparator<NodeRef> comparator) {

        return new SimpleTree(rotateTreeByComparator(tree, tree.getRoot(), comparator));
    }

    /**
     * Recursive function for constructing a newick tree representation in the given buffer.
     */
    private static SimpleNode rotateTreeByComparator(Tree tree, NodeRef node, Comparator<NodeRef> comparator) {

        SimpleNode newNode = new SimpleNode();
        newNode.setHeight(tree.getNodeHeight(node));
        newNode.setRate(tree.getNodeRate(node));
        newNode.setId(tree.getTaxonId(node.getNumber()));
        newNode.setNumber(node.getNumber());
        newNode.setTaxon(tree.getNodeTaxon(node));

        if (!tree.isExternal(node)) {

            NodeRef child1 = tree.getChild(node, 0);
            NodeRef child2 = tree.getChild(node, 1);

            if (comparator.compare(child1, child2) > 0) {
                newNode.addChild(rotateTreeByComparator(tree, child2, comparator));
                newNode.addChild(rotateTreeByComparator(tree, child1, comparator));
            } else {
                newNode.addChild(rotateTreeByComparator(tree, child1, comparator));
                newNode.addChild(rotateTreeByComparator(tree, child2, comparator));
            }
        }

        return newNode;
    }

    public static Comparator<NodeRef> createNodeDensityComparator(final Tree tree) {

        return new Comparator<NodeRef>() {

            public int compare(NodeRef node1, NodeRef node2) {
                return getLeafCount(tree, node2) - getLeafCount(tree, node1);
            }

            public boolean equals(NodeRef node1, NodeRef node2) {
                return getLeafCount(tree, node1) == getLeafCount(tree, node2);
            }
        };
    }

    public static Comparator<NodeRef> createNodeDensityMinNodeHeightComparator(final Tree tree) {

        return new Comparator<NodeRef>() {

            public int compare(NodeRef node1, NodeRef node2) {
                int larger = getLeafCount(tree, node1) - getLeafCount(tree, node2);

                if (larger != 0) return larger;

                double tipRecent = getMinNodeHeight(tree, node2) - getMinNodeHeight(tree, node1);
                if (tipRecent > 0.0) return 1;
                if (tipRecent < 0.0) return -1;
                return 0;
            }
        };
    }

    public static boolean allDisjoint(SimpleNode[] nodes) {

        // check with java 1.6
        Set<String>[] ids = new Set[nodes.length];
        for (int k = 0; k < nodes.length; ++k) {
            ids[k] = TreeUtils.getLeafSet(new SimpleTree(nodes[k]));
            for (int j = 0; j < k; ++j) {
                Set<String> intersection = new HashSet<String>(ids[j]);
                intersection.retainAll(ids[k]);
                if (intersection.size() > 0) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Compares 2 trees and returns true if they have the same topology (same taxon
     * order is assumed).
     */
    public static boolean equal(Tree tree1, Tree tree2) {

        return uniqueNewick(tree1, tree1.getRoot()).equals(uniqueNewick(tree2, tree2.getRoot()));
    }

    private static Node convertToJebl(Tree tree, NodeRef node, SimpleRootedTree jtree) {
        if (tree.isExternal(node)) {
            String taxonId = tree.getTaxonId(node.getNumber());
            Node externalNode = jtree.createExternalNode(jebl.evolution.taxa.Taxon.getTaxon(taxonId));
            jtree.setHeight(externalNode, tree.getNodeHeight(node));
            return externalNode;
        }
        List<Node> jchildren = new ArrayList<Node>();
        for (int nc = 0; nc < tree.getChildCount(node); ++nc) {
            NodeRef child = tree.getChild(node, nc);
            Node node1 = convertToJebl(tree, child, jtree);
            jtree.setHeight(node1, tree.getNodeHeight(child));
            jchildren.add(node1);
        }

        return jtree.createInternalNode(jchildren);
    }

    /**
     * Convert from beast tree to JEBL tree.
     * Note that currently only topology and branch lengths are preserved.
     * Can add attributes later if needed.
     *
     * @param tree beast
     * @return jebl tree
     */
    static public SimpleRootedTree asJeblTree(Tree tree) {
        SimpleRootedTree jtree = new SimpleRootedTree();

        convertToJebl(tree, tree.getRoot(), jtree);
        jtree.setHeight(jtree.getRootNode(), tree.getNodeHeight(tree.getRoot()));
        return jtree;
    }

    /**
     * Gets the set of clades in a tree
     *
     * @param tree the tree
     * @return the set of clades
     */
    public static Set<Set<String>> getClades(Tree tree) {
        Set<Set<String>> clades = new HashSet<Set<String>>();
        getClades(tree, tree.getRoot(), null, clades);

        return clades;
    }

    private static void getClades(Tree tree, NodeRef node, Set<String> leaves, Set<Set<String>> clades) {

        if (tree.isExternal(node)) {
            leaves.add(tree.getTaxonId(node.getNumber()));
        } else {

            Set<String> ls = new HashSet<String>();

            for (int i = 0; i < tree.getChildCount(node); i++) {

                NodeRef node1 = tree.getChild(node, i);

                getClades(tree, node1, ls, clades);
            }

            if (leaves != null) {
                // except for the root clade...
                leaves.addAll(ls);
                clades.add(ls);
            }

        }
    }

    /**
     * Tests whether the given tree is compatible with a set of clades
     *
     * @param tree   the test tree
     * @param clades the set of clades
     * @return
     */
    public static boolean isCompatible(Tree tree, Set<Set<String>> clades) {
        return isCompatible(tree, tree.getRoot(), null, clades);
    }

    private static boolean isCompatible(Tree tree, NodeRef node, Set<String> leaves, Set<Set<String>> clades) {
        if (tree.isExternal(node)) {
            leaves.add(tree.getTaxonId(node.getNumber()));
            return true;
        } else {

            Set<String> ls = new HashSet<String>();

            for (int i = 0; i < tree.getChildCount(node); i++) {

                NodeRef node1 = tree.getChild(node, i);

                if (!isCompatible(tree, node1, ls, clades)) {
                    // as soon as we have an incompatibility break out...
                    return false;
                }
            }

            if (leaves != null) {
                // except for the root clade...
                for (Set<String> clade : clades) {
                    Set<String> intersection = new HashSet<String>(clade);
                    intersection.retainAll(ls);

                    if (intersection.size() != 0 &&
                            intersection.size() != ls.size() &&
                            intersection.size() != clade.size()) {
                        return false;
                    }
                }

                leaves.addAll(ls);
            }
        }
        return true;
    }

    public static void correctBranchLengthToGetUltrametricTree(Tree tree, double givenLength) {


    }

    private void setHeight(Tree tree, NodeRef node, double givenLength) {
        if (tree.getChildCount(node) == 0) {

        }
        for (int i = 0; i < tree.getChildCount(node); i++) {

        }
    }

    public enum BranchLengthType {
        NO_BRANCH_LENGTHS, LENGTHS_AS_TIME, LENGTHS_AS_SUBSTITUTIONS
    }

    public static class MissingTaxonException extends Exception {
        /**
         *
         */
        private static final long serialVersionUID = 8468656622238269963L;

        public MissingTaxonException(Taxon taxon) {
            super(taxon.getId());
        }
    }
}
