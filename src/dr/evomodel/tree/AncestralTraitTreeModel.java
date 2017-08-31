/*
 * TransformedTreeModel.java
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

package dr.evomodel.tree;

import dr.evolution.tree.*;
import dr.evolution.util.MutableTaxonListListener;
import dr.evolution.util.Taxon;
import dr.evomodel.continuous.AncestralTaxonInTree;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.*;

/**
 * A tree model with additional ancestral taxon nodes for peeling
 *
 * @author Marc A. Suchard
 */
public class AncestralTraitTreeModel extends AbstractModel implements MutableTreeModel, Citable {

    private static final boolean DEBUG = false;

    private final int ancestorCount;

    private final int treeExternalCount;
    private final int treeInternalCount;

    private final int externalCount;
    private final int internalCount;

    private ShadowNode[] nodes;
    private ShadowNode[] storedNodes;

    private ShadowNode root;
    private ShadowNode storedRoot;

    public class ShadowNode implements NodeRef {

        private int number;
        private NodeRef originalNode;
        private AncestralTaxonInTree ancestor;

        protected ShadowNode child0 = null;
        protected ShadowNode child1 = null;
        protected ShadowNode parent = null;

        public ShadowNode() { /* Do nothing */ }

        public ShadowNode(int number, NodeRef originalNode, AncestralTaxonInTree ancestor) {
            this.number = number;
            this.originalNode = originalNode;
            this.ancestor = ancestor;
        }

        public void adoptValues(ShadowNode donor) {
            this.number = donor.number;
            this.originalNode = donor.originalNode;
            this.ancestor = donor.ancestor;
            this.child0 = donor.child0;
            this.child1 = donor.child1;
            this.parent = donor.parent;
        }

        @Override
        public int getNumber() { return number; }

        @Override
        public void setNumber(int n) {
            throw new RuntimeException("Node number is not modifiable");
        }

        public NodeRef getOriginalNode() { return originalNode; }

        public NodeRef getChild(int i) {
            if (i == 0) {
                return child0;
            } else if (i == 1) {
                return child1;
            } else {
                throw new IllegalArgumentException("Binary trees only!");
            }
        }

        public boolean isExternal() { return child0 == null && child1 == null; }

        public String toString() {
            return "node " + number + " new";
        }
    }

    public AncestralTraitTreeModel(String id, MutableTreeModel tree, List<AncestralTaxonInTree> ancestors) {
        super(id);
        this.treeModel = tree;
        this.ancestors = ancestors;

        ancestorCount = ancestors.size();

        treeExternalCount = treeModel.getExternalNodeCount();
        treeInternalCount = treeModel.getInternalNodeCount();

        externalCount = treeExternalCount + ancestorCount;
        internalCount = treeInternalCount + ancestorCount;

        addModel(tree);

        int index = 0;
        for (AncestralTaxonInTree ancestor : ancestors) {
            addRestrictedPartials(ancestor, index);
            ++index;
        }
    }

    private void checkShadowTree() {
        if (!validShadowTree) {
            buildShadowTree();
            validShadowTree = true;
        }
    }

    private void buildShadowTree() {

        if (DEBUG) {
            System.out.println("ATTM.bST");
        }

        setupClamps();

        nodes = new ShadowNode[externalCount + internalCount];

        root = buildRecursivelyShadowTree(treeModel.getRoot(), null,
                 0);

//        if (DEBUG) {
//            if (nodeToClampMap.size() != ancestors.size()) {
//                StringBuilder sb = new StringBuilder();
//                for (AncestralTaxonInTree ancestor : ancestors) {
//                    Set<String> names = new HashSet<String>();
//                    if (!nodeToClampMap.containsKey(ancestor.getNode())) {
//                        sb.append("Unable to find ancestor '" + ancestor.getTaxon() +
//                                "' in treeModel '" + treeModel.getId() + "'");
//                    }
//                    sb.append("\n                       ");
//                    for (Taxon taxon : ancestor.getTaxonList()) {
//                        sb.append(" " + taxon.getId());
//                        names.add(taxon.getId());
//                    }
//
//                    boolean mono = TreeUtils.isMonophyletic(treeModel, names);
//                    sb.append("\nis mono = " + mono);
//                }
//
//                System.out.println(sb.toString());
//
////            throw new RuntimeException(sb.toString());
//            }
//        }

//        // Check that all nodes are non-null
//        for (int i = 0; i < nodes.length; ++i) {
//            if (nodes[i] == null) {
//                System.err.println(treeModel.getExternalNodeCount());
//                System.err.println(treeModel.getInternalNodeCount());
//                throw new RuntimeException("Node " + i + " was uninitialized" +
//                "\n" + externalCount + "\n" + internalCount);
//            }
//        }

        validShadowTree = true;

//        if (DEBUG) {
//
//            if (!checkNegativeBranchLength(treeModel)) {
//                System.out.println("Error in TM");
//            }
//
//            if (!checkNegativeBranchLength(this)) {
//                System.out.printf("Error in ATTM");
//            }
//
//            String s1 = TreeUtils.newick(treeModel);
//            String s2 = TreeUtils.newick(this);
//
//            System.out.println("tm : " + s1);
//            System.out.println("     " + TreeUtils.uniqueNewick(treeModel, treeModel.getRoot()));
//            System.out.println("at : " + s2);
//            System.out.println("     " + TreeUtils.uniqueNewick(this, getRoot()));
//
//
//            if (s1.compareTo(s2) != 0.0) {
//                System.out.println("unequal trees ???");
//                System.out.println(s1);
//                System.out.println(s2);
//            } else {
//                System.out.println("equal trees");
//            }
//        }
    }

//    boolean checkNegativeBranchLength(Tree tree) {
//        for (int i = 0; i < tree.getNodeCount(); ++i) {
//            NodeRef node = tree.getNode(i);
//            if (!tree.isRoot(node)) {
//                double a = tree.getNodeHeight(node);
//                double b = tree.getNodeHeight(tree.getParent(node));
//                if (b < a) {
//                    return false;
//                }
//            }
//        }
//
//        return true;
//    }
    
    private void storeNode(ShadowNode node, NodeRef originalNode) {
        nodes[node.getNumber()] = node;
    }

    private ShadowNode buildRecursivelyShadowTree(NodeRef originalNode,
                                                  ShadowNode parentNode,
                                                  int extraInternal) {
        final int originalNumber = originalNode.getNumber();
        final int newNumber = mapOriginalToShadowNumber(originalNumber);

        ShadowNode newNode = new ShadowNode(newNumber, originalNode, null);
        newNode.parent = parentNode;
        storeNode(newNode, originalNode);

        ShadowNode returnNode = newNode;

        if (nodeToClampMap.containsKey(originalNode)) {

            // Add tip
            AncestralTaxonInTree ancestor = nodeToClampMap.get(originalNode);
            final int newTipNumber = treeExternalCount + ancestor.getIndex();

            ShadowNode newTipNode = new ShadowNode(newTipNumber, null, ancestor);

            newTipNode.parent = newNode;
            newNode.child1 = newTipNode;

            storeNode(newTipNode, null);

            ShadowNode newInternalNode = new ShadowNode(externalCount + treeInternalCount + extraInternal, null, null);

            newInternalNode.parent = newNode;
            newNode.child0 = newInternalNode;

            storeNode(newInternalNode, null);

            ++extraInternal;

            newNode = newInternalNode;
        }

        if (treeModel.isExternal(originalNode)) {
            // Do nothing
        } else {
            newNode.child0 = buildRecursivelyShadowTree(treeModel.getChild(originalNode, 0),
                    newNode, extraInternal);
            newNode.child1 = buildRecursivelyShadowTree(treeModel.getChild(originalNode, 1),
                    newNode, extraInternal);
        }

        return returnNode;
    }

    private int mapOriginalToShadowNumber(int originalNumber) {

        assert (originalNumber >= 0);
        assert (originalNumber < treeExternalCount + treeInternalCount);

        int newNumber = originalNumber;
        if (originalNumber >= treeExternalCount) {
            newNumber += ancestorCount;
        }
        return newNumber;
    }

    public String toString() {
        return TreeUtils.newick(this);
    }

    public boolean isVariable() {
        if (treeModel instanceof AbstractModel) {
            return ((AbstractModel) treeModel).isVariable();
        } else {
            return false;
        }
    }

    public double getNodeHeight(NodeRef inode) {

        assert (inode != null);
        checkShadowTree();

//        final int n = inode.getNumber();
//
//        return treeModel.getNodeHeight(getTreeModelNode(n));
        
        ShadowNode node = (ShadowNode) inode;
        NodeRef originalNode = node.getOriginalNode();
        if (originalNode != null) {
            return treeModel.getNodeHeight(originalNode);
        } else {

            double height = treeModel.getNodeHeight(node.parent.originalNode);
            if (node.isExternal()) {
                double diff = node.ancestor.getPseudoBranchLength();
                return height - diff;
            } else {
                return height;
            }
        }
    }

//    private NodeRef getTreeModelNode(int n) {
//        return treeModel.getNode(n);
//    }

    public double getBranchLength(NodeRef inode) {

        checkShadowTree();

        if (inode == null) {
            return 0.0;
        }

        ShadowNode node = (ShadowNode) inode;
        NodeRef originalNode = ((ShadowNode) node).getOriginalNode();
        if (originalNode != null) {
            return treeModel.getBranchLength(originalNode);
        } else {
            if (node.isExternal()) {
                return node.ancestor.getPseudoBranchLength();
            } else {
                return 0;
            }
        }
    }

    private void storeNodeStructure() {

        final int length = nodes.length;

        // Initialize once
        if (storedNodes == null) {
            storedNodes = new ShadowNode[length];
        }

        // Copy
        for (int i = 0; i < length; ++i) {

            if (nodes[i] == null) {
                storedNodes[i] = null;
            } else {
                if (storedNodes[i] == null) {
                    storedNodes[i] = new ShadowNode();
                }
                storedNodes[i].adoptValues(nodes[i]);
            }
        }
    }
    
    protected void storeState() {

        assert (nodes != null);

        savedValidShadowTree = validShadowTree;

        if (validShadowTree) {
            storeNodeStructure();
        }

        storedRoot = root;
    }

    /**
     * Restore the stored state
     */
    protected void restoreState() {

        assert (storedNodes != null);

//        validShadowTree = savedValidShadowTree;
//
//        if (validShadowTree) { // TODO Add back
//            ShadowNode[] tmp = storedNodes;
//            storedNodes = nodes;
//            nodes = tmp;
//        }
//
//        root = storedRoot;

        validShadowTree = false;
    }

    private class RemappedTreeChangeEvent implements TreeChangedEvent {

        final private TreeChangedEvent event;
        final private NodeRef node;

        private RemappedTreeChangeEvent(TreeChangedEvent event, NodeRef node) {
            this.event = event;
            this.node = node;
        }

        @Override public int getIndex() { return event.getIndex(); }

        @Override public NodeRef getNode() { return node; }

        @Override public Parameter getParameter() { return event.getParameter(); }

        @Override public boolean isNodeChanged() { return event.isNodeChanged(); }

        @Override public boolean isTreeChanged() { return event.isTreeChanged(); }
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

        if (DEBUG) {
            System.out.println("ATTM.hMCE hit with " + object.getClass().getCanonicalName());
        }

        if (model == treeModel) {

            if (object instanceof TreeChangedEvent) {
                final TreeChangedEvent treeChangedEvent = (TreeChangedEvent) object;

                if (treeChangedEvent.isTreeChanged()) {

                    if (DEBUG) {
                        System.out.println("ATTM.hMCE invalidate");
                    }

                    validShadowTree = false;
                    fireModelChanged(new TreeChangedEvent.WholeTree());

                } else if (treeChangedEvent.isNodeChanged()) {

                    if (DEBUG) {
                        System.out.println("ATTM.hMCE nodeChange");
                    }

                    final NodeRef originalNode = treeChangedEvent.getNode();
                    if (originalNode != null) { // Remap
                        ShadowNode shadow = nodes[mapOriginalToShadowNumber(originalNode.getNumber())];

                        if (DEBUG) {
                            System.out.println("ATTM.hMCE isRoot? " + isRoot(shadow) + " " + treeModel.isRoot(originalNode));
                        }

                        if (!shadow.isExternal()) {
                            if (shadow.child1.ancestor != null) {
                                /* If there is an ancestor (and zero-branch-length internal node),
                                 * we apparently need to hit the both the parent of the ancestor and the
                                 * zero-branch-length internal node.  MAS is unsure why we must hit the
                                 * ancestor parent.
                                 */
                                fireModelChanged(new RemappedTreeChangeEvent(treeChangedEvent, shadow), index);
                                shadow = shadow.child0;  // Get zero-branch-length internal node
                            }
                        }

                        object = new RemappedTreeChangeEvent(treeChangedEvent, shadow);

                    }
                    fireModelChanged(object, index);

                } else {
                    throw new IllegalArgumentException("TreeModel should not generate other events");
                }

            } else if (object instanceof Parameter) {
                // Do nothing
            } else {
                throw new IllegalArgumentException("TreeModel should not generate other objects");
            }

        } else if (ancestors.contains(model)) {

            AncestralTaxonInTree ancestor = (AncestralTaxonInTree) model;

            assert (nodes[treeExternalCount + ancestor.getIndex()].ancestor == ancestor);

            fireModelChanged(new TreeChangedEvent.WholeTree());

        } else {
            throw new IllegalArgumentException("Illegal model");
        }
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        // Do nothing; no variables
    }

    // Delegate the rest to treeModel

    public NodeRef getRoot() {
        checkShadowTree();
        return root;
    }

    public int getNodeCount() {
        return externalCount + internalCount;
    }

    public NodeRef getNode(int i) {
        assert (nodes[i].originalNode.getNumber() == nodes[i].getNumber()); // TODO Remove
        checkShadowTree();
        return nodes[i];
    }

    public NodeRef getInternalNode(int i) {
        checkShadowTree();
        NodeRef node = nodes[i + externalCount];
        return node;
    }

    public NodeRef getExternalNode(int i) {
        checkShadowTree();
        NodeRef node =  nodes[i];
        return node;
    }

    public int getExternalNodeCount() {
        return externalCount;
    }

    public int getInternalNodeCount() {
        return internalCount;
    }

    public Taxon getNodeTaxon(NodeRef node) {
        checkShadowTree();
        NodeRef originalNode = ((ShadowNode) node).getOriginalNode();
        if (originalNode != null) {
            return treeModel.getNodeTaxon(originalNode);
        } else {
            return getTaxonByTreeIndex(node.getNumber());
        }
    }

    public boolean hasNodeHeights() {
        return treeModel.hasNodeHeights();
    }

    public boolean hasBranchLengths() {
        return treeModel.hasBranchLengths();
    }

    public double getNodeRate(NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    public Object getNodeAttribute(NodeRef node, String name) {
        throw new RuntimeException("Not yet implemented");
    }

    public Iterator getNodeAttributeNames(NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    public boolean isExternal(NodeRef node) {
        checkShadowTree();
        return ((ShadowNode) node).isExternal();
    }

    public boolean isRoot(NodeRef node) {
        checkShadowTree();
        return node == root;
    }

    public int getChildCount(NodeRef node) {
        checkShadowTree();
        if (((ShadowNode) node).isExternal()) {
            return 0;
        } else {
            return 2;
        }
    }

    public NodeRef getChild(NodeRef node, int j) {
        assert (node != null);
        
        checkShadowTree();
        return ((ShadowNode) node).getChild(j);
    }

    public NodeRef getParent(NodeRef node) {
        assert (node != null);
        
        checkShadowTree();
        return ((ShadowNode) node).parent;
    }

    public Tree getCopy() {
        throw new RuntimeException("Not yet implemented");
    }

    protected void acceptState() {
        // Do nothing
    }

    private final MutableTreeModel treeModel;
    private final List<AncestralTaxonInTree> ancestors;

    public double[] getMultivariateNodeTrait(NodeRef node, String name) {
        throw new RuntimeException("Not yet implemented");
    }

    public void setMultivariateTrait(NodeRef n, String name, double[] value) {
        throw new RuntimeException("Not yet implemented");
    }

    public boolean beginTreeEdit() {
        throw new RuntimeException("Not yet implemented");
    }

    public void endTreeEdit() {
        throw new RuntimeException("Not yet implemented");
    }

    public void addChild(NodeRef parent, NodeRef child) {
        throw new RuntimeException("Not yet implemented");
    }

    public void removeChild(NodeRef parent, NodeRef child) {
        throw new RuntimeException("Not yet implemented");
    }

    public void replaceChild(NodeRef node, NodeRef child, NodeRef newChild) {
        throw new RuntimeException("Not yet implemented");
    }

    public void setRoot(NodeRef root) {
        throw new RuntimeException("Not yet implemented");
    }

    public void setNodeHeight(NodeRef node, double height) {
        throw new RuntimeException("Not yet implemented");
    }

    public void setNodeRate(NodeRef node, double height) {
        throw new RuntimeException("Not yet implemented");
    }

    public void setBranchLength(NodeRef node, double length) {
        throw new RuntimeException("Not yet implemented");
    }

    public void setNodeAttribute(NodeRef node, String name, Object value) {
        throw new RuntimeException("Not yet implemented");
    }

    public void addMutableTreeListener(MutableTreeListener listener) {
        throw new RuntimeException("Not yet implemented");
    }

    public void setAttribute(String name, Object value) {
        throw new RuntimeException("Not yet implemented");
    }

    public Object getAttribute(String name) {
        throw new RuntimeException("Not yet implemented");
    }

    public Iterator<String> getAttributeNames() {
        throw new RuntimeException("Not yet implemented");
    }

    public int addTaxon(Taxon taxon) {
        throw new RuntimeException("Not yet implemented");
    }

    public boolean removeTaxon(Taxon taxon) {
        throw new RuntimeException("Not yet implemented");
    }

    public void setTaxonId(int taxonIndex, String id) {
        treeModel.setTaxonId(taxonIndex, id);
        throw new RuntimeException("Not yet implemented");
    }

    public void setTaxonAttribute(int taxonIndex, String name, Object value) {
        throw new RuntimeException("Not yet implemented");
    }

    public void addMutableTaxonListListener(MutableTaxonListListener listener) {
        throw new RuntimeException("Not yet implemented");
    }

    public int getTaxonCount() {
        return treeModel.getTaxonCount() + ancestorCount;
    }

    public Taxon getTaxon(int taxonIndex) {
        Taxon taxon;
        if (taxonIndex < treeExternalCount) {
            taxon = treeModel.getTaxon(taxonIndex);
        } else {
            taxon = getTaxonByTreeIndex(taxonIndex);
        }
        return taxon;
    }

    public String getTaxonId(int taxonIndex) {
        String result;
        if (taxonIndex < treeExternalCount) {
            result = treeModel.getTaxonId(taxonIndex);
        } else {
            Taxon taxon = getTaxonByTreeIndex(taxonIndex);
            result = taxon.getId();
        }
//        System.err.println("gTI: " + result + " for " + taxonIndex);
        return result;
    }

    public int getTaxonIndex(String id) {
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            if (getTaxonId(i).equals(id)) return i;
        }
        return -1;
    }

    public int getTaxonIndex(Taxon taxon) {

        int index = treeModel.getTaxonIndex(taxon);
        if (index != -1) {
            return index;
        } else {
            throw new RuntimeException("Not yet implemented");
        }
    }

    public List<Taxon> asList() {
        throw new RuntimeException("Not yet implemented");
//        return treeModel.asList();
    }

    public Object getTaxonAttribute(int taxonIndex, String name) {
        if (taxonIndex < treeModel.getExternalNodeCount()) {
            return treeModel.getTaxonAttribute(taxonIndex, name);
        } else {
            Taxon taxon = getTaxonByTreeIndex(taxonIndex);
            if (taxon != null) {
                return taxon.getAttribute(name);
            }
            return null;
        }
    }

    public Iterator<Taxon> iterator() {
        throw new RuntimeException("Not yet implemented");
//        return treeModel.iterator();
    }

    public Type getUnits() {
        return treeModel.getUnits();
    }

    public void setUnits(Type units) {
        treeModel.setUnits(units);
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TRAIT_MODELS;
    }

    @Override
    public String getDescription() {
        return "Bayesian estimation of Pagel's lambda";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.VRANCKEN_2015_SIMULTANEOUSLY);
    }

    private Taxon getTaxonByTreeIndex(int index) {
        return ancestors.get(index - treeExternalCount).getTaxon();
    }

    public void setupClamps() {
        nodeToClampMap.clear();
        recursiveSetupClamp(treeModel, treeModel.getRoot(), new BitSet());
    }

    private void recursiveSetupClamp(Tree tree, NodeRef node, BitSet tips) {

        if (tree.isExternal(node)) {
            tips.set(node.getNumber());
        } else {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);

                BitSet childTips = new BitSet();
                recursiveSetupClamp(tree, child, childTips);
                tips.or(childTips);
            }

            if (clampList.containsKey(tips)) {
                AncestralTaxonInTree partials = clampList.get(tips);
                partials.setNode(node);
                nodeToClampMap.put(node, partials);
            }
        }
    }

    private void addRestrictedPartials(AncestralTaxonInTree nodeClamp, int index) {

        clampList.put(nodeClamp.getTipBitSet(), nodeClamp);
        addModel(nodeClamp);
        nodeClamp.setIndex(index);

//        System.err.println("Added a CLAMP #" + nodeClamp.getIndex() + " for " + nodeClamp.getTaxon().getId());
    }

    final private Map<BitSet, AncestralTaxonInTree> clampList = new HashMap<BitSet, AncestralTaxonInTree>();
    final private Map<NodeRef, AncestralTaxonInTree> nodeToClampMap = new HashMap<NodeRef, AncestralTaxonInTree>();

    private boolean validShadowTree = false;
    private boolean savedValidShadowTree;

}
