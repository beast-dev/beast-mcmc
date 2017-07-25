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
public class AncestralTraitTreeModel extends AbstractModel implements MultivariateTraitTree, Citable {

    private final int ancestorCount;

    private final int treeExternalCount;
    private final int treeInternalCount;

    private final int externalCount;
    private final int internalCount;

    private ShadowNode[] nodes;
    private ShadowNode[] storedNodes;

    private ShadowNode root;
    private int storedRootNumber;

//    private final Map<NodeRef, ShadowNode> originalToShadowMap = new HashMap<NodeRef, ShadowNode>();

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
    }

    public AncestralTraitTreeModel(String id, TreeModel tree, List<AncestralTaxonInTree> ancestors) {
        super(id);
        this.treeModel = tree;
        this.ancestors = ancestors;

        ancestorCount = ancestors.size();

        treeExternalCount = treeModel.getExternalNodeCount();
        treeInternalCount = treeModel.getInternalNodeCount();

        externalCount = treeExternalCount + ancestorCount;
        internalCount = treeInternalCount + ancestorCount;

        addModel(tree);

        nodes = new ShadowNode[externalCount + internalCount];

        int index = 0;
        for (AncestralTaxonInTree ancestor : ancestors) {
            addRestrictedPartials(ancestor, index);
            ++index;
        }

        buildShadowTree(treeModel);
    }

    private void checkShadowTree() {

    }

    private void buildShadowTree(Tree original) { // TODO Cache!

        setupClamps();

//        originalToShadowMap.clear();

        root = buildRecursivelyShadowTree(treeModel.getRoot(), null,
                 0);
    }

    private void storeNode(ShadowNode node, NodeRef originalNode) {
//        System.err.println("Storing node # " + node.getNumber());
        nodes[node.getNumber()] = node;
//        if (originalNode != null) {
//            originalToShadowMap.put(originalNode, node);
//            System.err.println("Stored node #" + originalNode.getNumber() + " into #" + node.getNumber());
//        }
    }

    private ShadowNode buildRecursivelyShadowTree(NodeRef originalNode,
                                                  ShadowNode parentNode,
                                                  int extraInternal) {
        final int originalNumber = originalNode.getNumber();
        final int newNumber = mapOriginalToShadowNumber(originalNumber);

//        Taxon taxon = (newNumber < externalCount) ?
//                (originalNumber < treeExternalCount ?
//                treeModel.getNodeTaxon(originalNode) :
//                        getTaxonByTreeIndex(originalNumber)) :
//                null;

        ShadowNode newNode = new ShadowNode(newNumber, originalNode, null);
        newNode.parent = parentNode;
        storeNode(newNode, originalNode);

        ShadowNode returnNode = newNode;

        if (nodeToClampMap.containsKey(originalNode)) {
//            System.err.println("Found a clamp!");

//            storeNode(newNode, originalNode);

            // Add tip
            AncestralTaxonInTree ancestor = nodeToClampMap.get(originalNode);
            final int newTipNumber = treeExternalCount + ancestor.getIndex();

            ShadowNode newTipNode = new ShadowNode(newTipNumber, null, ancestor);
//            newTipNode.taxon = ancestor.getTaxon();

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

        
//        nodes[newNumber] = newNode;
//        storeNode(newNode);

        if (treeModel.isExternal(originalNode)) {

//            newNode.taxon = treeModel.getTaxon(originalNumber);

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

//    private int mapShadowNumberToOriginal(int newNumber) {
//      return 0;
//    }

//    private buildAugmentedTree()

    public String toString() {
        return TreeUtils.newick(this);
    }

    public boolean isVariable() {
        return treeModel.isVariable();
    }

    public double getNodeHeight(NodeRef inode) {

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

    public double getBranchLength(NodeRef inode) {

        ShadowNode node = (ShadowNode) inode;
        NodeRef originalNode = ((ShadowNode) node).getOriginalNode();
        if (originalNode != null) {
            return treeModel.getBranchLength(originalNode);
        } else {
//            double length = treeModel.getBranchLength(node.parent.originalNode);
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
            for (int i = 0; i < length; ++i) {
                storedNodes[i] = new ShadowNode();
            }
        }

        // Copy
        for (int i = 0; i < length; ++i) {
            storedNodes[i].adoptValues(nodes[i]);
        }
    }
    
    protected void storeState() {

        assert (nodes != null);
        storeNodeStructure();

        assert (nodes[root.getNumber()] == root);
        storedRootNumber = root.getNumber();
    }

    /**
     * Restore the stored state
     */
    protected void restoreState() {

        assert (storedNodes != null);
        ShadowNode[] tmp = storedNodes;
        storedNodes = nodes;
        nodes = tmp;

        root = nodes[storedRootNumber];
    }

    class AncestralTreeChangeEvent implements TreeChangedEvent {

        final private TreeChangedEvent event;
        final private NodeRef node;

        private AncestralTreeChangeEvent(TreeChangedEvent event, NodeRef node) {
            this.event = event;
            this.node = node;
        }

//        private AncestralTreeChangeEvent(TreeChangedEvent event, NodeRef originalNode, boolean propogate) {
//            this.event = event;
////            System.err.println("Original # is " + originalNode.getNumber());
////            this.node = originalToShadowMap.get(originalNode);
//            ShadowNode shadow = nodes[mapOriginalToShadowNumber(originalNode.getNumber())];
////            System.err.println("Shadow # is " + shadow.getNumber());
//
//            if (!shadow.isExternal() && propogate) {
//                if (shadow.child1.ancestor != null) {
//                    shadow = shadow.child0;
//                }
//            }
////            if (shadow.child1.ancestor != null) {
////                System.err.println("Yep");
////            } else {
////                System.err.println("Nope");
////            }
////            if (shadow.ancestor != null) {
////                shadow = shadow.child0;
////            }
//
//            this.node = shadow;
////            System.err.println("New node # is " + node.getNumber());
////            System.exit(-1);
//        }

//        private AncestralTreeChangeEvent(TreeChangedEvent event, final int number) {
//            this.event = event;
//            this.node = new NodeRef() {
//
//                @Override public int getNumber() { return number; }
//
//                @Override public void setNumber(int n) {
//                    throw new IllegalArgumentException("Not mutable");
//                }
//            };
//        }

        @Override public int getIndex() { return event.getIndex(); }

        @Override public NodeRef getNode() { return node; }

        @Override public Parameter getParameter() { return event.getParameter(); }

        @Override public boolean isNodeChanged() { return event.isNodeChanged(); }

        @Override public boolean isTreeChanged() { return event.isTreeChanged(); }

//        @Override public boolean isHeightChanged() { return event.isHeightChanged(); }
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

        if (model == treeModel) {

            if (object instanceof TreeChangedEvent) {
                final TreeChangedEvent treeChangedEvent = (TreeChangedEvent) object;

                if (treeChangedEvent.isNodeChanged()) {
                    System.err.println("node changed");

                    final NodeRef originalNode = treeChangedEvent.getNode();
                    if (originalNode != null) { // Remap

                        ShadowNode shadow = nodes[mapOriginalToShadowNumber(originalNode.getNumber())];
                        if (!shadow.isExternal()) {
                            if (shadow.child1.ancestor != null) {

                                fireModelChanged(new AncestralTreeChangeEvent(treeChangedEvent, shadow), index);
                                shadow = shadow.child0;
                            }
                        }

                        object = new AncestralTreeChangeEvent(treeChangedEvent, shadow);

                    }

                    fireModelChanged(object, index);

                } else if (treeChangedEvent.isTreeChanged()) {
                    System.err.println("tree changed!");
                    fireModelChanged(object, index);
                } else {
                    throw new IllegalArgumentException("Illegal");
                }

            } else if (object instanceof Parameter) {
                // Do nothing
            } else {
                throw new IllegalArgumentException("TreeModel should not generate other events");
            }

        } else if (ancestors.contains(model)) {

            AncestralTaxonInTree ancestor = (AncestralTaxonInTree) model;
            NodeRef node = nodes[treeExternalCount + ancestor.getIndex()];

            assert (nodes[treeExternalCount + ancestor.getIndex()].ancestor == ancestor);
            fireModelChanged(treeModel.createTreeChangeEvent());

        } else {
            throw new IllegalArgumentException("Illegal model");
        }
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        // Do nothing; no variables
    }

    // Delegate the rest to treeModel

    public NodeRef getRoot() {
          return root;
//        return treeModel.getRoot();
    }

    public int getNodeCount() {
        return externalCount + internalCount;
    }

    public NodeRef getNode(int i) {
        return nodes[i];
//        return treeModel.getNode(i);
    }

    public NodeRef getInternalNode(int i) {
        return nodes[i + externalCount];
    }

    public NodeRef getExternalNode(int i) {
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
//        return treeModel.getNodeRate(node);
    }

    public Object getNodeAttribute(NodeRef node, String name) {
        throw new RuntimeException("Not yet implemented");
//        return treeModel.getNodeAttribute(node, name);
    }

    public Iterator getNodeAttributeNames(NodeRef node) {
        throw new RuntimeException("Not yet implemented");
//        return treeModel.getNodeAttributeNames(node);
    }

    public boolean isExternal(NodeRef node) {
        return ((ShadowNode) node).isExternal();
//        throw new RuntimeException("Not yet implemented");
//        return treeModel.isExternal(node);
    }

    public boolean isRoot(NodeRef node) {
        return node == root;
    }

    public int getChildCount(NodeRef node) {

        if (((ShadowNode) node).isExternal()) {
            return 0;
        } else {
            return 2;
        }
    }

    public NodeRef getChild(NodeRef node, int j) {
        return ((ShadowNode) node).getChild(j);
//        throw new RuntimeException("Not yet implemented");
//        return treeModel.getChild(node, j);
    }

    public NodeRef getParent(NodeRef node) {
        assert (node != null);

        return ((ShadowNode) node).parent;
//        throw new RuntimeException("Not yet implemented");
//        return treeModel.getParent(node);
    }

    public Tree getCopy() {
        throw new RuntimeException("Not yet implemented");
        //return treeModel.getCopy();
    }

    protected void acceptState() {
        // Do nothing
    }

    private final TreeModel treeModel;
    private final List<AncestralTaxonInTree> ancestors;

    public double[] getMultivariateNodeTrait(NodeRef node, String name) {
        throw new RuntimeException("Not yet implemented");
//        return treeModel.getMultivariateNodeTrait(node, name);
    }

    public void setMultivariateTrait(NodeRef n, String name, double[] value) {
        throw new RuntimeException("Not yet implemented");
//        treeModel.setMultivariateTrait(n, name, value);
    }

    public boolean beginTreeEdit() {
        return treeModel.beginTreeEdit();
    }

    public void endTreeEdit() {
        treeModel.endTreeEdit();
    }

    public void addChild(NodeRef parent, NodeRef child) {
        treeModel.addChild(parent, child);
    }

    public void removeChild(NodeRef parent, NodeRef child) {
        treeModel.removeChild(parent, child);
    }

    public void replaceChild(NodeRef node, NodeRef child, NodeRef newChild) {
        treeModel.replaceChild(node, child, newChild);
    }

    public void setRoot(NodeRef root) {
        treeModel.setRoot(root);
    }

    public void setNodeHeight(NodeRef node, double height) {
        treeModel.setNodeHeight(node, height);
    }

    public void setNodeRate(NodeRef node, double height) {
       treeModel.setNodeRate(node, height);
    }

    public void setBranchLength(NodeRef node, double length) {
        treeModel.setBranchLength(node, length);
    }

    public void setNodeAttribute(NodeRef node, String name, Object value) {
        treeModel.setNodeAttribute(node, name, value);
    }

    public void addMutableTreeListener(MutableTreeListener listener) {
        treeModel.addMutableTreeListener(listener);
    }

    public void setAttribute(String name, Object value) {
        treeModel.setAttribute(name, value);
    }

    public Object getAttribute(String name) {
        return treeModel.getAttribute(name);
    }

    public Iterator<String> getAttributeNames() {
        return treeModel.getAttributeNames();
    }

    public int addTaxon(Taxon taxon) {
        return treeModel.addTaxon(taxon);
    }

    public boolean removeTaxon(Taxon taxon) {
        return treeModel.removeTaxon(taxon);
    }

    public void setTaxonId(int taxonIndex, String id) {
        treeModel.setTaxonId(taxonIndex, id);
    }

    public void setTaxonAttribute(int taxonIndex, String name, Object value) {
        treeModel.setTaxonAttribute(taxonIndex, name, value);
    }

    public void addMutableTaxonListListener(MutableTaxonListListener listener) {
        treeModel.addMutableTaxonListListener(listener);
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
        int index = treeModel.getTaxonIndex(id);
        if (index != -1) {
            return index;
        } else {
            throw new RuntimeException("Not yet implemented");
        }
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

        System.err.println("Added a CLAMP #" + nodeClamp.getIndex() + " for " + nodeClamp.getTaxon().getId());
    }

    final private Map<BitSet, AncestralTaxonInTree> clampList = new HashMap<BitSet, AncestralTaxonInTree>();
    final private Map<NodeRef, AncestralTaxonInTree> nodeToClampMap = new HashMap<NodeRef, AncestralTaxonInTree>();

}
