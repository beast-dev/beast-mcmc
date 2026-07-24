/*
 * ExorcizedTreeModel.java
 *
 * Copyright © 2002-2026, the BEAST Development Team.
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
 */

package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.treedatalikelihood.RateRescalingScheme;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Citation;

import java.util.*;

/**
 * A model component for trees that extracts the 'corporeal' (non ghost) subtree from a TreeModel
 * that has ghost lineages in it. Ghost lineages are tips (and internal nodes) that don't contribute
 * to the tree data likelihood.
 *
 * @author Andrew Rambaut
 */
public class ExorcizedTreeModel extends TreeModel {
    private final List<Taxon> corporealLineages;
    private final Tree hauntedTree;
    private boolean hauntedTreeChanged;
    private NodeRef exorcizedRoot;
    private final List<NodeRef> exorcizedTips = new ArrayList<>();
    private final List<NodeRef> exorcizedNodes = new ArrayList<>();
    private final Map<NodeRef, NodeRef> exorcizedParentMap = new HashMap<>();
    private final Map<NodeRef, NodeRef[]> exorcizedChildMap = new HashMap<>();

    //
    // Public stuff
    //

    public static final String EXORCIZED_TREE_MODEL = "exorcizedTreeModel";

    public ExorcizedTreeModel(TreeModel ghostTreeModel, TaxonList corporealLineages) {
        this(EXORCIZED_TREE_MODEL, ghostTreeModel, corporealLineages);
    }

    public ExorcizedTreeModel(String name, TreeModel hauntedTree, TaxonList corporealLineages) {

        super(name, true);

        this.hauntedTree = hauntedTree;
        this.corporealLineages = corporealLineages.asList();

        for (int i = 0; i < hauntedTree.getExternalNodeCount(); i++) {
            NodeRef tip = hauntedTree.getExternalNode(i);
            if (corporealLineages.asList().contains(hauntedTree.getNodeTaxon(tip))) {
                exorcizedTips.add(tip);
            }
        }

        updateExorcizedTree();

        addModel(hauntedTree);
    }

    private void updateExorcizedTree() {
        exorcizedNodes.clear();
        updateExorcizedTree(hauntedTree, hauntedTree.getRoot());
        hauntedTreeChanged = false;
    }

    private NodeRef updateExorcizedTree(Tree hauntedTree, NodeRef node) {
        if (hauntedTree.isExternal(node)) {
            if (exorcizedTips.contains(node)) {
                return node;
            }
        } else {
            assert hauntedTree.getChildCount(node) == 2;

            NodeRef child1 = updateExorcizedTree(hauntedTree, hauntedTree.getChild(node, 0));
            NodeRef child2 = updateExorcizedTree(hauntedTree, hauntedTree.getChild(node, 1));

            if (child1 != null && child2 != null) {
                // both children have non-ghost tips so this is a exorcized internal node
                exorcizedNodes.add(node);
                exorcizedRoot = node;

                exorcizedParentMap.put(child1, node);
                exorcizedParentMap.put(child2, node);
                exorcizedChildMap.put(node, new NodeRef[]{child1, child2});

                return node;
            }

            if (child1 != null) {
                return child1;
            }

            if (child2 != null) {
                return child2;
            }

        }
        return null;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        assert model == hauntedTree;
        if (object instanceof TreeChangedEvent) {
            final TreeChangedEvent treeChangedEvent = (TreeChangedEvent) object;

            if (treeChangedEvent.isNodeChanged()) {
                fireModelChanged(treeChangedEvent);
            } else if (treeChangedEvent.isTreeChanged()) {
                fireModelChanged(treeChangedEvent);
            }
            // Other event types are ignored (probably trait changes).
        }

        hauntedTreeChanged = true;
    }

    @Override
    public NodeRef[] getNodes() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public NodeRef getRoot() {
        if (hauntedTreeChanged) {
            updateExorcizedTree();
        }
        return exorcizedRoot;
    }

    @Override
    public int getNodeCount() {
        // this needs to be the haunted tree's node count because
        // the node numbers are unchanged and thus go to hauntedTree.getNodeCount() - 1.
        return hauntedTree.getNodeCount();
    }

    @Override
    public NodeRef getNode(int i) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public NodeRef getInternalNode(int i) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public NodeRef getExternalNode(int i) {
        assert !hauntedTreeChanged;
        return exorcizedTips.get(i);
    }

    @Override
    public int getExternalNodeCount() {
        assert !hauntedTreeChanged;
        return exorcizedTips.size();
    }

    @Override
    public int getInternalNodeCount() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Taxon getNodeTaxon(NodeRef node) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getNodeHeight(NodeRef node) {
        assert !hauntedTreeChanged;
        return hauntedTree.getNodeHeight(node);
    }

    @Override
    public double getNodeRate(NodeRef node) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Object getNodeAttribute(NodeRef node, String name) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Iterator getNodeAttributeNames(NodeRef node) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean isExternal(NodeRef node) {
        assert !hauntedTreeChanged;
        return hauntedTree.isExternal(node);
    }

    @Override
    public boolean isRoot(NodeRef node) {
        assert !hauntedTreeChanged;
        return node == exorcizedRoot;
    }

    @Override
    public int getChildCount(NodeRef node) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public NodeRef getChild(NodeRef node, int j) {
        assert !hauntedTreeChanged;
        return exorcizedChildMap.get(node)[j];
    }

    @Override
    public NodeRef getParent(NodeRef node) {
        assert !hauntedTreeChanged;
        return exorcizedParentMap.get(node);
    }

    @Override
    public Taxon getTaxon(int taxonIndex) {
        return corporealLineages.get(taxonIndex);
    }

    // Unused methods
    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // nothing to do
    }

    @Override
    protected void storeState() {
        // nothing to do
    }

    @Override
    protected void restoreState() {
        // nothing to do
    }

    @Override
    protected void acceptState() {
        // nothing to do
    }

    // *****************************************************************
    // Interface MutableTree
    // *****************************************************************

    /**
     * Set a new node as root node.
     */
    @Override
    public void setRoot(NodeRef newRoot) {
        throw new UnsupportedOperationException("ExorcizedTreeModel is immutable");
    }

    @Override
    public void addChild(NodeRef p, NodeRef c) {
        throw new UnsupportedOperationException("ExorcizedTreeModel is immutable");
    }

    @Override
    public void removeChild(NodeRef p, NodeRef c) {
        throw new UnsupportedOperationException("ExorcizedTreeModel is immutable");
    }

    @Override
    public void replaceChild(NodeRef node, NodeRef child, NodeRef newChild) {

    }

    @Override
    public boolean beginTreeEdit() {
        throw new UnsupportedOperationException("ExorcizedTreeModel is immutable");
    }

    @Override
    public void endTreeEdit() {
        throw new UnsupportedOperationException("ExorcizedTreeModel is immutable");
    }

    @Override
    public void adoptTreeStructure(int[] edges, double[] nodeHeights, int[] childOrder, String[] taxaNames) {

    }

    @Override
    public boolean isTreeValid() {
        return false;
    }

    @Override
    public void setNodeHeight(NodeRef n, double height) {
        throw new UnsupportedOperationException("ExorcizedTreeModel is immutable");
    }

    @Override
    public void setNodeHeightQuietly(NodeRef n, double height) {
        throw new UnsupportedOperationException("ExorcizedTreeModel is immutable");
    }

    @Override
    public void setNodeRate(NodeRef n, double rate) {
        throw new UnsupportedOperationException("Function not available in GhostTreeModel");

    }

    @Override
    public void setBranchLength(NodeRef node, double length) {
        throw new UnsupportedOperationException("ExorcizedTreeModel is immutable");
    }

    @Override
    public void setNodeAttribute(NodeRef node, String name, Object value) {
        throw new UnsupportedOperationException("ExorcizedTreeModel is immutable");
    }

    @Override
    public double[] getMultivariateNodeTrait(NodeRef node, String name) {
        throw new UnsupportedOperationException("ExorcizedTreeModel is immutable");
    }

    @Override
    public void setMultivariateTrait(NodeRef n, String name, double[] value) {
        throw new UnsupportedOperationException("ExorcizedTreeModel is immutable");
    }


    @Override
    public List<Citation> getCitations() {
        // @todo add a citation in here
        return Collections.EMPTY_LIST;
    }

}