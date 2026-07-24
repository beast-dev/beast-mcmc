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

    //
    // Public stuff
    //

    public static final String EXORCIZED_TREE_MODEL = "exorcizedTreeModel";

    public ExorcizedTreeModel(TreeModel ghostTreeModel, TaxonList corporealLineages) {
        this(EXORCIZED_TREE_MODEL, ghostTreeModel, corporealLineages);
    }

    public ExorcizedTreeModel(String name, Tree hauntedTree, TaxonList corporealLineages) {

        super(name, true);

        this.hauntedTree = hauntedTree;
        this.corporealLineages = corporealLineages.asList();

//        createCorporealTree(ExorcizedTreeModel, ghostLineages);

        if (hauntedTree instanceof TreeModel) {
            addModel((TreeModel) hauntedTree);
        }

        hauntedTreeChanged = true;
    }

    private void updateExorcizedTree() {
        updateExorcizedTree(hauntedTree, hauntedTree.getRoot());
        hauntedTreeChanged = false;
    }

    private void updateExorcizedTree(Tree hauntedTree, NodeRef hauntedNode) {

    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        assert model == hauntedTree;
        hauntedTreeChanged = true;
        fireModelChanged();
    }

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


//    /**
//     * @param ExorcizedTreeModel
//     * @param ghostLineages
//     */
//    public void createCorporealTree(TreeModel ExorcizedTreeModel, TaxonList ghostLineages) {
//
//        //start with setting the external node heights
//        for (int i = 0; i < this.getExternalNodeCount(); i++) {
//            this.setNodeHeight(this.getExternalNode(nodeMap[i]), nodeHeights[i]);
//        }
//        //set the internal node heights
//        for (int i = 0; i < (this.getExternalNodeCount() - 1); i++) {
//            //No just restart counting, will fix later on in the code by adding additionalTaxa variable
//            this.setNodeHeight(this.getInternalNode(i), nodeHeights[this.getExternalNodeCount() + i]);
//        }
//
//        int newRootIndex = -1;
//        //now add the parent-child links again to ALL the nodes
//        for (int i = 0; i < edges.length; i++) {
//            if (edges[i] != -1) {
//                //make distinction between external nodes and internal nodes
//                if (i < this.getExternalNodeCount()) {
//                    //external node
//                    this.addChild(this.getNode(edges[i]), this.getExternalNode(nodeMap[i]));
//                    System.out.println("external: " + edges[i] + " > " + nodeMap[i]);
//                } else {
//                    //internal node
//                    this.addChild(this.getNode(edges[i]), this.getNode(i));
//                    System.out.println("internal: " + edges[i] + " > " + i);
//                }
//            } else {
//                newRootIndex = i;
//            }
//        }
//
//        //not possible to determine correct ordering of child nodes in the loop where they're being assigned
//        //hence perform possible swaps in a separate loop
//
//        for (int i = 0; i < edges.length; i++) {
//                if (edges[i] != -1) {
//                    if(i < this.externalNodeCount) {
//                        if (childOrder[i] == 0 && nodes[edges[i]].getChild(0) != nodes[nodeMap[i]]) {
//                            //swap child nodes
//                            Node childOne = nodes[edges[i]].removeChild(0);
//                            Node childTwo = nodes[edges[i]].removeChild(1);
//                            nodes[edges[i]].addChild(childTwo);
//                            nodes[edges[i]].addChild(childOne);
//                        }
//                    }else{
//                        if (childOrder[i] == 0 && nodes[edges[i]].getChild(0) != nodes[i]) {
//                            //swap child nodes
//                            Node childOne = nodes[edges[i]].removeChild(0);
//                            Node childTwo = nodes[edges[i]].removeChild(1);
//                            nodes[edges[i]].addChild(childTwo);
//                            nodes[edges[i]].addChild(childOne);
//                        }
//                    }
//                }
//
//        }
//
//        this.setRoot(nodes[newRootIndex]);
//    }

//    NodeRef createCorporealNodes(NodeRef node) {
//        if (isExternal(node)) {
//            Taxon taxon = getNodeTaxon(node);
//            if (ghostLineages.contains(taxon)) {
//                return null;
//            } else {
//                ExorcizedTreeModel.addTaxon(taxon);
//            }
//        } else {
//            for (int i = 0; i < getChildCount(node); i++) {
//                createCorporealNodes(getChild(node, i));
//            }
//        }
//    }

    @Override
    public List<Citation> getCitations() {
        // @todo add a citation in here
        return Collections.EMPTY_LIST;
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
        throw new UnsupportedOperationException("Not supported yet.");
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
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getExternalNodeCount() {
        return corporealLineages.size();
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
        throw new UnsupportedOperationException("Not supported yet.");
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
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean isRoot(NodeRef node) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getChildCount(NodeRef node) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public NodeRef getChild(NodeRef node, int j) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public NodeRef getParent(NodeRef node) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Taxon getTaxon(int taxonIndex) {
        return corporealLineages.get(taxonIndex);
    }
}
