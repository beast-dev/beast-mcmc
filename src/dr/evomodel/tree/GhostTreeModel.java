/*
 * TreeModel.java
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

package dr.evomodel.tree;

import dr.evolution.tree.*;
import dr.evolution.util.MutableTaxonListListener;
import dr.evolution.util.Taxon;
import dr.inference.model.*;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.Keywordable;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.*;

/**
 * A model component for trees that allows ghost lineages as a supertree of another TreeModel.
 *
 * @author Andrew Rambaut
 */
public class GhostTreeModel extends TreeModel {

    //
    // Public stuff
    //

    public static final String GHOST_TREE_MODEL = "ghostTreeModel";

    public GhostTreeModel(Tree tree, TreeModel corporealTreeModel) {
        this(TREE_MODEL, tree, corporealTreeModel);
    }

    public GhostTreeModel(String name, Tree tree, TreeModel corporealTreeModel) {

        super(name, tree, false, false, false);

        this.corporealTreeModel = corporealTreeModel;
    }

    /**
     * Called when a parameter changes.
     */
    public void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        final Node node = getNodeOfParameter((Parameter) variable);
        if (type == Parameter.ChangeType.ALL_VALUES_CHANGED) {
            //this signals events where values in all dimensions of a parameter is changed.
            pushTreeChangedEvent(new TreeChangedEvent(node, (Parameter) variable, TreeChangedEvent.CHANGE_IN_ALL_INTERNAL_NODES));
        } else {
            pushTreeChangedEvent(node, (Parameter) variable, index);
        }
    }

    // *****************************************************************
    // Interface MutableTree
    // *****************************************************************

    /**
     * Set a new node as root node.
     */
    public void setRoot(NodeRef newRoot) {
        super.setRoot(newRoot);
    }

    public void addChild(NodeRef p, NodeRef c) {
super.addChild(p, c);
    }

    public void removeChild(NodeRef p, NodeRef c) {
super.removeChild(p, c);
    }

    public boolean beginTreeEdit() {
        subTreeModel.beginTreeEdit();
        return super.beginTreeEdit();
    }

    public void endTreeEdit() {
        subTreeModel.endTreeEdit();
        super.endTreeEdit();
    }

    public void setNodeHeight(NodeRef n, double height) {
        super.setNodeHeight(n, height);
    }

    public void setNodeHeightQuietly(NodeRef n, double height) {
        super.setNodeHeightQuietly(n, height);
    }

    public void setNodeRate(NodeRef n, double rate) {
        throw new UnsupportedOperationException("Function not available in GhostTreeModel");

    }

    public void setNodeTrait(NodeRef n, String name, double value) {
        throw new UnsupportedOperationException("Function not available in GhostTreeModel");
    }

    public void setMultivariateTrait(NodeRef n, String name, double[] value) {
        throw new UnsupportedOperationException("Function not available in GhostTreeModel");
    }

    /**
     * Copies the node connections from this TreeModel's nodes array to the
     * destination array. Basically it connects up the nodes in destination
     * in the same way as this TreeModel is set up. This method is package
     * private.
     */
    private void copyNodeStructure(Node[] destination) {

        if (nodes.length != destination.length) {
            throw new IllegalArgumentException("Node arrays are of different lengths");
        }

        for (int i = 0, n = nodes.length; i < n; i++) {
            Node node0 = nodes[i];
            Node node1 = destination[i];

            // the parameter values are automatically stored and restored
            // just need to keep the links
            node1.heightParameter = node0.heightParameter;
            node1.rateParameter = node0.rateParameter;
            node1.traitParameters = node0.traitParameters;

            if (node0.parent != null) {
                node1.parent = storedNodes[node0.parent.getNumber()];
            } else {
                node1.parent = null;
            }

            if (node0.leftChild != null) {
                node1.leftChild = storedNodes[node0.leftChild.getNumber()];
            } else {
                node1.leftChild = null;
            }

            if (node0.rightChild != null) {
                node1.rightChild = storedNodes[node0.rightChild.getNumber()];
            } else {
                node1.rightChild = null;
            }
        }
    }

    /**
     * Copies a different tree into the current treeModel. Needs to reconnect
     * the existing internal and external nodes, taking into account that the
     * node numbers of the external nodes may differ between the two trees.
     */
    public void adoptTreeStructure(Tree donor) {

        /*System.err.println("internalNodeCount: " + this.internalNodeCount);
          System.err.println("externalNodeCount: " + this.externalNodeCount);
          for (int i = 0; i < this.nodeCount; i++) {
              System.err.println(nodes[i]);
          }*/

        //first remove all the child nodes of the internal nodes
        for (int i = this.externalNodeCount; i < this.nodeCount; i++) {
            int childCount = nodes[i].getChildCount();
            for (int j = 0; j < childCount; j++) {
                nodes[i].removeChild(j);
            }
        }

        // set-up nodes in this.nodes[] to mirror connectedness in donor via a simple recursion on donor.getRoot()
        addNodeStructure(donor, donor.getRoot());

        //Tree donor has no rates nor traits, only heights

    }

    /**
     * Modifies the current tree by adopting the provided collection of edges
     * @param edges Edges are provided as index: child number; parent: array entry
     * @param nodeHeights Also sets the node heights to the provided values
     * @param childOrder Array that contains whether a child node is left or right child
     */
    public void adoptTreeStructure(int[] edges, double[] nodeHeights, int[] childOrder, String[] taxaNames) {

        int[] nodeMap = createNodeMap(taxaNames);

        if (this.nodeCount != edges.length) {
            throw new RuntimeException("Incorrect number of edges provided: " + edges.length + " versus " + this.nodeCount + " nodes.");
        }

        //first remove all the child nodes of the internal nodes
        for (int i = this.externalNodeCount; i < this.nodeCount; i++) {
            int childCount = nodes[i].getChildCount();
            for (int j = 0; j < childCount; j++) {
                nodes[i].removeChild(j);
            }
        }

        //start with setting the external node heights
        for (int i = 0; i < this.getExternalNodeCount(); i++) {
            this.setNodeHeight(this.getExternalNode(nodeMap[i]), nodeHeights[i]);
        }
        //set the internal node heights
        for (int i = 0; i < (this.getExternalNodeCount() - 1); i++) {
            //No just restart counting, will fix later on in the code by adding additionalTaxa variable
            this.setNodeHeight(this.getInternalNode(i), nodeHeights[this.getExternalNodeCount() + i]);
        }

        int newRootIndex = -1;
        //now add the parent-child links again to ALL the nodes
        for (int i = 0; i < edges.length; i++) {
            if (edges[i] != -1) {
                //make distinction between external nodes and internal nodes
                if (i < this.getExternalNodeCount()) {
                    //external node
                    this.addChild(this.getNode(edges[i]), this.getExternalNode(nodeMap[i]));
                    System.out.println("external: " + edges[i] + " > " + nodeMap[i]);
                } else {
                    //internal node
                    this.addChild(this.getNode(edges[i]), this.getNode(i));
                    System.out.println("internal: " + edges[i] + " > " + i);
                }
            } else {
                newRootIndex = i;
            }
        }

        //not possible to determine correct ordering of child nodes in the loop where they're being assigned
        //hence perform possible swaps in a separate loop

        for (int i = 0; i < edges.length; i++) {
                if (edges[i] != -1) {
                    if(i < this.externalNodeCount) {
                        if (childOrder[i] == 0 && nodes[edges[i]].getChild(0) != nodes[nodeMap[i]]) {
                            //swap child nodes
                            Node childOne = nodes[edges[i]].removeChild(0);
                            Node childTwo = nodes[edges[i]].removeChild(1);
                            nodes[edges[i]].addChild(childTwo);
                            nodes[edges[i]].addChild(childOne);
                        }
                    }else{
                        if (childOrder[i] == 0 && nodes[edges[i]].getChild(0) != nodes[i]) {
                            //swap child nodes
                            Node childOne = nodes[edges[i]].removeChild(0);
                            Node childTwo = nodes[edges[i]].removeChild(1);
                            nodes[edges[i]].addChild(childTwo);
                            nodes[edges[i]].addChild(childOne);
                        }
                    }
                }

        }

        this.setRoot(nodes[newRootIndex]);
    }

    public TreeModel getCorporealTreeModel() {
        return corporealTreeModel;
    }

    private final TreeModel corporealTreeModel;

    @Override
    public List<Citation> getCitations() {
        // @todo add a citation in here
        return  Collections.EMPTY_LIST;
    }

}
