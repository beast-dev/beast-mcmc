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

package dr.evomodel.bigfasttree;

import dr.evolution.tree.*;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeModel;

import dr.util.Citation;


import java.util.*;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * A model component for trees that allows ghost lineages as a supertree of another TreeModel.
 *
 * @author JT McCrone
 * @author Andrew Rambaut
 */
public class GhostTreeModel extends BigFastTreeModel {

    //
    // Public stuff
    //

    public static final String GHOST_TREE_MODEL = "ghostTreeModel";

    public GhostTreeModel(Tree tree, TaxonList ghostLineages) {
        this(GHOST_TREE_MODEL, tree, ghostLineages);
    }

    public GhostTreeModel(String name, Tree tree, TaxonList ghostLineages) {

        super(name, tree, false, false);
        this.ghostLineages = ghostLineages.asList();
        this.corporealToGhostNodeMap = new int[2 * (super.getExternalNodeCount() - ghostLineages.getTaxonCount()) - 1];
        this.ghostToCorporealNodeMap = new int[super.getNodeCount()];
        this.ghostToNextOfKinMap = new int[super.getNodeCount()];

        for (int i = 0; i < super.getNodeCount(); i++) {
            this.ghostToCorporealNodeMap[i] = -1;
            this.ghostToNextOfKinMap[i] = -1;
        }
        createCorporealTree();
        this.updatedNodes = new boolean[super.getNodeCount()];
        Arrays.fill(updatedNodes, false);
    }

    /**
     * Set update flag for node in ghost tree.
     *
     * @param node - node in ghost tree
     */
    protected void updateNode(NodeRef node) {
        updatedNodes[node.getNumber()] = true;
        if (!super.isExternal(node)) {
            int corpNode = ghostToCorporealNodeMap[node.getNumber()];
            if (corpNode > -1) {
                availableNodeNumbers.add(corpNode);
                ghostToCorporealNodeMap[node.getNumber()] = -1;
                corporealToGhostNodeMap[corpNode] = -1;
            }
        }
        NodeRef parent = super.getParent(node);
        if (parent != null && !updatedNodes[parent.getNumber()]) {
            updateNode(parent);
        }
    }

    /**
     * Set update flag for all nodes
     */
    protected void updateAllNodes() {
        for (int i = 0; i < super.getNodeCount(); i++) {
            updatedNodes[i] = true;
        }
        availableNodeNumbers = IntStream.range(getExternalNodeCount(), getNodeCount())
                .boxed().collect(Collectors.toCollection(LinkedList::new));
    }


    /**
     * Return the ghost node that mirros a corporeal node
     *
     * @param corporealNode
     * @return ghostNode
     */
    private NodeRef getSpectralCounterPart(NodeRef corporealNode) {
        return super.getNode(corporealToGhostNodeMap[corporealNode.getNumber()]);
    }

    /**
     * get the corporeal node that mirros a ghost node
     *
     * @param ghostNode
     * @return corporealNode
     */
    private NodeRef getCorporealCounterPart(NodeRef ghostNode) {
        int corporealNodeIndex = ghostToCorporealNodeMap[ghostNode.getNumber()];
        if (corporealNodeIndex == -1) {
            return null;
        }
        return corporealTreeModel.getNode(corporealNodeIndex);
    }

    private NodeRef getNextCorporealDescendent(NodeRef ghostNode){
        if(hasCorporealCounterPart(ghostNode)){
            throw new IllegalArgumentException("Expected a node with degree 0 or 1 in corporeal tree but found degree 2");
        }
        int nodeNumber = ghostToNextOfKinMap[ghostNode.getNumber()];
        return nodeNumber==-1?null: super.getNode(nodeNumber);
    }

    /**
     * Helper function to check if a ghost node mirrors a corporeal node
     *
     * @param ghostNode
     * @return boolean
     */
    private boolean hasCorporealCounterPart(NodeRef ghostNode) {
        return ghostToCorporealNodeMap[ghostNode.getNumber()] != -1;
    }

    /**
     * Get the descendent ghost node(s) that mirror corporealNodes
     *
     * @param ghostNode
     * @return List of ghost nodes
     */
    //TODO cache
    private List<NodeRef> getDescendentsInCorporealTree(NodeRef ghostNode) {
        ArrayList<NodeRef> corporealNodes = new ArrayList<>();
        for (int i = 0; i < super.getChildCount(ghostNode); i++) {
            NodeRef child = super.getChild(ghostNode, i);
            if (super.isExternal(child)) {
                if (ghostToCorporealNodeMap[child.getNumber()] > -1) {
                    corporealNodes.add(child);
                }
            } else {
                List<NodeRef> childDescendents = getDescendentsInCorporealTree(child);
                if (childDescendents.size() == 2) {
                    corporealNodes.add(child);
                } else if (childDescendents.size() == 1) {
                    corporealNodes.addAll(childDescendents);
                }
            }
        }
        return corporealNodes;
    }
    // *****************************************************************
    // Interface MutableTree
    // *****************************************************************

    /**
     * Set a new node as root node.
     */
    public void setRoot(NodeRef newRoot) {
        super.setRoot(newRoot);
        updateNode(newRoot);
    }

    public void addChild(NodeRef p, NodeRef c) {
        super.addChild(p, c);
        updateNode(c);
    }

    public void removeChild(NodeRef p, NodeRef c) {
        updateNode(c);
        super.removeChild(p, c);

    }

    public boolean beginTreeEdit() {
        return super.beginTreeEdit();
    }

    public void endTreeEdit() {
        super.endTreeEdit();
        updateCorporealTreeModel();
    }


    public void setNodeHeight(NodeRef n, double height) {
        super.setNodeHeight(n, height);
        if (hasCorporealCounterPart(n)) {
            corporealTreeModel.adjustNodeHeight(getCorporealCounterPart(n), height);
        }
    }

    public void setNodeHeightQuietly(NodeRef n, double height) {
        super.setNodeHeightQuietly(n, height);
        if (hasCorporealCounterPart(n)) {
            corporealTreeModel.adjustNodeHeightQuietly(getSpectralCounterPart(n), height);
        }
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

//    public void makeDirty() {
//        updateCorporealTreeModelFromScratch();
//    }

    /**
     * Update the corporeal tree model topology
     */
    private void updateCorporealTreeModel() {
        //do the update;
        corporealTreeModel.beginTreeEdit();
        // Should return only 1 node
        List<NodeRef> relevantDescendents = updateCorporealTreeModel(super.getRoot());
        NodeRef corporealRoot;
        if(relevantDescendents.size()==2){
            //at the the root
            corporealRoot = getCorporealCounterPart(super.getRoot());
        }else if(relevantDescendents.size()==1){
            corporealRoot = getCorporealCounterPart(relevantDescendents.get(0));
        }else{
            throw new RuntimeException("Ghost tree root is out of sync. Root should be degree 1 or 2 in corporeal tree");
        }
        if (corporealRoot != corporealTreeModel.getRoot()) {
            corporealTreeModel.makeRoot(corporealRoot);
        }
        //TODO set root with output from above
        corporealTreeModel.endTreeEdit();
    }

    /**
     * Recursive version that traverses the ghost tree from the root and updates mappings between the trees and topology
     * in the corporeal tree.
     *
     * @param ghostNode
     */
    private List<NodeRef> updateCorporealTreeModel(NodeRef ghostNode) {
        ArrayList<NodeRef> relevantDescendents = new ArrayList<>();
        if (updatedNodes[ghostNode.getNumber()]) { //TODO update this

            for (int i = 0; i < super.getChildCount(ghostNode); i++) {
                NodeRef child = super.getChild(ghostNode, i);
                if (super.isExternal(child)) {
                    if (ghostToCorporealNodeMap[child.getNumber()] > -1) {
                        relevantDescendents.add(child);
                        //mark updated as false for child since wont visit?
                    }
                } else {
                    List<NodeRef> childDescendents = updateCorporealTreeModel(child);
                    if (childDescendents.size() == 2) {
                        relevantDescendents.add(child);
                    } else if (childDescendents.size() == 1) {
                        relevantDescendents.addAll(childDescendents);
                    }
                }
            }
            //vist me
            if(relevantDescendents.size()==0){
                ghostToNextOfKinMap[ghostNode.getNumber()]=-1;
            }else if(relevantDescendents.size()==1){
                ghostToNextOfKinMap[ghostNode.getNumber()]=relevantDescendents.get(0).getNumber();
            }else if (relevantDescendents.size() == 2) {
                Integer number = availableNodeNumbers.poll();
                updateNodeMaps(ghostNode.getNumber(), number);
                NodeRef me = corporealTreeModel.getNode(number);

                List<NodeRef> corporealDescendents = relevantDescendents
                        .stream()
                        .map(this::getCorporealCounterPart)
                        .collect(Collectors.toList());
                List<NodeRef> cachedCorporealDescendents = new ArrayList<>();
                for (int i = 0; i < corporealTreeModel.getChildCount(me); i++) {
                    cachedCorporealDescendents.add(corporealTreeModel.getChild(me, i));
                }
                for (NodeRef cachedDescendent :
                        cachedCorporealDescendents) {
                    if (!corporealDescendents.contains(cachedDescendent)) {
                        corporealTreeModel.disownChild(me, cachedDescendent);
                    }
                }
                for (NodeRef corporealDescendent : corporealDescendents) {
                    if (!cachedCorporealDescendents.contains(corporealDescendent)) {
                        NodeRef oldParent = corporealTreeModel.getParent(corporealDescendent);
                        if (oldParent != null) {
                            corporealTreeModel.disownChild(oldParent, corporealDescendent);
                        }
                        corporealTreeModel.adoptChild(me, corporealDescendent);
                    }
                }
                if(corporealTreeModel.getNodeHeight(me)!=super.getNodeHeight(ghostNode)){
                    corporealTreeModel.adjustNodeHeight(me, super.getNodeHeight(ghostNode));
                }
            }
        } else {
            if (hasCorporealCounterPart(ghostNode)) {
                relevantDescendents.add(ghostNode);
            } else {
                //either we are a degree 1 or degree 0 corpNode;
                NodeRef nextOfKin = getNextCorporealDescendent(ghostNode);
                if (nextOfKin != null) {
                    relevantDescendents.add(nextOfKin);
                }
            }
        }

        updatedNodes[ghostNode.getNumber()] = false;
        return relevantDescendents;

    }


    /**
     * Helper funciton to update node maps
     *
     * @param ghostNodeNumber
     * @param corpNodeNumber
     */
    private void updateNodeMaps(int ghostNodeNumber, int corpNodeNumber) {
        ghostToCorporealNodeMap[ghostNodeNumber] = corpNodeNumber;
        corporealToGhostNodeMap[corpNodeNumber] = ghostNodeNumber;
    }

    /**
     * Rebuild the corporeal tree from scratch
     */
    private void updateCorporealTreeModelFromScratch() {

        Map<NodeRef, NodeRef> temporaryMap = new HashMap<>();
        FlexibleNode corporealRootNode = createCorporealTree(this.getRoot(), temporaryMap);
        FlexibleTree flexibleTree = new FlexibleTree(corporealRootNode);
        flexibleTree.adoptTreeModelOrdering();
        corporealTreeModel.copyTopology(flexibleTree);

        for (Map.Entry<NodeRef, NodeRef> entry :
                temporaryMap.entrySet()) {
            int ghostNodeInt = entry.getKey().getNumber();
            int corporealNodeInt = entry.getValue().getNumber();
            ghostToCorporealNodeMap[ghostNodeInt] = corporealNodeInt;
            corporealToGhostNodeMap[corporealNodeInt] = ghostNodeInt;
        }
    }

    /**
     * Create the corporeal tree model and node maps
     */
    private void createCorporealTree() {
        //Traverse the tree and build the flexibleTree
        Map<NodeRef, NodeRef> temporaryMap = new HashMap<>();
        FlexibleNode corporealRootNode = createCorporealTree(this.getRoot(), temporaryMap);
        FlexibleTree flexibleTree = new FlexibleTree(corporealRootNode);
        flexibleTree.adoptTreeModelOrdering();
        corporealTreeModel = new CorporealTreeModel(this.getModelName() + "CorporealTree", flexibleTree);
        for (Map.Entry<NodeRef, NodeRef> entry :
                temporaryMap.entrySet()) {
            int ghostNodeInt = entry.getKey().getNumber();
            int corporealNodeInt = entry.getValue().getNumber();
            ghostToCorporealNodeMap[ghostNodeInt] = corporealNodeInt;
            corporealToGhostNodeMap[corporealNodeInt] = ghostNodeInt;
        }


        for (int j : corporealToGhostNodeMap) {
            NodeRef ghostNode = super.getNode(j);
            NodeRef ghostParent = super.getParent(ghostNode);
            while (ghostParent != null && !hasCorporealCounterPart(ghostParent)) {
                ghostToNextOfKinMap[ghostParent.getNumber()] = j;
                ghostNode = ghostParent;
                ghostParent = super.getParent(ghostNode);
            }
        }

    }

    /**
     * Recursive helper function
     *
     * @param nodeRef
     * @param nodeMap
     * @return
     */
    private FlexibleNode createCorporealTree(NodeRef nodeRef, Map<NodeRef, NodeRef> nodeMap) {
        if (this.isExternal(nodeRef)) {
            if (this.ghostLineages.contains(this.getNodeTaxon(nodeRef))) {
                return null;
            } else {
                FlexibleNode corporealNode = new FlexibleNode(this.getNodeTaxon(nodeRef));
                corporealNode.setHeight(this.getNodeHeight(nodeRef));
                nodeMap.put(nodeRef, corporealNode);
                return corporealNode;
            }
        } else {
            List<FlexibleNode> children = new ArrayList<>();
            for (int i = 0; i < this.getChildCount(nodeRef); i++) {
                FlexibleNode child = createCorporealTree(this.getChild(nodeRef, i), nodeMap);
                if (child != null) {
                    children.add(child);
                }
            }
            if (children.size() == 0) {
                return null;
            } else if (children.size() == 1) {
                return children.get(0);
            } else {
                FlexibleNode corporealNode = new FlexibleNode();
                corporealNode.setHeight(this.getNodeHeight(nodeRef));
                for (FlexibleNode child : children) {
                    child.setParent(corporealNode);
                    corporealNode.addChild(child);
                }
                nodeMap.put(nodeRef, corporealNode);
                return corporealNode;
            }

        }
    }

    public TreeModel getCorporealTreeModel() {
        return corporealTreeModel;
    }

    private final List<Taxon> ghostLineages;
    private CorporealTreeModel corporealTreeModel;
    private final boolean[] updatedNodes;
    private boolean seenRoot = false;
    Queue<Integer> availableNodeNumbers = new LinkedList<>();
    private final int[] corporealToGhostNodeMap;
    private final int[] ghostToCorporealNodeMap;
    private final int[] ghostToNextOfKinMap;

    @Override

    public List<Citation> getCitations() {
        // @todo add a citation in here
        return Collections.EMPTY_LIST;
    }

}
