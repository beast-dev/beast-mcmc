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
import java.util.stream.Stream;

/**
 * A model component for trees that allows ghost lineages as a supertree of another TreeModel.
 *
 * @author Andrew Rambaut
 * @author JT McCrone
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


        for (int i = 0; i < super.getNodeCount(); i++) {
            this.ghostToCorporealNodeMap[i] = -1;
        }
        createCorporealTree();
        this.updatedCorporealNodes= new boolean[corporealTreeModel.getNodeCount()];
        Arrays.fill(updatedCorporealNodes, false);

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

    private List<NodeRef> getDescendentsInCorporealTree(NodeRef n) {
        ArrayList<NodeRef> corporealNodes = new ArrayList<>();
        for (int i = 0; i < super.getChildCount(n); i++) {
            NodeRef child = getChild(n, i);
            if (isExternal(child)) {
                if (!ghostLineages.contains(super.getNodeTaxon(child))) {
                    corporealNodes.add(child);
                }
            }else{
                List<NodeRef> childDescendents = getDescendentsInCorporealTree(child);
                if(childDescendents.size()==2){
                    corporealNodes.add(child);
                }else if (childDescendents.size()==1){
                    corporealNodes.addAll(childDescendents);
                }
            }
        }
        return corporealNodes;
    }


    private NodeRef getSpectralCounterPart(NodeRef n) {
        return super.getNode(corporealToGhostNodeMap[n.getNumber()]);
    }

    private NodeRef getCorporealCounterPart(NodeRef n) {
        int corporealNodeIndex = ghostToCorporealNodeMap[n.getNumber()];
        if(corporealNodeIndex==-1){
            return null;
        }
        return corporealTreeModel.getNode(corporealNodeIndex);
    }

    private boolean hasCorporealCounterPart(NodeRef n) {
        return ghostToCorporealNodeMap[n.getNumber()] != -1;
    }


    public void removeChild(NodeRef p, NodeRef c) {
        super.removeChild(p, c);

    }

    public boolean beginTreeEdit() {
        return super.beginTreeEdit();
    }

    public void endTreeEdit() {
        // and cleanup
        super.endTreeEdit();
        corporealTreeModel.beginTreeEdit();
        updateCorporealTreeModel();
        corporealTreeModel.endTreeEdit();
    }


    public void setNodeHeight(NodeRef n, double height) {
        super.setNodeHeight(n, height);
        if (hasCorporealCounterPart(n)) {
            corporealTreeModel.setNodeHeight(getCorporealCounterPart(n), height);
        }
    }

    public void setNodeHeightQuietly(NodeRef n, double height) {
        super.setNodeHeightQuietly(n, height);
        if (hasCorporealCounterPart(n)) {
            corporealTreeModel.setNodeHeightQuietly(getSpectralCounterPart(n), height);
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

    public int getCorporealDegree(NodeRef nodeRef) {
       return getDescendentsInCorporealTree(nodeRef).size();
    }

    public void makeDirty(){
        updateCorporealTreeModelFromScratch();
    }

    private void updateCorporealTreeModel(){
        vistNode(super.getRoot());
    }

    private void vistNode(NodeRef nodeRef) {
        int corporealDegree = getCorporealDegree(nodeRef);
        boolean inCorporealTree = hasCorporealCounterPart(nodeRef);



        if (corporealDegree == 2) {
            List<NodeRef> corporealDescendents = getDescendentsInCorporealTree(nodeRef)
                    .stream()
                    .map(this::getCorporealCounterPart)
                    .collect(Collectors.toList());
            if (inCorporealTree) {
                NodeRef corporealNode = getCorporealCounterPart(nodeRef);
                List<NodeRef> cachedCorporealDescendents = new ArrayList<>();
                for (int i = 0; i < corporealTreeModel.getChildCount(corporealNode); i++) {
                    cachedCorporealDescendents.add(corporealTreeModel.getChild(corporealNode, i));
                }
                boolean updated = false;
                for (NodeRef cachedDescendent :
                        cachedCorporealDescendents) {
                    if(!corporealDescendents.contains(cachedDescendent)){
                        corporealTreeModel.removeChild(corporealNode,cachedDescendent);
                        updated =true;
                    }
                }
                for (NodeRef corporealDescendent :corporealDescendents) {
                    if(!cachedCorporealDescendents.contains(corporealDescendent)){

                        NodeRef oldParent = corporealTreeModel.getParent(corporealDescendent);
                        if(oldParent!=null){
                            corporealTreeModel.removeChild(oldParent,corporealDescendent);

                        }
                        corporealTreeModel.addChild(corporealNode,corporealDescendent);

                        if(corporealDescendent==corporealTreeModel.getRoot()){
                            corporealTreeModel.setRoot(corporealNode);
                        }
                        updated =true;
                    }
                }
                if(updated){
                    corporealTreeModel.setNodeHeight(corporealNode,super.getNodeHeight(nodeRef));
                }
            }else{
                boolean isNewNodeRoot = false;
// presumably we "looked past" a degree 1 corporeal node when we looked for the children
                    List<NodeRef> cachedDegree1Parents = corporealDescendents
                            .stream()
                            .map(n->corporealTreeModel.getParent(n))
                            .filter(Objects::nonNull)
                            .map(this::getSpectralCounterPart)
                            .filter(n->getCorporealDegree(n)==1)
                            .map(this::getCorporealCounterPart)
                            .collect(Collectors.toList());
                    assert cachedDegree1Parents.size()==1;

                NodeRef corporealNode =  cachedDegree1Parents.get(0);
                corporealTreeModel.setNodeHeight(corporealNode,super.getNodeHeight(nodeRef));
                updateNodeMaps(nodeRef,corporealNode);

                List<NodeRef> cachedCorporealDescendents = new ArrayList<>();
                for (int i = 0; i < corporealTreeModel.getChildCount(corporealNode); i++) {
                    cachedCorporealDescendents.add(corporealTreeModel.getChild(corporealNode, i));
                }
                for (NodeRef cachedDescendent :
                        cachedCorporealDescendents) {
                    if(!corporealDescendents.contains(cachedDescendent)){
                        corporealTreeModel.removeChild(corporealNode,cachedDescendent);
                    }
                }
                for (NodeRef corporealDescendent :corporealDescendents) {
                    if(!cachedCorporealDescendents.contains(corporealDescendent)){

                        NodeRef oldParent = corporealTreeModel.getParent(corporealDescendent);
                        if(oldParent!=null){
                            corporealTreeModel.removeChild(oldParent,corporealDescendent);

                        }
                        corporealTreeModel.addChild(corporealNode,corporealDescendent);

                        if(corporealDescendent==corporealTreeModel.getRoot()){
                            corporealTreeModel.setRoot(corporealNode);
                        }
                    }
                }


            }
        }else if(corporealDegree==1){
            if(inCorporealTree){
                NodeRef corpNode = getCorporealCounterPart(nodeRef);
                if(corpNode==corporealTreeModel.getRoot()){
                    //move up
                    NodeRef lineageToCorpRoot = getDescendentsInCorporealTree(nodeRef).get(0);
                    corporealTreeModel.setNodeHeight(corpNode,super.getNodeHeight(lineageToCorpRoot));
                    updateNodeMaps(lineageToCorpRoot,corpNode);
                }else { //free for use later
                    NodeRef corpChild = getCorporealCounterPart(getDescendentsInCorporealTree(nodeRef).get(0));
                    NodeRef corpParent = corporealTreeModel.getParent(corpNode);

                    corporealTreeModel.removeChild(corpParent, corpNode);
                    corporealTreeModel.removeChild(corpNode, corpChild);
                    corporealTreeModel.addChild(corpParent, corpChild);
                }
            }
        }else if(corporealDegree==0){
            if(inCorporealTree){
                NodeRef corpNode = getCorporealCounterPart(nodeRef);
                NodeRef corpP = corporealTreeModel.getParent(corpNode);
                corporealTreeModel.removeChild(corpP, corpNode);

            }
        }

        for (int i = 0; i < getChildCount(nodeRef); i++) {
            NodeRef child = getChild(nodeRef, i);
            if(!isExternal(child)){
                vistNode(child);
            }
        }

    }

    private void updateNodeMaps(NodeRef ghostNode,NodeRef corpNode){
        ghostToCorporealNodeMap[corporealToGhostNodeMap[corpNode.getNumber()]]=-1;
        ghostToCorporealNodeMap[ghostNode.getNumber()] = corpNode.getNumber();
        corporealToGhostNodeMap[corpNode.getNumber()] = ghostNode.getNumber();
    }

    private NodeRef getSibling(Tree tree,NodeRef node){
        NodeRef parent = tree.getParent(node);
        if(tree.getChild(parent,0)==node){
            return tree.getChild(parent,1);
        }
        return tree.getChild(parent, 0);
    }


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


    private void createCorporealTree() {
        //Traverse the tree and build the flexibleTree
        Map<NodeRef, NodeRef> temporaryMap = new HashMap<>();
        FlexibleNode corporealRootNode = createCorporealTree(this.getRoot(), temporaryMap);
        FlexibleTree flexibleTree = new FlexibleTree(corporealRootNode);
        flexibleTree.adoptTreeModelOrdering();
        corporealTreeModel = new BigFastTreeModel(this.getModelName() + "CorporealTree", flexibleTree);
        for (Map.Entry<NodeRef, NodeRef> entry :
                temporaryMap.entrySet()) {
            int ghostNodeInt = entry.getKey().getNumber();
            int corporealNodeInt = entry.getValue().getNumber();
            ghostToCorporealNodeMap[ghostNodeInt] = corporealNodeInt;
            corporealToGhostNodeMap[corporealNodeInt] = ghostNodeInt;
        }
    }

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
    private BigFastTreeModel corporealTreeModel;
    private boolean[] updatedCorporealNodes;
    private final int[] corporealToGhostNodeMap;
    private final int[] ghostToCorporealNodeMap;


    @Override

    public List<Citation> getCitations() {
        // @todo add a citation in here
        return Collections.EMPTY_LIST;
    }

}
