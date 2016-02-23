/*
 * PartitionedTreeModel.java
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

package dr.evomodel.epidemiology.casetocase;

import dr.app.beauti.util.TreeUtils;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.tree.TreeModelParser;
import dr.inference.model.*;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;

import java.util.*;

/**
 * TreeModel plus partition information
 *
 * todo a lot of methods should eventually move here
 */
public class PartitionedTreeModel extends TreeModel {

    private BranchMapModel branchMap;

    public final static String PARTITIONED_TREE_MODEL = "partitionedTreeModel";
    Set<NodeRef> partitionsQueue = new HashSet<NodeRef>();

    public PartitionedTreeModel(String id, Tree tree, BranchMapModel branchMap){
        super(id, tree);
        this.branchMap = branchMap;
    }

    public PartitionedTreeModel(String id, Tree tree){
        super(id, tree);
        branchMap = new BranchMapModel(this);
    }

    public PartitionedTreeModel(TreeModel treeModel, BranchMapModel branchMap){
        this(PARTITIONED_TREE_MODEL, treeModel, branchMap);
    }

    public PartitionedTreeModel(TreeModel treeModel){
        this(PARTITIONED_TREE_MODEL, treeModel);
    }

    public void partitionsChangingAlert(HashSet<AbstractCase> casesToRecalculate){
        // TreeLikelihood and TreeParameter listeners are irrelevant
        listenerHelper.fireModelChanged(this, new PartitionsChangedEvent(casesToRecalculate));
    }

    public void partitionChangingAlert(AbstractCase caseToRecalculate){
        HashSet<AbstractCase> out = new HashSet<AbstractCase>();
        out.add(caseToRecalculate);

        partitionsChangingAlert(out);
    }

    public void universalAlert(){
        HashSet<AbstractCase> allCases = new HashSet<AbstractCase>(Arrays.asList(branchMap.getArrayCopy()));
        partitionsChangingAlert(allCases);
    }

    public BranchMapModel getBranchMap(){
        return branchMap;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

        // shouldn't be any

    }

    public class PartitionsChangedEvent{

        private final HashSet<AbstractCase> casesToRecalculate;

        public PartitionsChangedEvent(HashSet<AbstractCase> casesToRecalculate){
            this.casesToRecalculate = casesToRecalculate;
        }

        public HashSet<AbstractCase> getCasesToRecalculate(){
            return casesToRecalculate;
        }
    }


    public void pushNodePartitionsChangedEvent(NodeRef node){
        int nodeNumber = node.getNumber();

        if(!inTreeEdit()){
            partitionChangingAlert(branchMap.get(nodeNumber));
        } else {
            partitionsQueue.add(node);
        }
    }


    public void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        super.handleVariableChangedEvent(variable, index, type);

        if (type == Parameter.ChangeType.ALL_VALUES_CHANGED) {
            //this signals events where values in all dimensions of a parameter is changed.
            universalAlert();
        } else {

            final NodeRef node = getNodeOfParameter((Parameter) variable);

            partitionsChangingAlert(adjacentPartitions(node));
        }
    }

    public HashSet<AbstractCase> adjacentPartitions(NodeRef node){
        HashSet<AbstractCase> changedCases = new HashSet<AbstractCase>();
        ArrayList<NodeRef> affectedNodes = new ArrayList<NodeRef>();

        affectedNodes.add(node);
        affectedNodes.add(getParent(node));
        affectedNodes.add(getChild(node, 0));
        affectedNodes.add(getChild(node, 1));

        for(NodeRef aNode : affectedNodes){
            if(aNode!=null){
                changedCases.add(branchMap.get(aNode.getNumber()));
            }
        }

        return changedCases;
    }

    private void flushQueue(){
        if(inTreeEdit()){
            throw new RuntimeException("Wait until you've finished editing the tree before flushing the partition" +
                    "queue");
        }

        for(NodeRef node : partitionsQueue){
            AbstractCase nodePartition = branchMap.get(node.getNumber());
            partitionChangingAlert(nodePartition);
            NodeRef parent = getParent(node);
            if(parent!=null && branchMap.get(node.getNumber())!=branchMap.get(parent.getNumber())){
                partitionChangingAlert(branchMap.get(parent.getNumber()));
            }
        }

        partitionsQueue.clear();

    }

    public void addChild(NodeRef p, NodeRef c) {

        pushNodePartitionsChangedEvent(p);
        pushNodePartitionsChangedEvent(c);

        super.addChild(p, c);
    }

    public void removeChild(NodeRef p, NodeRef c) {

        pushNodePartitionsChangedEvent(p);
        pushNodePartitionsChangedEvent(c);

        super.removeChild(p, c);

    }

    public void setNodeHeight(NodeRef n, double height) {

        partitionsChangingAlert(adjacentPartitions(n));

        super.setNodeHeight(n, height);
    }



    // anything you do to the partitions must be finished before you call this.

    public void endTreeEdit() {
        super.endTreeEdit();

        // todo in the end, want to check the tree partitions are sane here before flushing the queue

        flushQueue();
    }

    public boolean checkPartitions(){
        return checkPartitions(branchMap, true);
    }

    protected boolean checkPartitions(BranchMapModel map, boolean verbose){
        boolean foundProblem = false;
        for(int i=0; i<getInternalNodeCount(); i++){
            boolean foundTip = false;
            for(Integer nodeNumber : samePartitionElement(getInternalNode(i))){
                if(isExternal(getNode(nodeNumber))){
                    foundTip = true;
                }
            }
            if(!foundProblem && !foundTip){
                foundProblem = true;
                if(verbose){
                    System.out.println("Node "+(i+getExternalNodeCount()) + " is not connected to a tip");
                }
            }

        }

        // @todo wasteful - something accessible should keep a list of cases

        for(int i=0; i<getExternalNodeCount(); i++){
            AbstractCase aCase = branchMap.get(i);

            int[] tips = allTipsForThisCase(aCase);

            NodeRef tipMRCA = Tree.Utils.getCommonAncestor(this, tips);

            if(branchMap.get(tipMRCA.getNumber())!=aCase){
                throw new BadPartitionException("Node partition disconnected");
            }


        }


        return !foundProblem;
    }

    //Return a set of nodes that are not descendants of (or equal to) the current node and are in the same partition as
    // it. If flagForRecalc is true, then this also sets the flags for likelihood recalculation for all these nodes
    // to true


    public HashSet<Integer> samePartitionDownTree(NodeRef node){

        HashSet<Integer> out = new HashSet<Integer>();
        AbstractCase painting = branchMap.get(node.getNumber());
        NodeRef currentNode = node;
        NodeRef parentNode = getParent(node);
        while(parentNode!=null && branchMap.get(parentNode.getNumber())==painting){
            out.add(parentNode.getNumber());
            if(countChildrenInSamePartition(parentNode)==2){
                NodeRef otherChild = sibling(this, currentNode);
                out.add(otherChild.getNumber());
                out.addAll(samePartitionUpTree(otherChild));
            }
            currentNode = parentNode;
            parentNode = getParent(currentNode);
        }
        return out;
    }

    //Return a set of nodes that are descendants (and not equal to) the current node and are in the same partition as
    // it.



    public HashSet<Integer> samePartitionUpTree(NodeRef node){
        HashSet<Integer> out = new HashSet<Integer>();
        AbstractCase painting = branchMap.get(node.getNumber());
        for(int i=0; i< getChildCount(node); i++){
            if(branchMap.get(getChild(node,i).getNumber())==painting){
                out.add(getChild(node,i).getNumber());
                out.addAll(samePartitionUpTree(getChild(node, i)));
            }
        }
        return out;
    }


    public Integer[] samePartitionElement(NodeRef node){
        HashSet<Integer> out = new HashSet<Integer>();
        out.add(node.getNumber());
        out.addAll(samePartitionDownTree(node));
        out.addAll(samePartitionUpTree(node));
        return out.toArray(new Integer[out.size()]);
    }

    private int[] allTipsForThisCase(AbstractCase thisCase){
        ArrayList<Integer> listOfRefs = new ArrayList<Integer>();

        for(int i=0; i<getExternalNodeCount(); i++){
            if(branchMap.get(i)==thisCase){
                listOfRefs.add(i);
            }

        }

        int[] out = new int[listOfRefs.size()];

        for(int i=0; i<out.length; i++){out[i] = listOfRefs.get(i);}

        return out;

    }


    public NodeRef getEarliestNodeInPartition(AbstractCase thisCase){
        if(thisCase.wasEverInfected()) {

            int[] tips = allTipsForThisCase(thisCase);

            NodeRef tipMRCA = Tree.Utils.getCommonAncestor(this, tips);

            if(branchMap.get(tipMRCA.getNumber())!=thisCase){
                throw new BadPartitionException("Node partition disconnected");
            }

            NodeRef child = tipMRCA;
            NodeRef parent = getParent(child);
            boolean transmissionFound = false;
            while (!transmissionFound) {
                if (branchMap.get(child.getNumber()) != branchMap.get(parent.getNumber())) {
                    transmissionFound = true;
                } else {
                    child = parent;
                    parent = getParent(child);
                    if (parent == null) {
                        transmissionFound = true;
                    }
                }
            }
            return child;
        }
        return null;
    }

    public HashSet<AbstractCase> getDescendants(AbstractCase thisCase){
        HashSet<AbstractCase> out = new HashSet<AbstractCase>(getInfectees(thisCase));

        if(thisCase.wasEverInfected()) {
            for (AbstractCase child : out) {
                out.addAll(getDescendants(child));
            }
        }
        return out;
    }




    /* Return the case that infected this case */

    /* Return the case which was the infector in the infection event represented by this node */

    public AbstractCase getInfector(AbstractCase thisCase){
        if(thisCase.wasEverInfected()) {
            int[] tips = allTipsForThisCase(thisCase);

            NodeRef tipMRCA = Tree.Utils.getCommonAncestor(this, tips);

            if(branchMap.get(tipMRCA.getNumber())!=thisCase){
                throw new BadPartitionException("Node partition disconnected");
            }

            NodeRef currentNode = tipMRCA;

            while(branchMap.get(currentNode.getNumber())==thisCase){
                currentNode = getParent(currentNode);
                if(currentNode==null){
                    return null;
                }
            }
            return branchMap.get(currentNode.getNumber());


        }
        return null;
    }

    public AbstractCase getRootCase(){
        return branchMap.get(getRoot().getNumber());
    }



    public HashSet<AbstractCase> getInfectees(AbstractCase thisCase){
        if(thisCase.wasEverInfected()) {
            return getInfecteesInClade(getEarliestNodeInPartition(thisCase));
        }
        return new HashSet<AbstractCase>();
    }

    public HashSet<AbstractCase> getInfecteesInClade(NodeRef node){
        HashSet<AbstractCase> out = new HashSet<AbstractCase>();
        if(isExternal(node)){
            return out;
        } else {
            AbstractCase thisCase = branchMap.get(node.getNumber());
            for(int i=0; i<getChildCount(node); i++){
                NodeRef child = getChild(node, i);
                AbstractCase childCase = branchMap.get(child.getNumber());
                if(childCase!=thisCase){
                    out.add(childCase);
                } else {
                    out.addAll(getInfecteesInClade(child));
                }
            }
            return out;
        }
    }

    //infector of the case assigned to this node

    public AbstractCase getInfector(NodeRef node){
        if(isRoot(node) || node.getNumber() == getRoot().getNumber()){
            return null;
        } else {
            AbstractCase nodeCase = branchMap.get(node.getNumber());
            if(branchMap.get(getParent(node).getNumber())!=nodeCase){
                return branchMap.get(getParent(node).getNumber());
            } else {
                return getInfector(getParent(node));
            }
        }
    }


    /* Return the partition of the parent of this node */

    public AbstractCase getParentCase(NodeRef node){
        return branchMap.get(getParent(node).getNumber());
    }


    //Counts the children of the current node which are in the same partition element as itself



    public int countChildrenInSamePartition(NodeRef node){
        if(isExternal(node)){
            return -1;
        } else {
            int count = 0;
            AbstractCase parentCase = branchMap.get(node.getNumber());
            for(int i=0; i< getChildCount(node); i++){
                if(branchMap.get(getChild(node,i).getNumber())==parentCase){
                    count++;
                }
            }
            return count;
        }
    }



    public static NodeRef sibling(TreeModel tree, NodeRef node){
        if(tree.isRoot(node)){
            return null;
        } else {
            NodeRef parent = tree.getParent(node);
            for(int i=0; i<tree.getChildCount(parent); i++){
                if(tree.getChild(parent,i)!=node){
                    return tree.getChild(parent,i);
                }
            }
        }
        return null;
    }

}
