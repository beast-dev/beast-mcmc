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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.math.MathUtils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * TreeModel plus partition information
 *
 * todo a lot of methods should eventually move here
 */
public class PartitionedTreeModel extends TreeModel {

    private final AbstractOutbreak outbreak;
    private BranchMapModel branchMap;
    private final int elementCount;

    public final static String PARTITIONED_TREE_MODEL = "partitionedTreeModel";
    Set<NodeRef> partitionsQueue = new HashSet<NodeRef>();

    public PartitionedTreeModel(String id, Tree tree, AbstractOutbreak outbreak){
        super(id, tree);
        this.outbreak = outbreak;
        elementCount = outbreak.infectedSize();
        branchMap = new BranchMapModel(this);
        partitionAccordingToRandomTT(false);
    }

    public PartitionedTreeModel(String id, Tree tree, AbstractOutbreak outbreak, String startingTTFileName){
        super(id, tree);
        this.outbreak = outbreak;
        elementCount = outbreak.infectedSize();
        branchMap = new BranchMapModel(this);
        partitionAccordingToSpecificTT(startingTTFileName);
    }

    public PartitionedTreeModel(TreeModel treeModel, AbstractOutbreak outbreak){
        this(PARTITIONED_TREE_MODEL, treeModel, outbreak);
    }

    public PartitionedTreeModel(TreeModel treeModel, AbstractOutbreak outbreak, String startingTTFileName){
        this(PARTITIONED_TREE_MODEL, treeModel, outbreak, startingTTFileName);
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

            partitionsChangingAlert(adjacentElements(node));
        }
    }

    public HashSet<AbstractCase> adjacentElements(NodeRef node){
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
            AbstractCase nodeElement = branchMap.get(node.getNumber());
            partitionChangingAlert(nodeElement);
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

        partitionsChangingAlert(adjacentElements(n));

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


            NodeRef tipMRCA = caseMRCA(aCase);

            if(branchMap.get(tipMRCA.getNumber())!=aCase){
                throw new BadPartitionException("Node partition disconnected");
            }


        }


        return !foundProblem;
    }

    //Return a set of nodes that are not descendants of (or equal to) the current node and are in the same partition as
    // it.


    public HashSet<Integer> samePartitionElementUpTree(NodeRef node){
        HashSet<Integer> out = new HashSet<Integer>();
        AbstractCase elementCase = branchMap.get(node.getNumber());
        NodeRef currentNode = node;
        NodeRef parentNode = getParent(node);
        while(parentNode!=null && branchMap.get(parentNode.getNumber())==elementCase){
            out.add(parentNode.getNumber());
            if(countChildrenInSameElement(parentNode)==2){
                NodeRef otherChild = sibling(this, currentNode);
                out.add(otherChild.getNumber());
                out.addAll(samePartitionElementDownTree(otherChild));
            }
            currentNode = parentNode;
            parentNode = getParent(currentNode);
        }
        return out;
    }

    //Return a set of nodes that are descendants (and not equal to) the current node and are in the same partition as
    // it.

    public HashSet<Integer> samePartitionElementDownTree(NodeRef node){
        HashSet<Integer> out = new HashSet<Integer>();
        AbstractCase elementCase = branchMap.get(node.getNumber());
        for(int i=0; i< getChildCount(node); i++){
            if(branchMap.get(getChild(node,i).getNumber())==elementCase){
                out.add(getChild(node,i).getNumber());
                out.addAll(samePartitionElementDownTree(getChild(node, i)));
            }
        }
        return out;
    }


    public Integer[] samePartitionElement(NodeRef node){
        HashSet<Integer> out = new HashSet<Integer>();
        out.add(node.getNumber());
        out.addAll(samePartitionElementUpTree(node));
        out.addAll(samePartitionElementDownTree(node));
        return out.toArray(new Integer[out.size()]);
    }

    public int[] allTipsForThisCase(AbstractCase thisCase){
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


    public NodeRef getEarliestNodeInElement(AbstractCase thisCase){
        if(thisCase.wasEverInfected()) {

            NodeRef tipMRCA = caseMRCA(thisCase);

            if(branchMap.get(tipMRCA.getNumber())!=thisCase){
                throw new BadPartitionException("Node partition element disconnected");
            }

            NodeRef child = tipMRCA;
            NodeRef parent = getParent(child);
            boolean transmissionFound = parent == null;
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

    public AbstractCase getInfector(AbstractCase thisCase){
        if(thisCase.wasEverInfected()) {

            NodeRef tipMRCA = caseMRCA(thisCase);

            if(branchMap.get(tipMRCA.getNumber())!=thisCase){
                throw new BadPartitionException("Node partition element disconnected");
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
            return getInfecteesInClade(getEarliestNodeInElement(thisCase));
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


    /* Return the partition element of the parent of this node */

    public AbstractCase getParentCase(NodeRef node){
        return branchMap.get(getParent(node).getNumber());
    }


    public int getElementCount(){
        return elementCount;
    }


    //Counts the children of the current node which are in the same partition element as itself



    public int countChildrenInSameElement(NodeRef node){
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


    public NodeRef caseMRCA(AbstractCase aCase, boolean checkConnectedness){
        int[] caseTips = allTipsForThisCase(aCase);
        NodeRef mrca =  TreeUtils.getCommonAncestor(this, caseTips);

        if(checkConnectedness) {
            if (branchMap.get(mrca.getNumber()) != aCase) {
                throw new BadPartitionException("A partition element is disconnected");
            }
        }

        return mrca;
    }

    public NodeRef caseMRCA(AbstractCase aCase){
        return caseMRCA(aCase, true);
    }

    private HashSet<NodeRef> getDescendantTips(NodeRef node){
        HashSet<NodeRef> out = new HashSet<NodeRef>();
        if(isExternal(node)){
            out.add(node);
            return out;
        } else {
            out.addAll(getDescendantTips(getChild(node, 0)));
            out.addAll(getDescendantTips(getChild(node, 1)));
        }
        return out;
    }

    public boolean isAncestral(NodeRef node){
        AbstractCase currentCase = branchMap.get(node.getNumber());

        for(NodeRef tip : getDescendantTips(node)){
            if(branchMap.get(tip.getNumber())==currentCase){
                return true;
            }
        }

        return false;
    }

    public boolean isRootBlockedBy(AbstractCase aCase, AbstractCase potentialBlocker){
        return directDescendant(caseMRCA(aCase), caseMRCA(potentialBlocker));
    }

    public boolean isRootBlocked(AbstractCase aCase){
        for(AbstractCase anotherCase : outbreak.getCases()){
            if(anotherCase.wasEverInfected && anotherCase!=aCase){
                if(isRootBlockedBy(aCase, anotherCase)){
                    return true;
                }
            }
        }
        return false;
    }

    private HashSet<NodeRef> getTipsInThisPartitionElement(AbstractCase aCase){
        HashSet<NodeRef> out = new HashSet<NodeRef>();
        // todo check that external nodes come first

        for(int i=0; i<getExternalNodeCount(); i++){
            if(branchMap.get(i)==aCase){
                out.add(getExternalNode(i));
            }
        }

        return out;
    }

    private boolean directDescendant(NodeRef node, NodeRef possibleAncestor){
        NodeRef currentNode = node;

        while(currentNode!=null){
            if(currentNode==possibleAncestor){
                return true;
            }
            currentNode = getParent(currentNode);
        }
        return false;
    }

    private boolean directRelationship(NodeRef node1, NodeRef node2){
        return directDescendant(node1, node2) || directDescendant(node2, node1);
    }

       /* Populates the branch map for external nodes */

    private AbstractCase[] prepareExternalNodeMap(AbstractCase[] map){
        for(int i=0; i< getExternalNodeCount(); i++){
            TreeModel.Node currentExternalNode = (TreeModel.Node) getExternalNode(i);
            Taxon currentTaxon = currentExternalNode.taxon;
            for(AbstractCase thisCase : outbreak.getCases()){
                if(thisCase.wasEverInfected()) {
                    for (Taxon caseTaxon : thisCase.getAssociatedTaxa()) {
                        if (caseTaxon.equals(currentTaxon)) {
                            map[currentExternalNode.getNumber()] = thisCase;
                        }
                    }
                }
            }
        }
        return map;
    }

/*  The CSV file should have a header, and then lines matching each case to its infector*/

    private void partitionAccordingToSpecificTT(String networkFileName){
        System.out.println("Using specified starting transmission tree.");
        try{
            BufferedReader reader = new BufferedReader (new FileReader(networkFileName));
            HashMap<AbstractCase, AbstractCase> specificParentMap = new HashMap<AbstractCase, AbstractCase>();
            // skip header line
            reader.readLine();
            String currentLine = reader.readLine();
            while(currentLine!=null){
                currentLine = currentLine.replace("\"", "");
                String[] splitLine = currentLine.split("\\,");
                if(!splitLine[1].equals("Start")){
                    specificParentMap.put(outbreak.getCase(splitLine[0]), outbreak.getCase(splitLine[1]));
                } else {
                    specificParentMap.put(outbreak.getCase(splitLine[0]), null);
                }
                currentLine = reader.readLine();
            }
            reader.close();
            partitionAccordingToSpecificTT(specificParentMap);
        } catch(IOException e){
            throw new RuntimeException("Cannot read file: " + networkFileName );
        }
    }


    private void partitionAccordingToSpecificTT(HashMap<AbstractCase, AbstractCase> map){
        branchMap.setAll(prepareExternalNodeMap(new AbstractCase[getNodeCount()]), true);

        //various sanity checks

        for(AbstractCase aCase : map.keySet()){
            if(!aCase.wasEverInfected){
                throw new RuntimeException("This starting transmission tree involves never-infected cases");
            }
        }

        AbstractCase firstCase=null;
        int indexCaseCount = 0;

        for(AbstractCase aCase : outbreak.getCases()){
            if(aCase.wasEverInfected()) {
                if (map.get(aCase) == null) {
                    firstCase = aCase;
                    indexCaseCount++;
                }
            }
        }
        if(indexCaseCount==0){
            throw new RuntimeException("Given starting transmission tree appears to have a cycle");
        }
        if(indexCaseCount>1){
            throw new RuntimeException("Given starting transmission tree appears not to be connected");
        }


        NodeRef root = getRoot();
        specificallyPartitionDownwards(root, firstCase, map);
        if(!checkPartitions()){
            throw new RuntimeException("Given starting transmission tree is not compatible with the starting tree");
        }

    }

    private void specificallyPartitionDownwards(NodeRef node, AbstractCase thisCase,
                                                HashMap<AbstractCase, AbstractCase> map){
        if(isExternal(node)){
            return;
        }
        branchMap.set(node.getNumber(), thisCase, true);
        if(isAncestral(node)){
            for(int i=0; i<getChildCount(node); i++){
                specificallyPartitionDownwards(getChild(node, i), thisCase, map);
            }
        } else {
            branchMap.set(node.getNumber(), null, true);
            HashSet<AbstractCase> children = new HashSet<AbstractCase>();
            for(AbstractCase aCase : outbreak.getCases()){
                if(map.get(aCase)==thisCase){
                    children.add(aCase);
                }
            }
            HashSet<AbstractCase> relevantChildren = new HashSet<AbstractCase>(children);
            for(AbstractCase child: children){

                NodeRef caseMRCA = caseMRCA(child);

                //either ALL the tips need to be a descendant of this node, or none. Otherwise not compatible.

                if(directDescendant(node, caseMRCA)){
                    throw new RuntimeException("Starting transmission tree is incompatible with starting phylogeny");
                }

                if(caseMRCA==node){
                    //I'm afraid I must insist...
                    relevantChildren = new HashSet<AbstractCase>();
                    relevantChildren.add(child);
                    break;
                }

                NodeRef currentNode = caseMRCA;
                while(currentNode!=node && currentNode!=null){
                    currentNode = getParent(currentNode);
                }
                if(currentNode==null){
                    relevantChildren.remove(child);
                }
            }
            if(relevantChildren.size()==1){
                //this ends an infection branch
                AbstractCase child = relevantChildren.iterator().next();
                branchMap.set(node.getNumber(), child, true);
            } else {

                //this can't end an infection branch
                branchMap.set(node.getNumber(), thisCase, true);
            }
            for(int i=0; i<getChildCount(node); i++){
                specificallyPartitionDownwards(getChild(node, i), branchMap.get(node.getNumber()), map);
            }
        }

    }


    /*
     todo - The trouble with initialising this without the likelihood class is that lots of starting trees might
     todo - fail. Need to think about how best to deal with this.

     Generally allowCreep is a bad idea, since it tends to place infections after tip times and tip times
     are frequently noninfectiousness times. Might be useful for some pathogens, however.
    */


    private void partitionAccordingToRandomTT(boolean allowCreep){

        System.out.println("Generating a random starting partition of the tree");

        branchMap.setAll(prepareExternalNodeMap(new AbstractCase[getNodeCount()]), true);

        NodeRef root = getRoot();
        randomlyAssignNode(root, allowCreep);

    }


    private AbstractCase randomlyAssignNode(NodeRef node, boolean allowCreep){

        if(isExternal(node)){
            return branchMap.get(node.getNumber());
        } else {

            //If this is a descendant of a case MRCA and an ancestor of one of that case's tips, it must be
            //assigned that case. If it is that of two cases then this tree is incompatible

            ArrayList<AbstractCase> forcedByTopology = new ArrayList<AbstractCase>();

            for(AbstractCase aCase : outbreak.getCases()){
                if(aCase.wasEverInfected) {
                    NodeRef caseMRCA = caseMRCA(aCase, false);
                    HashSet<NodeRef> caseTips = getTipsInThisPartitionElement(aCase);

                    for (NodeRef caseTip : caseTips) {
                        if (directDescendant(node, caseMRCA) && directDescendant(caseTip, node)) {
                            if(!forcedByTopology.contains(aCase)) {
                                forcedByTopology.add(aCase);
                            }
                        }
                    }
                }
            }

            if(forcedByTopology.size()>1){
                throw new RuntimeException("Starting phylogeny is incompatible with this tip partition");
            } else if(forcedByTopology.size()==1){
                branchMap.set(node.getNumber(), forcedByTopology.get(0), true);

                for (int i = 0; i < getChildCount(node); i++) {
                    if(!isExternal(getChild(node, i))){
                        randomlyAssignNode(getChild(node, i), allowCreep);
                    }
                }

                return forcedByTopology.get(0);
            } else {
                //not mandated by the topology
                //three choices - case of child 1, case of child 2, case of parent, unless this is the root

                AbstractCase[] choices = new AbstractCase[2];



                for (int i = 0; i < getChildCount(node); i++) {
                    if(!isExternal(getChild(node, i))){
                        choices[i] = randomlyAssignNode(getChild(node, i), allowCreep);
                    } else {
                        choices[i] = branchMap.get(getChild(node,i).getNumber());
                    }
                }
                //if both choices are null and we're at the root, try again

                while(isRoot(node) && choices[0]==null && choices[1]==null){
                    for (int i = 0; i < getChildCount(node); i++) {
                        if(!isExternal(getChild(node, i))){
                            choices[i] = randomlyAssignNode(getChild(node, i), allowCreep);
                        } else {
                            choices[i] = branchMap.get(getChild(node,i).getNumber());
                        }
                    }
                }

                int randomSelection;
                if (isRoot(node)) {
                    //must make a choice at this point
                    randomSelection = MathUtils.nextInt(2);
                    //they can't both be null
                    if(choices[randomSelection]==null){
                        randomSelection = 1-randomSelection;
                    }
                    AbstractCase winner = choices[randomSelection];
                    fillDownTree(node, winner);
                    return winner;

                } else {
                    randomSelection = MathUtils.nextInt(allowCreep ? 3 : 2);
                }
                if (randomSelection != 2) {
                    AbstractCase winner = choices[randomSelection];
                    AbstractCase loser = choices[1-randomSelection];

                    // check that this isn't going to cause a timings problem

                    if(getNodeHeight(getChild(node, randomSelection)) >
                            loser.getInfectionBranchPosition().getParameterValue(0)
                                    *getBranchLength(getChild(node, 1-randomSelection))
                                    + getNodeHeight(getChild(node, 1-randomSelection))) {
                        winner = loser;
                    }

                    if(winner!=null) {
                        fillDownTree(node, winner);
                    } else {
                        branchMap.set(node.getNumber(), null, true);
                    }

                    return winner;

                } else {
                    //parent partition element will creep to here, but we don't know what that is yet
                    return null;
                }
            }
        }
    }

    private void fillDownTree(NodeRef node, AbstractCase aCase){
        if(branchMap.get(node.getNumber())==null){
            branchMap.set(node.getNumber(), aCase, true);
            for(int i=0; i<2; i++){
                fillDownTree(getChild(node, i), aCase);
            }
        }
    }

}
