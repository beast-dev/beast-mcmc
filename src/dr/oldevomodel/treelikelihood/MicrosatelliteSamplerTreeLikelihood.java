/*
 * MicrosatelliteSamplerTreeLikelihood.java
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

package dr.oldevomodel.treelikelihood;

import dr.evomodel.tree.MicrosatelliteSamplerTreeModel;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.oldevomodel.substmodel.MicrosatelliteModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evolution.tree.NodeRef;
import dr.inference.model.Model;

import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 *
 * Treelikelihood that allows ancestral state sampling for microsatellite data.
 *
 */
@Deprecated // Switching to BEAGLE
public class MicrosatelliteSamplerTreeLikelihood extends AbstractTreeLikelihood{

    MicrosatelliteSamplerTreeModel treeMicrosatSamplerModel;
    MicrosatelliteModel microsatelliteModel;
    BranchRateModel branchRateModel;
    double logL = 0.0;
    double storedLogL = 0.0;
    boolean modelChanged = false;
    private ArrayList<Integer> updateAllList;
    private ArrayList<Integer> updatedNodeList;
    public MicrosatelliteSamplerTreeLikelihood(MicrosatelliteSamplerTreeModel treeMicrosatSamplerModel,
                                               MicrosatelliteModel microsatelliteModel,
                                               BranchRateModel branchRate){
        super("MicrosatelliteSamplerTreeLikelihood",
                treeMicrosatSamplerModel.getMicrosatPattern(),
                treeMicrosatSamplerModel.getTreeModel());
        this.treeMicrosatSamplerModel = treeMicrosatSamplerModel;
        this.microsatelliteModel = microsatelliteModel;
        this.branchRateModel = branchRate;

        updateAllList = new ArrayList<Integer>();

        for(int i = 0; i < nodeCount; i++){
            updateAllList.add(i);
        }
        updatedNodeList = updateAllList;

        addModel(this.branchRateModel);
        addModel(this.treeMicrosatSamplerModel);
        addModel(this.microsatelliteModel);


    }


    protected void handleModelChangedEvent(Model model, Object object, int index) {

        if(updatedNodeList == updateAllList){

        }else if (model == treeModel) {

            if (object instanceof TreeModel.TreeChangedEvent) {

                if(((TreeModel.TreeChangedEvent) object).areAllInternalHeightsChanged()){
                    updateAllNodes();
                } else if (((TreeChangedEvent) object).isNodeChanged()) {
                    // If a node event occurs the node and its two child nodes
                    // are flagged for updating (this will result in everything
                    // above being updated as well. Node events occur when a node
                    // is added to a branch, removed from a branch or its height or
                    // rate changes.
                    updateNodeAndChildren(((TreeChangedEvent) object).getNode());
                } else if (((TreeChangedEvent) object).isTreeChanged()) {
                    // Full tree events result in a complete updating of the tree likelihood
                    // Currently this event type is not used.
                    System.err.println("Full tree update event - these events currently aren't used\n" +
                            "so either this is in error or a new feature is using them so remove this message.");
                    updateAllNodes();
                } else {
                    // Other event types are ignored (probably trait changes).
                    System.err.println("Another tree event has occured (possibly a trait change).");
                }
            }

        } else if (model == treeMicrosatSamplerModel){

            if(treeMicrosatSamplerModel.areInternalNodesChanged()) {

                updateNodeAndChildren((treeMicrosatSamplerModel.getTreeModel()).getNode(index));
                treeMicrosatSamplerModel.setInternalNodesChanged(false);
            }

        } else if (model == microsatelliteModel) {
            if(microsatelliteModel.isModelUpdated()){
                updateAllNodes();
            }
        } else if (model == branchRateModel) {
            if (index == -1) {
                updateAllNodes();
            } else {
                updateNode(treeModel.getNode(index));
            }
        } else {
            throw new RuntimeException("Unknown componentChangedEvent");
        }

        modelChanged = true;
        super.handleModelChangedEvent(model, object, index);
    }



    public boolean hasModelChanged(){
        return modelChanged;
    }

    /**
     * Nodes are added to upadteNodeList if it is to be updated.
     */
    protected void updateNodeAndChildren(NodeRef node) {
        int nodeNum = node.getNumber();

        if(!updateNode[nodeNum]){
            updateNode[nodeNum] = true;
            updatedNodeList.add(nodeNum);
        }

        for (int i = 0; i < treeModel.getChildCount(node); i++) {
            int childNodeNum = (treeModel.getChild(node, i)).getNumber();
            if(!updateNode[childNodeNum]){
                updateNode[childNodeNum] = true;
                updatedNodeList.add(childNodeNum);
            }
        }
        likelihoodKnown = false;
    }

    /**
     *
     * If an event requires all the nodes to be updated
     * then the updateNodeList will be referenced to
     * updateAllList, which is a list that contains all the node
     * in the tree.
     *
     */
    protected void updateAllNodes() {
        updatedNodeList = updateAllList;
        likelihoodKnown = false;
    }


    public double calculateLogLikelihood(){

        traverse();

        //used for likelihood calculation.
        double temp = 0.0;
        double temp1 = 0.0;
        double temp2 = 0.0;

        //get the number of nodes that needs to be updated.
        int updateNum = updatedNodeList.size();

        //iterate through the nodes to be updated
        for(int i = 0; i < updateNum; i++){

            //get the node number of the node to be updated
            int nodeNum = updatedNodeList.get(i);

            //calculate the sum of the previous log probabilities of the nodes flagged for update
            //only calculate this if the number of nodes to be updated is less than the total number of nodes in the tree and
            //when there is a change in the tree.
            if(modelChanged && updatedNodeList != updateAllList)
                temp1 += treeMicrosatSamplerModel.getStoredLogBranchLikelihood(nodeNum);

            //the sum of the probabilities of the nodes that have been updated.
            temp2 += treeMicrosatSamplerModel.getLogBranchLikelihood(nodeNum);

            //indicate that the update is complete.
            updateNode[nodeNum] = false;
        }

        //subtract the sum of the previous log probabilities of the updated nodes from the previous likelihood.
        if(modelChanged && updatedNodeList != updateAllList)
            temp = logL - temp1;


        //add the new likelihoods of the updated node to get the current likelihood.
        logL = temp + temp2;

        //indicate that the completion of calculating the likelihood after a modelChange
        modelChanged = false;

        //clear the list of the nodes to be updated
        updatedNodeList = new ArrayList<Integer>();


      
        return logL;

    }

    protected void storeState() {
        storedLogL = logL;

        super.storeState();

    }

    /**
     * Restore the additional stored state
     */
    protected void restoreState() {
        logL = storedLogL;
        super.restoreState();

    }


    private void traverse(){
        TreeModel tree = treeMicrosatSamplerModel.getTreeModel();

        int updateNum = updatedNodeList.size();
        for(int i = 0; i < updateNum; i++){
                NodeRef node = tree.getNode(updatedNodeList.get(i));
                int nodeState =  treeMicrosatSamplerModel.getNodeValue(node);
                if(!tree.isRoot(node)){

                    NodeRef parent = tree.getParent(node);
                    int parentState = treeMicrosatSamplerModel.getNodeValue(parent);

                    double branchLength = tree.getBranchLength(node)*branchRateModel.getBranchRate(tree,node);
                    double nodePr = microsatelliteModel.getLogOneTransitionProbabilityEntry(branchLength, parentState, nodeState);

                    treeMicrosatSamplerModel.setLogBranchLikelihood(node, nodePr);
                }else{
                    double logEqFreq = Math.log(microsatelliteModel.getStationaryDistribution()[nodeState]);
                    treeMicrosatSamplerModel.setLogBranchLikelihood(node, logEqFreq);

                }


        }

    }




}
