/*
 * CheckPointUpdater.java
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

package dr.app.realtime;

import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.Patterns;
import dr.evolution.tree.BranchRates;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.MultiPartitionDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Likelihood;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Guy Baele
 */
public class CheckPointTreeModifier {

    public final static String TREE_UPDATE_OPTION = "JC69Distance";
    public final static Double EPSILON = 0.10;
    public final static Double MIN_DIST = 0.00000001;

    public final static boolean CURRENT_APPROACH = true;

    private TreeModel treeModel;
    private ArrayList<String> newTaxaNames;
    private int[] nodeMap;
    private int additionalTaxa;

    public CheckPointTreeModifier(TreeModel treeModel) {
        this.treeModel = treeModel;
    }

    /**
     * Modifies the current tree by adopting the provided collection of edges and the stored update criterion.
     * @param edges Edges are provided as index: child number; parent: array entry
     * @param nodeHeights Also sets the node heights to the provided values
     * @param childOrder Array that contains whether a child node is left or right child
     * @param taxaNames Taxa names as stored in the checkpoint file
     */
    public void adoptTreeStructure(int[] edges, double[] nodeHeights, int[] childOrder, String[] taxaNames) {

        this.additionalTaxa = treeModel.getExternalNodeCount() - taxaNames.length;

        //check if the taxa are in the right order, i.e. for now only allow adding taxa at the end of the list
        boolean correctOrder = checkTaxaOrder(taxaNames);
        if (!correctOrder) {
            //create a map between the old node order and the new node order
            createNodeMap(taxaNames);
        }

        treeModel.beginTreeEdit();

        int currentTaxa = treeModel.getExternalNodeCount();

        System.out.println("#external nodes = " + currentTaxa);

        //first remove all the child nodes of the internal nodes
        for (int i = treeModel.getExternalNodeCount(); i < treeModel.getNodeCount(); i++) {
            int childCount = treeModel.getChildCount(treeModel.getNodes()[i]);
            for (int j = 0; j < childCount; j++) {
                treeModel.removeChild(treeModel.getNodes()[i], treeModel.getChild(treeModel.getNodes()[i], j));
            }
        }

        //use node mapping to reconstruct the checkpointed tree
        //disregard the added taxa for now

        //start with setting the external node heights
        //this.additionalTaxa can function as the offset in the internal node array
        for (int i = 0; i < (treeModel.getExternalNodeCount()-additionalTaxa); i++) {
            treeModel.setNodeHeight(treeModel.getExternalNode(nodeMap[i]), nodeHeights[i]);
        }
        //set the internal node heights
        //this.additionalTaxa can function as the offset in the internal node array
        for (int i = 0; i < (treeModel.getExternalNodeCount()-additionalTaxa-1); i++) {
            //No just restart counting, will fix later on in the code by adding additionalTaxa variable
            treeModel.setNodeHeight(treeModel.getInternalNode(i), nodeHeights[treeModel.getExternalNodeCount()-additionalTaxa+i]);
        }

        int newRootIndex = -1;
        //now add the parent-child links again to ALL the nodes
        for (int i = 0; i < edges.length; i++) {
            if (edges[i] != -1) {
                //make distinction between external nodes and internal nodes
                if (i < (treeModel.getExternalNodeCount()-additionalTaxa)) {
                    //external node
                    treeModel.addChild(treeModel.getNode(edges[i]+additionalTaxa), treeModel.getExternalNode(nodeMap[i]));
                    System.out.println("external: " + (edges[i]+additionalTaxa) + " > " + nodeMap[i]);
                } else {
                    //internal node
                    treeModel.addChild(treeModel.getNode(edges[i]+additionalTaxa), treeModel.getNode(i+additionalTaxa));
                    System.out.println("internal: " + (edges[i]+additionalTaxa) + " > " + (i+additionalTaxa));
                }
            } else {
                newRootIndex = i;
            }
        }

        //not possible to determine correct ordering of child nodes in the loop where they're being assigned
        //hence perform possible swaps in a separate loop
        for (int i = 0; i < edges.length; i++) {
            if (edges[i] != -1) {
                if (i < (treeModel.getExternalNodeCount()-additionalTaxa)) {
                    if (childOrder[i] == 0 && treeModel.getChild(treeModel.getNode(edges[i] + additionalTaxa), 0) != treeModel.getExternalNode(nodeMap[i])) {
                        NodeRef childOne = treeModel.getChild(treeModel.getNode(edges[i] + additionalTaxa), 0);
                        NodeRef childTwo = treeModel.getChild(treeModel.getNode(edges[i] + additionalTaxa), 1);
                        treeModel.removeChild(treeModel.getNode(edges[i] + additionalTaxa), childOne);
                        treeModel.removeChild(treeModel.getNode(edges[i] + additionalTaxa), childTwo);
                        treeModel.addChild(treeModel.getNode(edges[i] + additionalTaxa), childTwo);
                        treeModel.addChild(treeModel.getNode(edges[i] + additionalTaxa), childOne);
                        System.out.println(i + ": " + edges[i]);
                        System.out.println("  " + treeModel.getChild(treeModel.getNode(edges[i] + additionalTaxa), 0));
                        System.out.println("  " + treeModel.getChild(treeModel.getNode(edges[i] + additionalTaxa), 1));
                    }
                } else {
                    if (childOrder[i] == 0 && treeModel.getChild(treeModel.getNode(edges[i] + additionalTaxa), 0) != treeModel.getNode(i + additionalTaxa)) {
                        NodeRef childOne = treeModel.getChild(treeModel.getNode(edges[i] + additionalTaxa), 0);
                        NodeRef childTwo = treeModel.getChild(treeModel.getNode(edges[i] + additionalTaxa), 1);
                        treeModel.removeChild(treeModel.getNode(edges[i] + additionalTaxa), childOne);
                        treeModel.removeChild(treeModel.getNode(edges[i] + additionalTaxa), childTwo);
                        treeModel.addChild(treeModel.getNode(edges[i] + additionalTaxa), childTwo);
                        treeModel.addChild(treeModel.getNode(edges[i] + additionalTaxa), childOne);
                        System.out.println(i + ": " + edges[i]);
                        System.out.println("  " + treeModel.getChild(treeModel.getNode(edges[i] + additionalTaxa), 0));
                        System.out.println("  " + treeModel.getChild(treeModel.getNode(edges[i] + additionalTaxa), 1));
                    }
                }
            }
        }

        System.out.println("new root index: " + newRootIndex);
        treeModel.setRoot(treeModel.getNode(newRootIndex+additionalTaxa));

        System.out.println(treeModel.toString());

        for (int i = 0; i < edges.length; i++) {
            if (edges[i] != -1) {
                System.out.println(i + ": " + edges[i]);
                System.out.println("  " + treeModel.getChild(treeModel.getNode(edges[i] + additionalTaxa), 0));
                System.out.println("  " + treeModel.getChild(treeModel.getNode(edges[i] + additionalTaxa), 1));
            }
        }

        treeModel.endTreeEdit();

    }

    /**
     * Imports trait information from a file
     * @param edges Edges are provided as index: child number; parent: array entry
     * @param traitModels List of TreeParameterModel object that contain trait information
     * @param traitValues Values to be copied into the List of TreeParameterModel objects
     */
    //TODO Small difference in reconstructed log likelihood, probably due to increased number of rate categories
    public void adoptTraitData(int[] edges, ArrayList<TreeParameterModel> traitModels, double[][] traitValues) {
        int index = 0;
        for (TreeParameterModel tpm : traitModels) {
            int k = 0;
            for (int i = 0; i < edges.length; i++) {
                System.out.println(i + "   " + edges[i]);
                if (edges[i] != -1) {
                    //TODO Seems like I messed up here
                    if (i < (treeModel.getExternalNodeCount()-additionalTaxa)) {
                        tpm.setNodeValue(this.treeModel, this.treeModel.getExternalNode(nodeMap[i]), traitValues[index][k]);
                        System.out.println("Setting external node " + this.treeModel.getExternalNode(nodeMap[i]) + " to " + traitValues[index][k]);
                    } else {
                        tpm.setNodeValue(this.treeModel, this.treeModel.getNode(i + additionalTaxa), traitValues[index][k]);
                        System.out.println("Setting internal node " + this.treeModel.getNode(i + additionalTaxa) + " to " + traitValues[index][k]);
                    }
                } else {
                    k--;
                }
                k++;
            }
            //set trait of remaining internal nodes to -1.0
            for (int i = 0; i < additionalTaxa; i++) {
                tpm.setNodeValue(this.treeModel, treeModel.getNode(treeModel.getNodeCount()-1-i), -1.0);
            }
            //set trait of newly added external taxa to -1.0
            int shift = 0;

            for (int i = 0; i < (treeModel.getExternalNodeCount()-additionalTaxa); i++) {
                System.out.println("i = " + i + " ; nodeMap[i] = " + nodeMap[i]);
            }

            System.out.println();
            for (String name : newTaxaNames) {
                System.out.println("new taxon: " + name);
            }

            int externalMissing = 0;
            for (String name : newTaxaNames) {
                for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
                    if (treeModel.getNodeTaxon(treeModel.getExternalNode(i)).getId().equals(name)) {
                        externalMissing++;
                        tpm.setNodeValue(this.treeModel, this.treeModel.getExternalNode(i), -1.0);
                    }
                }
            }
            System.out.println("External nodes with trait set to -1.0 = " + externalMissing + "\n");


            /*for (int i = 0; i < (treeModel.getExternalNodeCount()-additionalTaxa); i++) {
                System.out.println("i = " + i + " ; nodeMap[i] = " + nodeMap[i]);
                if (i != (nodeMap[i]-shift)) {
                    int difference = nodeMap[i] - i;
                    shift = difference;
                    tpm.setNodeValue(this.treeModel, this.treeModel.getExternalNode(nodeMap[i]-1), -1.0);
                    System.out.println("Setting external node: " + (nodeMap[i]-1));
                    externalMissing++;
                }
            }
            System.out.println("External node with trait set to -1.0 = " + externalMissing);
            index++;*/
        }
    }

    /**
     * The newly added taxa still need to be provided with trait values if there are any.
     * @param traitModels List of the trait models for which trait values need to be imputed / interpolated.
     */
    public void interpolateTraitValues(ArrayList<TreeParameterModel> traitModels) {
        System.out.println();
        for (TreeParameterModel tpm : traitModels) {
            int numberOfInterpolations = 0;
            for (int i = 0; i < treeModel.getNodeCount(); i++) {
                if (tpm.getNodeValue(treeModel, treeModel.getNode(i)) == -1.0) {
                    System.out.println("> Current trait = -1.0 for node: " + treeModel.getNode(i));
                    numberOfInterpolations++;
                    double newValue = -1.0;
                    //get trait value from sibling
                    NodeRef parent = treeModel.getParent(treeModel.getNode(i));
                    if (parent != null) {
                        for (int j = 0; j < treeModel.getChildCount(parent); j++) {
                            NodeRef child = treeModel.getChild(parent, j);
                            if (tpm.getNodeValue(treeModel, child) != -1.0) {
                                tpm.setNodeValue(treeModel, treeModel.getNode(i), tpm.getNodeValue(treeModel, child) + 1.0);
                                System.out.println("Checking sibling trait.");
                                System.out.println("Setting node trait for node " + treeModel.getNode(i) + " to " + tpm.getNodeValue(treeModel, treeModel.getNode(i)));
                                break;
                            }
                        }
                    }
                    //if not successful, get trait from its parent
                    if (tpm.getNodeValue(treeModel, treeModel.getNode(i)) == -1.0) {
                        NodeRef currentNode = treeModel.getNode(i);
                        //TODO Check for situations where no proper trait can be found
                        while (currentNode != treeModel.getRoot() && tpm.getNodeValue(treeModel, currentNode) == -1.0) {
                            currentNode = treeModel.getParent(currentNode);
                        }
                        tpm.setNodeValue(treeModel, treeModel.getNode(i), tpm.getNodeValue(treeModel, currentNode) + 1.0);
                        System.out.println("Checking parent trait.");
                        System.out.println("Setting node trait for node " + treeModel.getNode(i) + " to " + tpm.getNodeValue(treeModel, currentNode));
                    }
                    //adjust the other trait values after a trait has been imputed
                    System.out.println("tree root number: " + treeModel.getRoot().getNumber());
                    double maxValue = -1.0;
                    for (int j = 0; j < treeModel.getNodeCount(); j++) {
                        if (treeModel.getNode(j) != treeModel.getNode(i)) {
                            double nodeJValue = tpm.getNodeValue(treeModel, treeModel.getNode(j));
                            //special case in
                            if (nodeJValue >= tpm.getNodeValue(treeModel, treeModel.getNode(i)) && j != treeModel.getRoot().getNumber()) {
                                System.out.print("Updating trait from " + nodeJValue);
                                tpm.setNodeValue(treeModel, treeModel.getNode(j), nodeJValue + 1.0);
                                System.out.println(" to " + (nodeJValue+1.0) + " (j = " + j + ")");
                                if ((nodeJValue + 1.0) > maxValue) {
                                    maxValue = nodeJValue + 1.0;
                                }
                            }
                        }
                    }

                    System.out.println("Maximal trait value = " + maxValue);
                }
                /*if (tpm.getNodeValue(treeModel, treeModel.getNode(i)) == -1.0) {
                    for (int j = 0; j < treeModel.getNodeCount(); j++) {
                        if (treeModel.getNode(j) != treeModel.getNode(i)) {
                            if (tpm.getNodeValue(treeModel, treeModel.getNode(j)) >= tpm.getNodeValue(treeModel, treeModel.getNode(i))) {
                                System.out.print("Updating trait from " + tpm.getNodeValue(treeModel, treeModel.getNode(j)));
                                tpm.setNodeValue(treeModel, treeModel.getNode(j), tpm.getNodeValue(treeModel, treeModel.getNode(j)) + 1.0);
                                System.out.println(" to " + tpm.getNodeValue(treeModel, treeModel.getNode(j)));
                            }
                        }
                    }
                }*/
            }
            System.out.println("Number of interpolations: " + numberOfInterpolations);
        }
        System.out.println("Done.\n");
    }

    /**
     * Add the remaining taxa, which can be identified through the TreeDataLikelihood XML elements.
     */
    public ArrayList<NodeRef> incorporateAdditionalTaxa(CheckPointUpdaterApp.UpdateChoice choice, BranchRates rateModel) {

        System.out.println("Tree before adding taxa:\n" + treeModel.toString() + "\n");

        ArrayList<NodeRef> newTaxaNodes = new ArrayList<NodeRef>();
        for (String str : newTaxaNames) {
            for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
                if (treeModel.getNodeTaxon(treeModel.getExternalNode(i)).getId().equals(str)) {
                    newTaxaNodes.add(treeModel.getExternalNode(i));
                    //always take into account Taxon dates vs. dates set through a TreeModel
                    System.out.println(treeModel.getNodeTaxon(treeModel.getExternalNode(i)).getId() + " with height "
                            + treeModel.getNodeHeight(treeModel.getExternalNode(i)) + " or " + treeModel.getNodeTaxon(treeModel.getExternalNode(i)).getHeight());
                }
            }
        }
        System.out.println("newTaxaNodes length = " + newTaxaNodes.size());

        ArrayList<Taxon> currentTaxa = new ArrayList<Taxon>();
        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
            boolean taxonFound = false;
            for (String str : newTaxaNames) {
                if (str.equals((treeModel.getNodeTaxon(treeModel.getExternalNode(i))).getId())) {
                    taxonFound = true;
                }
            }
            if (!taxonFound) {
                System.out.println("Adding " + treeModel.getNodeTaxon(treeModel.getExternalNode(i)).getId() + " to list of current taxa");
                currentTaxa.add(treeModel.getNodeTaxon(treeModel.getExternalNode(i)));
            }
        }
        System.out.println("Current taxa count = " + currentTaxa.size());

        //iterate over both current taxa and to be added taxa
        boolean originTaxon = true;
        for (Taxon taxon : currentTaxa) {
            if (taxon.getHeight() == 0.0) {
                originTaxon = false;
                System.out.println("Current taxon " + taxon.getId() + " has node height 0.0");
            }
        }
        for (NodeRef newTaxon : newTaxaNodes) {
            if (treeModel.getNodeTaxon(newTaxon).getHeight() == 0.0) {
                originTaxon = false;
                System.out.println("New taxon " + treeModel.getNodeTaxon(newTaxon).getId() + " has node height 0.0");
            }
        }

        //check the Tree(Data)Likelihoods in the connected set of likelihoods
        //focus on TreeDataLikelihood, which has getTree() to get the tree for each likelihood
        //also get the DataLikelihoodDelegate from TreeDataLikelihood
        ArrayList<TreeDataLikelihood> likelihoods = new ArrayList<TreeDataLikelihood>();
        ArrayList<Tree> trees = new ArrayList<Tree>();
        ArrayList<DataLikelihoodDelegate> delegates = new ArrayList<DataLikelihoodDelegate>();
        for (Likelihood likelihood : Likelihood.CONNECTED_LIKELIHOOD_SET) {
            if (likelihood instanceof TreeDataLikelihood) {
                likelihoods.add((TreeDataLikelihood)likelihood);
                trees.add(((TreeDataLikelihood) likelihood).getTree());
                delegates.add(((TreeDataLikelihood) likelihood).getDataLikelihoodDelegate());
            }
        }

        //suggested to go through TreeDataLikelihoodParser and give it an extra option to create a HashMap
        //keyed by the tree; am currently not overly fond of this approach
        ArrayList<PatternList> patternLists = new ArrayList<PatternList>();
        for (DataLikelihoodDelegate del : delegates) {
            if (del instanceof BeagleDataLikelihoodDelegate) {
                patternLists.add(((BeagleDataLikelihoodDelegate) del).getPatternList());
            } else if (del instanceof MultiPartitionDataLikelihoodDelegate) {
                MultiPartitionDataLikelihoodDelegate mpdld = (MultiPartitionDataLikelihoodDelegate)del;
                List<PatternList> list = mpdld.getPatternLists();
                for (PatternList pList : list) {
                    patternLists.add(pList);
                }
            }
        }

        if (patternLists.size() == 0) {
            throw new RuntimeException("No patterns detected. Please make sure the XML file is BEAST 1.9 compatible.");
        }

        //aggregate all patterns to create distance matrix
        //TODO What about different trees for different partitions?
        Patterns patterns = new Patterns(patternLists.get(0));
        if (patternLists.size() > 1) {
            for (int i = 1; i < patternLists.size(); i++) {
                patterns.addPatterns(patternLists.get(i));
            }
        }

        //set the patterns for the distance matrix computations
        choice.setPatterns(patterns);

        //add new taxa one at a time
        System.out.println("Adding " + newTaxaNodes.size() + " taxa ...");

        if (CURRENT_APPROACH) {
            for (NodeRef newTaxon : newTaxaNodes) {
                treeModel.setNodeHeight(newTaxon, treeModel.getNodeTaxon(newTaxon).getHeight());
                System.out.println("\nadding Taxon: " + newTaxon + " (height = " + treeModel.getNodeHeight(newTaxon) + ")");
                //check if this taxon has a more recent sampling date than all other nodes in the current TreeModel
                double offset = checkCurrentTreeNodes(newTaxon, treeModel.getRoot());
                System.out.println("Sampling date offset when adding " + newTaxon + " = " + offset);
                //if so, update all nodes current in the tree (i.e. recursively from the root)
                //AND set its current node height to 0.0 IF no originTaxon has been found
                if (offset < 0.0) {
                    if (!originTaxon) {
                        System.out.println("Updating all node heights with offset " + Math.abs(offset));
                        updateAllTreeNodes(Math.abs(offset), treeModel.getRoot());
                        treeModel.setNodeHeight(newTaxon, 0.0);
                    }
                } else if (offset == 0.0) {
                    if (!originTaxon) {
                        treeModel.setNodeHeight(newTaxon, 0.0);
                    }
                }
                //get the closest Taxon to the Taxon that needs to be added
                //take into account which taxa can currently be chosen
                Taxon closest = choice.getClosestTaxon(treeModel.getNodeTaxon(newTaxon), currentTaxa);
                System.out.println("\nclosest Taxon: " + closest + " with original height: " + closest.getHeight());
                //get the distance between these two taxa
                double distance = choice.getDistance(treeModel.getNodeTaxon(newTaxon), closest);
                if(distance == 0.0){
                    distance = MIN_DIST;
                }
                System.out.println("at distance: " + distance);
                //find the NodeRef for the closest Taxon (do not rely on node numbering)
                NodeRef closestRef = null;
                //careful with node numbering and subtract number of new taxa
                for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
                    if (treeModel.getNodeTaxon(treeModel.getExternalNode(i)) == closest) {
                        closestRef = treeModel.getExternalNode(i);
                    }
                }
                System.out.println(closestRef + " with height " + treeModel.getNodeHeight(closestRef));
                //System.out.println("trying to set node height: " + closestRef + " from " + treeModel.getNodeHeight(closestRef) + " to " + closest.getHeight());
                //treeModel.setNodeHeight(closestRef, closest.getHeight());
                double timeForDistance = distance / rateModel.getBranchRate(treeModel, closestRef);
                System.out.println("timeForDistance = " + timeForDistance);
                //get parent node of branch that will be split
                NodeRef parent = treeModel.getParent(closestRef);
                //child node of branch that will be split (will be assigned value later)
                NodeRef splitBranchChild;

                //determine height of new node
                double insertHeight;
                if (treeModel.getNodeHeight(closestRef) == treeModel.getNodeHeight(newTaxon)) {
                    insertHeight = treeModel.getNodeHeight(closestRef) + timeForDistance / 2.0;
                    System.out.println("treeModel.getNodeHeight(closestRef) == treeModel.getNodeHeight(newTaxon): " + insertHeight);
                    splitBranchChild = closestRef;

                    if (insertHeight >= treeModel.getNodeHeight(parent)) {
                        while (insertHeight >= treeModel.getNodeHeight(parent)) {
                            if(treeModel.getParent(parent)==null){
                                // Use this insertHeight value in case parent doesn't have parent
                                // Otherwise, move up tree
                                insertHeight = treeModel.getNodeHeight(splitBranchChild) + EPSILON * (treeModel.getNodeHeight(parent) - treeModel.getNodeHeight(splitBranchChild));
                                break;
                            }else{
                                splitBranchChild = parent;
                                parent = treeModel.getParent(splitBranchChild);
                            }
                        }
                    }
                } else {
                    double remainder = (timeForDistance - Math.abs(treeModel.getNodeHeight(closestRef) - treeModel.getNodeHeight(newTaxon))) / 2.0;
                    if (remainder > 0) {
                        insertHeight = Math.max(treeModel.getNodeHeight(closestRef), treeModel.getNodeHeight(newTaxon)) + remainder;
                        System.out.println("remainder > 0: " + insertHeight);
                        splitBranchChild = closestRef;

                        if (insertHeight >= treeModel.getNodeHeight(parent)) {
                            while (insertHeight >= treeModel.getNodeHeight(parent)) {
                                if(treeModel.getParent(parent)==null){
                                    // Use this insertHeight value in case parent doesn't have parent
                                    // Otherwise, move up tree
                                    insertHeight = treeModel.getNodeHeight(splitBranchChild) + EPSILON * (treeModel.getNodeHeight(parent) - treeModel.getNodeHeight(splitBranchChild));
                                    break;
                                }else {
                                    splitBranchChild = parent;
                                    parent = treeModel.getParent(splitBranchChild);
                                }
                            }
                        }
                    } else {
                        // Come up with better way to handle this?
                        insertHeight = EPSILON * (treeModel.getNodeHeight(parent) - Math.max(treeModel.getNodeHeight(closestRef), treeModel.getNodeHeight(newTaxon)));
                        insertHeight += Math.max(treeModel.getNodeHeight(closestRef), treeModel.getNodeHeight(newTaxon));
                        System.out.println("remainder <= 0: " + insertHeight);
                        splitBranchChild = closestRef;
                    }
                }

                System.out.println("insert at height: " + insertHeight);
                //pass on all the necessary variables to a method that adds the new taxon to the tree
               // addTaxonAlongBranch(newTaxon, parent, closestRef, insertHeight);
                addTaxonAlongBranch(newTaxon, parent, splitBranchChild, insertHeight);
                //option to print tree after each taxon addition
                System.out.println("\nTree after adding taxon " + newTaxon + ":\n" + treeModel.toString());
                //add newly added Taxon to list of current taxa
                currentTaxa.add(treeModel.getNodeTaxon(newTaxon));
            }
        } else {

            //test section for Ebola problem
            throw new RuntimeException("Not yet implemented ...");

        }

        //System.out.println(treeModel.toString());

        return newTaxaNodes;
    }

    /**
     * Add a given offset to all node height (both internal and external) of a tree
     * @param offset The offset to add to the node height
     * @param currentRoot Root of the tree
     */
    private void updateAllTreeNodes(double offset, NodeRef currentRoot) {
        treeModel.setNodeHeight(currentRoot, treeModel.getNodeHeight(currentRoot) + offset);
        for (int i = 0; i < treeModel.getChildCount(currentRoot); i++) {
            updateAllTreeNodes(offset, treeModel.getChild(currentRoot, i));
        }
    }

    /**
     * Check if the taxon to be added has a more recent sampling date than all external nodes currently in the tree.
     * @param newTaxon The taxon to be added
     * @param currentRoot The current connected tree that does not yet contain the @newTaxon taxon
     * @return Offset of the node height of the new taxon versus the current most recent sampling date
     */
    private double checkCurrentTreeNodes(NodeRef newTaxon, NodeRef currentRoot) {
        /*System.out.println(currentRoot);
        if (treeModel.getChildCount(currentRoot) == 0) {
            System.out.println("root date = " + treeModel.getNodeTaxon(currentRoot).getDate().getTimeValue());
        }*/
        if (treeModel.getChildCount(currentRoot) == 0) {
            return treeModel.getNodeTaxon(currentRoot).getDate().getTimeValue() - treeModel.getNodeTaxon(newTaxon).getDate().getTimeValue();
        } else {
            double maximum = -Double.MAX_VALUE;
            for (int i = 0; i < treeModel.getChildCount(currentRoot); i++) {
                maximum = Math.max(checkCurrentTreeNodes(newTaxon, treeModel.getChild(currentRoot, i)), maximum);
            }
            return maximum;
        }
    }

    /**
     * Adds a taxon along a branch specified by its parent and child node, at a specified height
     * @param newTaxon The new Taxon to be added to the tree
     * @param parentNode The parent node of the branch along which the new taxon is to be added
     * @param childNode The child of the branch along which the new taxon is to be added
     * @param insertHeight The height of the new internal node that will be created
     */
    private void addTaxonAlongBranch(NodeRef newTaxon, NodeRef parentNode, NodeRef childNode, double insertHeight) {
        treeModel.beginTreeEdit();

        //remove child node that corresponds to childNode argument
        treeModel.removeChild(parentNode, childNode);
        //add a currently unused internal node as the new child node
        NodeRef internalNode = null;
        for (int i = 0; i < treeModel.getInternalNodeCount(); i++) {
            if (treeModel.getChildCount(treeModel.getInternalNode(i)) == 0 && treeModel.getParent(treeModel.getInternalNode(i)) == null) {
                internalNode = treeModel.getInternalNode(i);
                System.out.println("\ninternal node found: " + internalNode.getNumber());
                break;
            }
        }
        if (internalNode == null) {
            throw new RuntimeException("No free internal node found for adding a new taxon.");
        }
        //add internal node as a child of the old parent node
        treeModel.addChild(parentNode, internalNode);
        //add the old child node as a child of the new internal node
        treeModel.addChild(internalNode, childNode);
        //add the new taxon as a child of the new internal node
        treeModel.addChild(internalNode, newTaxon);
        //still need to set the height of the new internal node
        treeModel.setNodeHeight(internalNode, insertHeight);

        System.out.println("\nparent node (" + parentNode + ") height: " + treeModel.getNodeHeight(parentNode));
        System.out.println("child node (" + childNode + ") height: " + treeModel.getNodeHeight(childNode));
        System.out.println("internal node (" + internalNode + ") height: " + treeModel.getNodeHeight(internalNode));

        treeModel.endTreeEdit();
    }

    /**
     * Check whether all the old taxa are still present in the new analysis. Also check whether the
     * additional taxa have been put at the end of the taxa list.
     * @param taxaNames The taxa names (and order) from the old XML, retrieved from the checkpoint file
     * @return An ArrayList containing the names of the additional taxa.
     */
    private boolean checkTaxaOrder(String[] taxaNames) {

        int external = treeModel.getExternalNodeCount();

        newTaxaNames = new ArrayList<String>();

        for (int i = 0; i < external; i++) {
            int j = 0;
            boolean taxonFound = false;
            while (!taxonFound) {
                if (treeModel.getNodeTaxon(treeModel.getExternalNode(i)).getId().equals(taxaNames[j])) {
                    taxonFound = true;
                } else {
                    j++;
                    if (j >= taxaNames.length) {
                        newTaxaNames.add(treeModel.getNodeTaxon(treeModel.getExternalNode(i)).getId());
                        taxonFound = true;
                    }
                }
            }
        }

        System.out.println();
        for (String str : newTaxaNames) {
            System.out.println("New taxon found: " + str);
        }
        System.out.println();

        //With added taxa the order of the external nodes changes
        //Aim here is to build a starting tree containing only the taxa from the previous analysis.

        for (int i = 0; i < taxaNames.length; i++) {
            if (!taxaNames[i].equals(treeModel.getNodeTaxon(treeModel.getExternalNode(i)).getId())) {
                String error = taxaNames[i] + " vs. " + treeModel.getNodeTaxon(treeModel.getExternalNode(i)).getId();
                System.out.println("Taxa listed in the wrong order; append new taxa at the end:\n" + error);
                //throw new RuntimeException("Taxa listed in the wrong order; append new taxa at the end:\n" + error);
                return false;
            } else {
                System.out.println(taxaNames[i] + " vs. " + treeModel.getNodeTaxon(treeModel.getExternalNode(i)).getId());
            }
        }

        return true;

    }

    private void createNodeMap(String[] taxaNames) {

        System.out.println("Creating a node mapping:");

        int external = treeModel.getExternalNodeCount();

        nodeMap = new int[external];
        for (int i = 0; i < taxaNames.length; i++) {
            for (int j = 0; j < external; j++) {
                if (taxaNames[i].equals(treeModel.getNodeTaxon(treeModel.getExternalNode(j)).getId())) {
                    //taxon found
                    nodeMap[i] = j;
                }
            }
        }

    }

}
