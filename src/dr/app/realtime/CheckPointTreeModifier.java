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

import dr.evomodel.tree.TreeModel;
import java.util.ArrayList;

/**
 * @author Guy Baele
 */
public class CheckPointTreeModifier {

    public final static String TREE_UPDATE_OPTION = "JC69Distance";

    private TreeModel treeModel;

    public CheckPointTreeModifier(TreeModel treeModel) {
        this.treeModel = treeModel;
    }

    /**
     * Modifies the current tree by adopting the provided collection of edges and the stored update criterion.
     * @param edges Edges are provided as index: child number; parent: array entry
     * @param nodeHeights Also sets the node heights to the provided values
     * @param childOrder Array that contains whether a child node is left or right child
     */
    public void adoptTreeStructure(int[] edges, double[] nodeHeights, int[] childOrder, String[] taxaNames) {

        //check if the taxa are in the right order, i.e. for now only allow adding taxa at the end of the list
        checkTaxaOrder(taxaNames);

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






        treeModel.endTreeEdit();

        /*

        //set the node heights
        for (int i = 0; i < nodeHeights.length; i++) {
            setNodeHeight(nodes[i], nodeHeights[i]);
        }

        int newRootIndex = -1;
        //now add the parent-child links again to ALL the nodes
        for (int i = 0; i < edges.length; i++) {
            if (edges[i] != -1) {
                nodes[edges[i]].addChild(nodes[i]);
            } else {
                newRootIndex = i;
            }
        }

        //not possible to determine correct ordering of child nodes in the loop where they're being assigned
        //hence perform possible swaps in a separate loop
        for (int i = 0; i < edges.length; i++) {
            if (edges[i] != -1) {
                if (childOrder[i] == 0 && nodes[edges[i]].getChild(0) != nodes[i]) {
                    //swap child nodes
                    TreeModel.Node childOne = nodes[edges[i]].removeChild(0);
                    TreeModel.Node childTwo = nodes[edges[i]].removeChild(1);
                    nodes[edges[i]].addChild(childTwo);
                    nodes[edges[i]].addChild(childOne);
                }
            }
        }

        this.setRoot(nodes[newRootIndex]);

        */

    }

    /**
     * Check whether all the old taxa are still present in the new analysis. Also check whether the
     * additional taxa have been put at the end of the taxa list.
     * @param taxaNames The taxa names (and order) from the old XML, retrieved from the checkpoint file
     * @return An ArrayList containing the names of the additional taxa.
     */
    private ArrayList<String> checkTaxaOrder(String[] taxaNames) {

        int external = treeModel.getExternalNodeCount();

        //TODO With added taxa the order of the external nodes changes?
        //TODO The trick can be to build a starting tree containing only the taxa from the previous analysis.
        //TODO We may need a tool to extract this from the previous analysis and load in into the new BEAST XML
        for (int i = 0; i < external; i++) {
            System.out.println(treeModel.getExternalNode(i));
        }

        for (int i = 0; i < taxaNames.length; i++) {
            if (!taxaNames[i].equals(treeModel.getNodeTaxon(treeModel.getExternalNode(i)).getId())) {
                String error = taxaNames[i] + " vs. " + treeModel.getNodeTaxon(treeModel.getExternalNode(i)).getId();
                throw new RuntimeException("Taxa listed in the wrong order; append new taxa at the end:\n" + error);
            } else {
                System.out.println(taxaNames[i] + " vs. " + treeModel.getNodeTaxon(treeModel.getExternalNode(i)).getId());
            }
        }

        ArrayList<String> newTaxaNames = new ArrayList<String>();

        for (int i = taxaNames.length; i < external; i++) {
            newTaxaNames.add(taxaNames[i]);
        }

        return newTaxaNames;

    }

}
