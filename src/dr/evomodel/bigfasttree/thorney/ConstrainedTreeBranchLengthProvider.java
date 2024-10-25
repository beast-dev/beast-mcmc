/*
 * ConstrainedTreeBranchLengthProvider.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
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
 *
 */

package dr.evomodel.bigfasttree.thorney;

import dr.evolution.datatype.ContinuousDataType;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.IntegerDataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;

import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;


public class ConstrainedTreeBranchLengthProvider  extends MutationBranchMap.AbstractMutationBranchMap {
    public static final String CONSTRAINED_TREE_BRANCHLENGTH_PROVIDER = "ConstrainedTreeBranchMutationProvider";
    

    public ConstrainedTreeBranchLengthProvider(ConstrainedTreeModel constrainedTreeModel,Tree dataTree,Double scale,double minBranchLength,DataType dataType){
        super(dataType);
        
        boolean discrete = dataType instanceof ContinuousDataType? false:true;

        this.minBranchLength = minBranchLength;

        this.constrainedTreeModel = constrainedTreeModel;

        externalBranchLengths = new double[dataTree.getExternalNodeCount()];
        cladeBranchLengths = new double[dataTree.getInternalNodeCount()];

        Map<String, NodeRef> taxonIdNodeMap = new HashMap<>();
        for (int i = 0; i < dataTree.getExternalNodeCount(); i++) {
            NodeRef node = dataTree.getExternalNode(i);
            taxonIdNodeMap.put(dataTree.getNodeTaxon(node).getId(), node);
        }
        //set up external branchLengths
        for (int i = 0; i < constrainedTreeModel.getExternalNodeCount(); i++) {
            NodeRef node = constrainedTreeModel.getExternalNode(i);
            String taxonId =  constrainedTreeModel.getNodeTaxon(node).getId();
            NodeRef dataNode = taxonIdNodeMap.get(taxonId);

            externalBranchLengths[node.getNumber()] = discrete? Math.round(dataTree.getBranchLength(dataNode)*scale):dataTree.getBranchLength(dataNode)*scale ;
        }
        // need to make subtrees to nodes in the dataTree Subtrees
        Map<BitSet, NodeRef> dataTreeMap = getBitSetNodeMap(dataTree,dataTree);
        Map<BitSet, NodeRef> treeModelMap = getBitSetNodeMap(dataTree,constrainedTreeModel);
        HashMap <NodeRef, NodeRef> dataTreeNodeMap = new HashMap<>();

        for (Map.Entry<BitSet, NodeRef> entry: dataTreeMap.entrySet()){
            dataTreeNodeMap.put(entry.getValue(), treeModelMap.get(entry.getKey()));
        }

        for(int i=0; i<dataTree.getInternalNodeCount(); i++){
            NodeRef dataNode = dataTree.getInternalNode(i);
            NodeRef constrainedNode = dataTreeNodeMap.get(dataNode);
            cladeBranchLengths[constrainedTreeModel.getSubtreeIndex(constrainedNode)] = discrete? Math.round(dataTree.getBranchLength(dataNode)*scale):dataTree.getBranchLength(dataNode)*scale ;
        }
    }
    public ConstrainedTreeBranchLengthProvider(ConstrainedTreeModel constrainedTreeModel,Tree dataTree){
        this(constrainedTreeModel,dataTree,1.0,0.0,IntegerDataType.INSTANCE);
    }

    public double getBranchLength( NodeRef node) {
        if (constrainedTreeModel.isExternal(node)) {
            return externalBranchLengths[node.getNumber()];
        }
        TreeModel subtree = constrainedTreeModel.getSubtree(node);
        NodeRef nodeInSubtree = constrainedTreeModel.getNodeInSubtree(subtree,node);

        if (subtree.isRoot(nodeInSubtree)) {
            int subtreeIndex =  constrainedTreeModel.getSubtreeIndex(node);
            return cladeBranchLengths[subtreeIndex];
        }else{
            return minBranchLength;
        }
    }

    public MutationList getMutations(NodeRef node){
        MutationList mutations = new MutationList.SimpleMutationList(getBranchLength(node));
        return mutations;
    }



    /**
     * Gets a HashMap of clade bitsets to nodes in tree. This is useful for comparing the topology of trees
     * @param referenceTree  the tree that will be used to define taxa and tip numbers
     * @param tree the tree for which clades are being defined
     * @return A HashMap with a BitSet of descendent taxa as the key and a node as value
     */
    private HashMap<BitSet, NodeRef> getBitSetNodeMap(Tree referenceTree,Tree  tree) {
        HashMap<BitSet, NodeRef> map = new HashMap<>();
        addBits(referenceTree,tree,tree.getRoot(),map);
        return map;
    }

    /**
     *  A private recursive function used by getBitSetNodeMap
     *  This is modeled after the addClades in CladeSet and getClades in compatibility statistic
     * @param referenceTree  the tree that will be used to define taxa and tip numbers
     * @param tree the tree for which clades are being defined
     * @param node current node
     * @param map map that is being appended to
     */
    private BitSet addBits(Tree referenceTree, Tree tree, NodeRef node, HashMap map) {
        BitSet bits = new BitSet();
        if (tree.isExternal(node)) {
            String taxonId = tree.getNodeTaxon(node).getId();
            bits.set(referenceTree.getTaxonIndex(taxonId));

        } else {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef node1 = tree.getChild(node, i);
                bits.or(addBits(referenceTree,tree, node1, map));
            }
        }

        map.put(bits, node);
        return bits;
    }

    private final double[] cladeBranchLengths;
    private final double[] externalBranchLengths;

    private final double minBranchLength;
    private ConstrainedTreeModel constrainedTreeModel;



}
