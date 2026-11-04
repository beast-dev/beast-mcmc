/*
 * FixedTreeBranchLengthProvider.java
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

//
// Source code recreated from a .class file by IntelliJ IDEA
// (powered by FernFlower decompiler)
//

package dr.evomodel.bigfasttree.thorney;

import dr.evolution.datatype.ContinuousDataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

public class FixedTreeBranchLengthProvider extends MutationBranchMap.AbstractMutationBranchMap {
    public static final String FIXED_TREE_BRANCHLENGTH_PROVIDER = "FixedTreeBranchLengthProvider";
    private final double[] branchLengths;
    private final Tree tree;


    public FixedTreeBranchLengthProvider(Tree fixedTree, Tree dataTree, Double scale, Double minBranchLength, boolean discrete) {
        super(ContinuousDataType.INSTANCE);

        this.tree = fixedTree;
        this.branchLengths = new double[dataTree.getNodeCount()];
        if (this.tree.getNodeCount() != dataTree.getNodeCount()) {
            throw new RuntimeException("The number of nodes in the datatree does not match that in input tree");
        } else {
            Map<String, NodeRef> taxonIdNodeMap = new HashMap();


            for(int i = 0; i < dataTree.getExternalNodeCount(); ++i) {
                NodeRef node = dataTree.getExternalNode(i);
                taxonIdNodeMap.put(dataTree.getNodeTaxon(node).getId(), node);
            }

            for(int i = 0; i < fixedTree.getExternalNodeCount(); ++i) {
                NodeRef node = fixedTree.getExternalNode(i);
                String taxonId = fixedTree.getNodeTaxon(node).getId();
                NodeRef dataNode = taxonIdNodeMap.get(taxonId);
                double observedLength = discrete ? (double) Math.round(dataTree.getBranchLength(dataNode) * scale) : dataTree.getBranchLength(dataNode) * scale;
                this.branchLengths[node.getNumber()] = Math.max(minBranchLength, observedLength);
            }

            Map<BitSet, NodeRef> dataTreeMap = this.getBitSetNodeMap(dataTree, dataTree);
            Map<BitSet, NodeRef> treeModelMap = this.getBitSetNodeMap(dataTree, fixedTree);
            HashMap<NodeRef, NodeRef> dataTreeNodeMap = new HashMap();

            for (Entry <BitSet, NodeRef> entry:
                    dataTreeMap.entrySet()) {
                dataTreeNodeMap.put(entry.getValue(), treeModelMap.get(entry.getKey()));
            }


            for(int i = 0; i < dataTree.getInternalNodeCount(); ++i) {
                NodeRef dataNode = dataTree.getInternalNode(i);
                NodeRef node = dataTreeNodeMap.get(dataNode);
                double observedLength = discrete ? (double)Math.round(dataTree.getBranchLength(dataNode) * scale) : dataTree.getBranchLength(dataNode) * scale;
                this.branchLengths[node.getNumber()] =  Math.max(minBranchLength, observedLength);
            }

        }
    }

    public FixedTreeBranchLengthProvider(Tree tree, Tree dataTree) {
        this(tree, dataTree, 1.0D,0.0D,  true);
    }

    public double getBranchLength(Tree tree, NodeRef node) {
        if (this.tree == tree) {
            return this.branchLengths[node.getNumber()];
        } else {
            throw new RuntimeException("Unrecognized Tree");
        }
    }
        public double getBranchLength (NodeRef node) {
            return getBranchLength(this.tree, node);
       
    }
    public MutationList getMutations(NodeRef node){
        MutationList mutations = new MutationList.SimpleMutationList(getBranchLength(node));
        return mutations;
    }

    private HashMap<BitSet, NodeRef> getBitSetNodeMap(Tree referenceTree, Tree tree) {
        HashMap<BitSet, NodeRef> map = new HashMap();
        this.addBits(referenceTree, tree, tree.getRoot(), map);
        return map;
    }

    private BitSet addBits(Tree referenceTree, Tree tree, NodeRef node, HashMap map) {
        BitSet bits = new BitSet();
        if (tree.isExternal(node)) {
            String taxonId = tree.getNodeTaxon(node).getId();
            bits.set(referenceTree.getTaxonIndex(taxonId));
        } else {
            for(int i = 0; i < tree.getChildCount(node); ++i) {
                NodeRef node1 = tree.getChild(node, i);
                bits.or(this.addBits(referenceTree, tree, node1, map));
            }
        }

        map.put(bits, node);
        return bits;
    }


}
