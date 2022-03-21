/*
 * DuplicationTreeModel.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.TaxonList;
import dr.evomodel.continuous.AncestralTaxonInTree;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.*;

/**
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class DuplicationTreeModel extends AncestralTraitTreeModel{

    private final Map<Integer, Parameter> duplicationTimeRatioMap;

    public DuplicationTreeModel(String id,
                                MutableTreeModel tree,
                                List<AncestralTaxonInTree> duplications,
                                List<Parameter> duplicaitonTimeRatios) {
        super(id, tree, duplications);
        this.duplicationTimeRatioMap = new HashMap<>();
        for (int i = 0; i < duplications.size(); i++) {
            duplicationTimeRatioMap.put(getCommonAncestor(tree, duplications.get(i).getTaxonList()).getNumber(), duplicaitonTimeRatios.get(i));
            addVariable(duplicaitonTimeRatios.get(i));
        }
    }

    private NodeRef getCommonAncestor(MutableTreeModel tree, TaxonList descendents) {
        Set<String> leafNodes = new HashSet<>();

        for (int i = 0; i < descendents.getTaxonCount(); i++) {
            leafNodes.add(descendents.getTaxon(i).getId());
        }

        NodeRef mrca = TreeUtils.getCommonAncestorNode(tree, leafNodes);
        return mrca;
    }

    public double getNodeHeight(NodeRef iNode) {
        assert (iNode != null);

        checkShadowTree();

        double height;

        ShadowNode node = (ShadowNode) iNode;
        int originalNumber = node.getOriginalNumber();
        if (originalNumber >= 0) {
            height = treeModel.getNodeHeight(node.getOriginalNode());
        } else {

            final AncestralTaxonInTree ancestor = node.getAncestor();

            if (node.isExternal()) {
                height = getNodeHeight(node.getParent());
            } else {
                final double timeRatio = duplicationTimeRatioMap.get(getInternalChildNode(node).getOriginalNode().getNumber()).getParameterValue(0);
                height = treeModel.getNodeHeight(treeModel.getParent(getInternalChildNode(node).getOriginalNode())) * timeRatio +
                        treeModel.getNodeHeight(getInternalChildNode(node).getOriginalNode()) * (1.0 - timeRatio);
            }
        }
        return height;
    }

    private ShadowNode getInternalChildNode(ShadowNode node) {
        if (node.getChild0().isExternal()) {
            return node.getChild1();
        } else {
            return node.getChild0();
        }
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        fireModelChanged(new TreeChangedEvent.WholeTree());
    }

    protected ShadowNode buildRecursivelyShadowTree(NodeRef originalNode,
                                                    ShadowNode parentNode) {

        final int originalNumber = originalNode.getNumber();
        final int newNumber = mapOriginalToShadowNumber(originalNumber);
        ShadowNode newNode = new ShadowNode(newNumber, originalNode, null);

        ShadowNode recurse = parentNode;
        if (parentNode != null) {
            final boolean isLeftChild = treeModel.getChild(treeModel.getNode(parentNode.getOriginalNumber()), 0) == originalNode;

            if (nodeToClampMap.containsKey(originalNode.getNumber())) {

                List<AncestralTaxonInTree> duplications = nodeToClampMap.get(originalNode.getNumber());
                assert(duplications.size() == 1);


                for (AncestralTaxonInTree duplication : duplications) {

                    ShadowNode newInternalNode = new ShadowNode(externalCount + treeInternalCount + duplication.getIndex(), null, duplication);
                    newInternalNode.setParent(recurse);

                    final int newTipNumber = treeExternalCount + duplication.getIndex();
                    ShadowNode newTipNode = new ShadowNode(newTipNumber, null, duplication);
                    newTipNode.setParent(newInternalNode);
                    newInternalNode.setChild1(newTipNode);

                    if (isLeftChild) {
                        recurse.setChild0(newInternalNode);
                    } else {
                        recurse.setChild1(newInternalNode);
                    }

                    recurse = newInternalNode;

                    storeNode(newTipNode);
                    storeNode(newInternalNode);

                }
            }

            if (isLeftChild) {
                recurse.setChild0(newNode);
            } else {
                recurse.setChild1(newNode);
            }
        }

        newNode.setParent(recurse);
        storeNode(newNode);


        if (!treeModel.isExternal(originalNode)) {
            NodeRef originalChild0 = treeModel.getChild(originalNode, 0);
            NodeRef originalChild1 = treeModel.getChild(originalNode, 1);

            buildRecursivelyShadowTree(originalChild0, newNode);
            buildRecursivelyShadowTree(originalChild1, newNode);
        }

        return newNode;
    }

}
