package dr.evomodel.tree;

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evomodel.continuous.AncestralTaxonInTree;

import java.util.List;

public class DuplicationTreeModel extends AncestralTraitTreeModel{

    public DuplicationTreeModel(String id,
                                MutableTreeModel tree,
                                List<AncestralTaxonInTree> duplications) {
        super(id, tree, duplications);
    }

    protected ShadowNode buildRecursivelyShadowTree(NodeRef originalNode,
                                                    ShadowNode parentNode) {

        final int originalNumber = originalNode.getNumber();
        final int newNumber = mapOriginalToShadowNumber(originalNumber);
        ShadowNode newNode = new ShadowNode(newNumber, originalNode, null);

        ShadowNode recurse = parentNode;
        if (parentNode == null) {

        } else {
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
