package dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.evolution.tree.NodeRef;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.AdaptationMode;
import dr.math.MathUtils;

public class ConstrainedTreeOperator extends AbstractTreeOperator {
    public ConstrainedTreeOperator(ConstrainedTreeModel tree,double weight,ConstrainableTreeOperator operator){
        setWeight(weight);
        constrainedTreeModel = tree;
        this.operator = operator;
    }
    @Override
    public String getOperatorName() {
        return "Constrained"+ operator.getOperatorName();
    }

    /**
     * Called by operate(), does the actual operation.
     *
     * @return the hastings ratio
     */
    @Override
    public double doOperation() {
            NodeRef node = constrainedTreeModel.getNode(MathUtils.nextInt(constrainedTreeModel.getNodeCount()));
            TreeModel subtree = constrainedTreeModel.getSubtree(node);

            while(subtree.getExternalNodeCount()<3){
                node = constrainedTreeModel.getNode(MathUtils.nextInt(constrainedTreeModel.getNodeCount()));
                subtree = constrainedTreeModel.getSubtree(node);
            }
        return operator.doOperation(subtree);
    }
    private final ConstrainedTreeModel constrainedTreeModel;
    private final ConstrainableTreeOperator operator;

}
