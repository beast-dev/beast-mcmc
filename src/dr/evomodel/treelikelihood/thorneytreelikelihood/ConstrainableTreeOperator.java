package dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.evomodel.tree.TreeModel;

public interface ConstrainableTreeOperator{
    double doOperation(TreeModel tree);
    String getOperatorName();
}
