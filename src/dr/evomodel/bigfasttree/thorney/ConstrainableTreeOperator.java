package dr.evomodel.bigfasttree.thorney;

import dr.evomodel.tree.TreeModel;

public interface ConstrainableTreeOperator{
    double doOperation(TreeModel tree);
    String getOperatorName();
}
