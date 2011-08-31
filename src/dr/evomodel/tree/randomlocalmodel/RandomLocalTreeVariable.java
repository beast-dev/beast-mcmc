package dr.evomodel.tree.randomlocalmodel;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;

/**
 * This interface provides for a stochastic variable selection model on a tree.
 * The model has two values at each node in the tree: an indicator and a real-valued parameter.
 * The real-valued parameter is only used if the indicator is true.
 *
 * @author Alexei Drummond
 */
public interface RandomLocalTreeVariable {

    /**
     * @param tree the tree
     * @param node the node to retrieve the variable of
     * @return the raw real-valued variable at this node
     */
    double getVariable(TreeModel tree, NodeRef node);

    /**
     * @param tree the tree
     * @param node the node
     * @return true of the variable at this node is included in function, thus representing a change in the
     *         function looking down the tree.
     */
    boolean isVariableSelected(TreeModel tree, NodeRef node);
}
