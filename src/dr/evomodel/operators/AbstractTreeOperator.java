package dr.evomodel.operators;

import dr.evomodel.tree.TreeModel;
import dr.evolution.tree.*;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public abstract class AbstractTreeOperator extends SimpleMCMCOperator {

	/* exchange subtrees whose root are i and j */
	protected void exchangeNodes(TreeModel tree, NodeRef i, NodeRef j,
	                             NodeRef iP, NodeRef jP) throws OperatorFailedException {

	    tree.beginTreeEdit();
	    tree.removeChild(iP, i);
	    tree.removeChild(jP, j);
	    tree.addChild(jP, i);
	    tree.addChild(iP, j);

	    try {
	        tree.endTreeEdit();
	    } catch (MutableTree.InvalidTreeException ite) {
	        throw new OperatorFailedException(ite.toString());
	    }
	}


	/**
	 * @param tree   the tree
	 * @param parent the parent
	 * @param child  the child that you want the sister of
	 * @return the other child of the given parent.
	 */
	protected NodeRef getOtherChild(Tree tree, NodeRef parent, NodeRef child) {

	    if (tree.getChild(parent, 0) == child) {
	        return tree.getChild(parent, 1);
	    } else {
	        return tree.getChild(parent, 0);
	    }
	}

	/**
	 * Scales the subtree by the given factor starting from the node
	 * subtreeRoot. The tips stay unchanged.
	 *
	 * @param tree        the tree on which the operation is transformed
	 * @param subtreeRoot the root of the subtree to scale
	 * @param factor      the scaling factor
	 */
	protected void scaleSubtree(TreeModel tree, NodeRef subtreeRoot,
	                            double factor) {
	    if (tree.getChildCount(subtreeRoot) > 0) {
	        double height = tree.getNodeHeight(subtreeRoot);
	        double newHeight = height * factor;
	        tree.setNodeHeight(subtreeRoot, newHeight);

	        for (int i = 0; i < tree.getChildCount(subtreeRoot); i++) {
	            if (tree.getChild(subtreeRoot, i) != null) {
	                scaleSubtree(tree, tree.getChild(subtreeRoot, i), factor);
	            } else {
	                scaleSubtree(tree, tree.getChild(subtreeRoot, 1), factor);
	            }
	        }

	        assert (newHeight != Double.NaN);
	        if (factor < 1) {
	            assert (height >= newHeight);
	        } else {
	            assert (height <= newHeight);
	        }
	    }
	}
}
