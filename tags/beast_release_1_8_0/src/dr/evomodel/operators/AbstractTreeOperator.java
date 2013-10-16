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

	private int transitions = 0;

	/**
     * @return the number of transitions since last call to reset().
     */
    public int getTransitions() {
    	return transitions;
    }

    /**
     * Set the number of transitions since last call to reset(). This is used
     * to restore the state of the operator
     *
     * @param transitions number of transition
     */
    public void setTransitions(int transitions) {
    	this.transitions = transitions;
    }

    public double getTransistionProbability() {
        final int accepted = getAcceptCount();
        final int rejected = getRejectCount();
        final int transition = getTransitions();
        return (double) transition / (double) (accepted + rejected);
    }

	/* exchange sub-trees whose root are i and j */
	protected void exchangeNodes(TreeModel tree, NodeRef i, NodeRef j,
	                             NodeRef iP, NodeRef jP) throws OperatorFailedException {

	    tree.beginTreeEdit();
	    tree.removeChild(iP, i);
	    tree.removeChild(jP, j);
	    tree.addChild(jP, i);
	    tree.addChild(iP, j);

        tree.endTreeEdit();
	}

	public void reset() {
        super.reset();
        transitions = 0;
    }

	/**
	 * @param tree   the tree
	 * @param parent the parent
	 * @param child  the child that you want the sister of
	 * @return the other child of the given parent.
	 */
    protected NodeRef getOtherChild(Tree tree, NodeRef parent, NodeRef child) {
        if( tree.getChild(parent, 0) == child ) {
            return tree.getChild(parent, 1);
        } else {
            return tree.getChild(parent, 0);
        }
    }
}
