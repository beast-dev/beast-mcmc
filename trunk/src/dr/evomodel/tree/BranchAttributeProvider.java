package dr.evomodel.tree;

import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface BranchAttributeProvider {

	String getBranchAttributeLabel();

	String getAttributeForBranch(Tree tree, NodeRef node);
}
