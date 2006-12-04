package dr.evomodel.tree;

import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface NodeAttributeProvider {

	String getAttributeLabel();

	String getAttributeForNode(Tree tree, NodeRef node);
}
