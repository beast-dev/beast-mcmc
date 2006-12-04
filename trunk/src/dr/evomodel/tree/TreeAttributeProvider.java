package dr.evomodel.tree;

import dr.evolution.tree.Tree;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface TreeAttributeProvider {

	String getAttributeLabel();

	String getAttributeForTree(Tree tree);
}
