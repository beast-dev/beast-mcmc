package dr.evolution.tree;

/**
 * Provides a named "per branch" attrinute.
 * 
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface BranchAttributeProvider {

    /**
     *
     * @return  Descriptive name for the attribute.
     */
    String getBranchAttributeLabel();

    /**
     * @param tree
     * @param node
     * @return    The value of attribute in 'node' inside 'tree'.
     */
    String getAttributeForBranch(Tree tree, NodeRef node);
}
