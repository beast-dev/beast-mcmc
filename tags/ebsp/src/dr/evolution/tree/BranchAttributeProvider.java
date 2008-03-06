package dr.evolution.tree;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface BranchAttributeProvider {

    String getBranchAttributeLabel();

    String getAttributeForBranch(Tree tree, NodeRef node);
}
