package dr.evolution.tree;

/**
 * Provides a named "per branch" attribute.
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


    class Wrapper implements NodeAttributeProvider {
        public Wrapper(TreeTrait treeTrait) {
            if (treeTrait.getIntent() != TreeTrait.Intent.BRANCH) {
                throw new RuntimeException("TreeTrait without BRANCH intent wrapped by BranchAttributeProvider");
            }
            this.treeTrait = treeTrait;
        }

        public String[] getNodeAttributeLabel() {
            return new String[] { treeTrait.getTraitName() };
        }

        public String[] getAttributeForNode(Tree tree, NodeRef node) {
            return treeTrait.getTraitString(tree, node);
        }

        private final TreeTrait treeTrait;
    }
}
