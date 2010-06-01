package dr.evolution.tree;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface NodeAttributeProvider {

    /**
     *
     * @return All attributes names
     */
    String[] getNodeAttributeLabel();

    /**
     *
     * @param tree
     * @param node
     * @return  All values of node inside tree
     */
    String[] getAttributeForNode(Tree tree, NodeRef node);

    class Wrapper implements NodeAttributeProvider {
        public Wrapper(TreeTrait treeTrait) {
            if (treeTrait.getIntent() != TreeTrait.Intent.NODE) {
                throw new RuntimeException("TreeTrait without NODE intent wrapped by NodeAttributeProvider");
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
