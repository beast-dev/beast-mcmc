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
}
