package dr.evolution.tree;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface NodeAttributeProvider {

    String[] getNodeAttributeLabel();

    String[] getAttributeForNode(Tree tree, NodeRef node);
}
