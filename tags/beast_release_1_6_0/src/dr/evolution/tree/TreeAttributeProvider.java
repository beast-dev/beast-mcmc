package dr.evolution.tree;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public interface TreeAttributeProvider {

	String[] getTreeAttributeLabel();

	String[] getAttributeForTree(Tree tree);
}
