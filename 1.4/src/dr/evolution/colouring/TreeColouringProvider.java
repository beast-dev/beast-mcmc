package dr.evolution.colouring;

import dr.evolution.tree.Tree;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public interface TreeColouringProvider {
    TreeColouring getTreeColouring(Tree tree);
}
