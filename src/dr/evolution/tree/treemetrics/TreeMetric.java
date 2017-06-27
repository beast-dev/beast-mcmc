package dr.evolution.tree.treemetrics;

import dr.evolution.tree.Tree;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface TreeMetric {

    double getMetric(Tree tree1, Tree tree2);

}
