package dr.app.tools.treeannotator;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

import java.util.Set;

/**
 * @author Andrew Rambaut
 * @version $
 */
public interface CladeSystem {
    void add(Tree tree);

    Clade getRootClade();

    void collectAttributes(Set<String> attributeNames, Tree tree);

    void calculateCladeCredibilities(int totalTreesUsed);

    double getLogCladeCredibility(Tree tree);

    int getCladeCount();
}
