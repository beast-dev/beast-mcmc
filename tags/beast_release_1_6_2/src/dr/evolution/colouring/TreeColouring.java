package dr.evolution.colouring;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public interface TreeColouring {
    int getColourCount();

    Tree getTree();

    int getNodeColour(NodeRef node);

    BranchColouring getBranchColouring(NodeRef node);

    int getColourChangeCount();

    boolean hasProbability();

    void setLogProbabilityDensity(double p);
}
