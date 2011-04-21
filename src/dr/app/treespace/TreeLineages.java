package dr.app.treespace;

import jebl.evolution.graphs.Node;
import jebl.evolution.trees.RootedTree;

import java.awt.*;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class TreeLineages {

    public TreeLineages() {
    }

    public TreeLineages(final List<RootedTree> trees) {
        addTrees(trees);
    }

    public void addTrees(final List<RootedTree> trees) {
        for (RootedTree tree : trees) {
        }
    }

    public void addTree(RootedTree tree) {
        addTree(tree, tree.getRootNode());
    }

    public long addTree(RootedTree tree, Node node) {
        if (!tree.isExternal(node)) {
            for (Node child : tree.getChildren(node)) {
                long l = addTree(tree, node);

            }
        }
        return 0;
    }

    class Lineage {
        double x1, y1, x2, y2;
        long child1, child2;
        boolean isRotated;
        Paint color;
    }
}
