package dr.app.treespace;

import jebl.evolution.graphs.Node;
import jebl.evolution.trees.RootedTree;

import java.awt.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
            addTree(tree);
        }
    }

    public double getMaxWidth() {
        return maxWidth;
    }

    public double getMaxHeight() {
        return maxHeight;
    }

    public List<Lineage> getRootLineages() {
        return rootLineages;
    }

    public void addTree(RootedTree tree) {
        // set the tip count to zero for this traversal
        currentY = 0.0;

        // the offset is the distance to the right most tip - used to align the tips of the tree
        offsetX = 0.0;
        Lineage lineage = new Lineage();
        addNode(tree, tree.getRootNode(), lineage, 0.0);


        lineage.dx = -offsetX;

        if (offsetX > maxWidth) {
            maxWidth = offsetX;
        }
        if (currentY > maxHeight) {
            maxHeight = currentY;
        }

        rootLineages.add(lineage);
    }

    public double addNode(RootedTree tree, Node node, Lineage lineage, double cumulativeX) {
        lineage.dx = tree.getLength(node);
        cumulativeX += lineage.dx;

        if (!tree.isExternal(node)) {
            List<Node> children = tree.getChildren(node);
            if (children.size() != 2) {
                throw new RuntimeException("Tree is not binary");
            }

            lineage.child1 = new Lineage();
            lineage.child2 = new Lineage();

            lineage.dy = addNode(tree, children.get(0), lineage.child1, cumulativeX);
            lineage.dy += addNode(tree, children.get(1), lineage.child2, cumulativeX);

            // the y of this node is the average of the two children
            lineage.dy /= 2;

            // now change the children to relative y positions.
            lineage.child1.dy = lineage.child1.dy - lineage.dy;
            lineage.child2.dy = lineage.child2.dy - lineage.dy;

//            if (lineage.child1.tipCount > lineage.child2.tipCount ||
//                    (lineage.child1.tipCount == lineage.child2.tipCount &&
//                            lineage.child1.tipNumber > lineage.child2.tipNumber)) {
//                Lineage tmp = lineage.child1;
//                lineage.child1 = lineage.child2;
//                lineage.child2 = tmp;
//
//                lineage.child1.dy = -lineage.child1.dy;
//                lineage.child2.dy = -lineage.child2.dy;
//            }

            lineage.tipCount = lineage.child1.tipCount + lineage.child2.tipCount;

        } else {
            Integer tipNumber = tipNumbers.get(tree.getTaxon(node).getName());
            if (tipNumber == null) {
                tipNumber = tipNumbers.size();
                tipNumbers.put(tree.getTaxon(node).getName(), tipNumber);
            }

            // the initial (absolute) y position of a tip is its count in the traversal
            lineage.tipNumber = tipNumber;
            lineage.tipCount = 1;
            lineage.dy = tipNumber;
//            lineage.dy = currentY;
            currentY += 1.0;

            if (cumulativeX > offsetX) {
                offsetX = cumulativeX;
            }
        }

        return lineage.dy;
    }

    class Lineage {
        double dx = 0;
        double dy = 0;
        Lineage child1 = null;
        Lineage child2 = null;
        int tipNumber = 0;
        int tipCount = 0;
    }

    private List<Lineage> rootLineages = new ArrayList<Lineage>();
    private Map<String, Integer> tipNumbers = new HashMap<String, Integer>();
    private double maxWidth;
    private double maxHeight;

    private double offsetX;
    private double currentY;
}
