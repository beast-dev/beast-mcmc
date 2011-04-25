package dr.app.treespace;

import jebl.evolution.graphs.Node;
import jebl.evolution.trees.RootedTree;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class TreeLineages {

    private final static String[] LOCATIONS = {
            "Africa",
            "USA",
            "Taiwan",
            "China",
            "Russia",
            "Oceania",
            "Asia",
            "Japan",
            "Mexico",
            "South America",
            "Canada",
            "Europe",
            "Southeast Asia",
            "South Korea"
    };

    private final static String[] AIR_COMMUNITIES = {
            "AC1",
            "AC2",
            "AC3",
            "AC4",
            "AC5",
            "AC6",
            "AC7",
            "AC8",
            "AC9",
            "AC10",
            "AC11",
            "AC12",
            "AC13",
            "AC14"
    };

    Map<String, Integer> locationMap = new HashMap<String, Integer>();

    public TreeLineages() {
        for (int i = 0; i < AIR_COMMUNITIES.length; i++) {
            locationMap.put(AIR_COMMUNITIES[i], i);
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
        positionNode(lineage);

        lineage.dx = -offsetX;

        if (offsetX > maxWidth) {
            maxWidth = offsetX;
        }
        if (currentY > maxHeight) {
            maxHeight = currentY;
        }

        rootLineages.add(lineage);
    }

    public void addNode(RootedTree tree, Node node, Lineage lineage, double cumulativeX) {
        lineage.dx = tree.getLength(node);
        cumulativeX += lineage.dx;

        String location = (String)node.getAttribute("states");
        if (location != null) {
            lineage.state = locationMap.get(location);
        }

        if (!tree.isExternal(node)) {
            List<Node> children = tree.getChildren(node);
            if (children.size() != 2) {
                throw new RuntimeException("Tree is not binary");
            }

            lineage.child1 = new Lineage();
            lineage.child2 = new Lineage();

            addNode(tree, children.get(0), lineage.child1, cumulativeX);
            addNode(tree, children.get(1), lineage.child2, cumulativeX);

            lineage.tipCount = lineage.child1.tipCount + lineage.child2.tipCount;

            if (lineage.child1.tipCount > lineage.child2.tipCount ||
                    (lineage.child1.tipCount == lineage.child2.tipCount &&
                            lineage.child1.tipNumber > lineage.child2.tipNumber)) {
                Lineage tmp = lineage.child1;
                lineage.child1 = lineage.child2;
                lineage.child2 = tmp;
            }
        } else {
            Integer tipNumber = taxonNumbers.get(tree.getTaxon(node).getName());
            if (tipNumber == null) {
                tipNumber = taxonNumbers.size();
                taxonNumbers.put(tree.getTaxon(node).getName(), tipNumber);
            }

            lineage.tipNumber = tipNumber;
            lineage.tipCount = 1;

            if (cumulativeX > offsetX) {
                offsetX = cumulativeX;
            }
        }
    }

    public double positionNode(Lineage lineage) {
        if (lineage.child1 != null) {
            lineage.dy = positionNode(lineage.child1);
            lineage.dy += positionNode(lineage.child2);

            // the y of this node is the average of the two children
            lineage.dy /= 2.0;
//            lineage.dy = ((double)lineage.tipCount) / 2.0;

            // now change the children to relative y positions.
            lineage.child1.dy = lineage.child1.dy - lineage.dy;
            lineage.child2.dy = lineage.child2.dy - lineage.dy;

        } else {
            Integer orderedTipNumber = orderedTipNumbers.get(lineage.tipNumber);
            if (orderedTipNumber == null) {
                orderedTipNumber = orderedTipNumbers.size();
                orderedTipNumbers.put(lineage.tipNumber, orderedTipNumber);
            }

            lineage.tipNumber = orderedTipNumber;
            lineage.dy = lineage.tipNumber;
//            lineage.dy = currentY;
            currentY += 1.0;
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
        int state;
    }

    private List<Lineage> rootLineages = new ArrayList<Lineage>();
    private Map<String, Integer> taxonNumbers = new HashMap<String, Integer>();
    private Map<Integer, Integer> orderedTipNumbers = new HashMap<Integer, Integer>();
    private double maxWidth;
    private double maxHeight;

    private double offsetX;
    private double currentY;
}
