package dr.app.pathogen;

import java.util.*;

import dr.stats.DiscreteStatistics;
import dr.math.UnivariateFunction;
import dr.math.UnivariateMinimum;
import dr.evolution.tree.*;
import dr.evolution.util.Taxon;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class TemporalStress {

    public static Set<Taxon> annotateStress(MutableTree tree, NodeRef node) {
        Set<Taxon> taxa = new HashSet<Taxon>();

        if (!tree.isExternal(node)) {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                taxa.addAll(annotateStress(tree, child));
            }

            if (taxa.size() > 2) {
                Tree subtree = null;

                double stress = findGlobalRoot(subtree);
                tree.setNodeAttribute(node, "stress", stress);

            }
        } else {
            taxa.add(tree.getNodeTaxon(node));
        }

        return taxa;
    }


    private static double findGlobalRoot(Tree source) {

        FlexibleTree bestTree = new FlexibleTree(source);
        double minF = findLocalRoot(bestTree);

        for (int i = 0; i < source.getNodeCount(); i++) {
            FlexibleTree tmpTree = new FlexibleTree(source);
            NodeRef node = tmpTree.getNode(i);
            if (!tmpTree.isRoot(node)) {
                double length = tmpTree.getBranchLength(node);
                tmpTree.changeRoot(node, length * 0.5, length * 0.5);

                double f = findLocalRoot(tmpTree);
                if (f < minF) {
                    minF = f;
                    bestTree = tmpTree;
                }
            }
        }
        return minF;
    }

    private static double findLocalRoot(final FlexibleTree tree) {

        NodeRef node1 = tree.getChild(tree.getRoot(), 0);
        NodeRef node2 = tree.getChild(tree.getRoot(), 1);

        final double length1 = tree.getBranchLength(node1);
        final double length2 = tree.getBranchLength(node2);

        final double sumLength = length1 + length2;

        final Set<NodeRef> tipSet1 = Tree.Utils.getExternalNodes(tree, node1);
        final Set<NodeRef> tipSet2 = Tree.Utils.getExternalNodes(tree, node2);

        final double[] y = new double[tree.getExternalNodeCount()];

        UnivariateFunction f = new UnivariateFunction() {
            public double evaluate(double argument) {
                double l1 = argument * sumLength;

                for (NodeRef tip : tipSet1) {
                    y[tip.getNumber()] = getRootToTipDistance(tree, tip) - length1 + l1;
                }

                double l2 = (1.0 - argument) * sumLength;

                for (NodeRef tip : tipSet2) {
                    y[tip.getNumber()] = getRootToTipDistance(tree, tip) - length2 + l2;
                }

                return DiscreteStatistics.variance(y);
            }

            public double getLowerBound() { return 0; }
            public double getUpperBound() { return 1.0; }
        };

        UnivariateMinimum minimum = new UnivariateMinimum();

        double x = minimum.findMinimum(f);

        double fminx = minimum.fminx;

        double l1 = x * sumLength;
        double l2 = (1.0 - x) * sumLength;

        tree.setBranchLength(node1, l1);
        tree.setBranchLength(node2, l2);

        return fminx;
    }

    private static double getRootToTipDistance(Tree tree, NodeRef node) {
        double distance = 0;
        while (node != null) {
            distance += tree.getBranchLength(node);
            node = tree.getParent(node);
        }
        return distance;
    }


}
