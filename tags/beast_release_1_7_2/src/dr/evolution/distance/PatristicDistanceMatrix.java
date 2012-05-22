package dr.evolution.distance;

import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.Importer;
import dr.matrix.Matrix;

import java.util.HashSet;
import java.util.Set;
import java.io.*;


/**
 * @author Alexei Drummond
 */
public class PatristicDistanceMatrix extends DistanceMatrix {

    Tree tree;

    public PatristicDistanceMatrix(Tree tree) {

        this.tree = tree;

        dimension = tree.getExternalNodeCount();
        distancesKnown = false;

        calculateDistances();
    }

    protected double calculatePairwiseDistance(final int taxon1, int taxon2) {

        NodeRef node1 = tree.getExternalNode(taxon1);
        NodeRef node2 = tree.getExternalNode(taxon2);

        Set<String> nodes = new HashSet<String>();
        nodes.add(tree.getNodeTaxon(node1).getId());
        nodes.add(tree.getNodeTaxon(node2).getId());

        NodeRef ancestor = Tree.Utils.getCommonAncestorNode(tree, nodes);

        if (ancestor == null) System.out.println("common ancestor is null!");

        double height1 = height(tree, node1, ancestor);
        double height2 = height(tree, node2, ancestor);

        return height1 + height2;

    }

    public double height(Tree tree, NodeRef node, NodeRef ancestor) {

        if (node == ancestor) return 0.0;

        return tree.getBranchLength(node) + height(tree, tree.getParent(node), ancestor);

    }

    public String toString() {
        try {
            double[] dists = getUpperTriangle();

            StringBuffer buffer = new StringBuffer();

            buffer.append(" \t");
            for (int k = 0; k < dimension; k++) {
                buffer.append(tree.getNodeTaxon(tree.getExternalNode(k))).append("\t");
            }

            buffer.append("\n");
            int i = 0;
            for (int k = 0; k < dimension; k++) {
                buffer.append(tree.getNodeTaxon(tree.getExternalNode(k)));
                for (int j = 0; j <= k; j++) {
                    buffer.append("\t");
                }
                for (int j = 0; j < dimension-k-1; j++) {
                    buffer.append(String.valueOf(dists[i])).append("\t");
                    i += 1;
                }
                buffer.append("\n");
            }
            return buffer.toString();
        } catch(Matrix.NotSquareException e) {
            return e.toString();
        }

    }


    public static void main(String[] args) throws IOException, Importer.ImportException {

        File file = new File(args[0]);

        File outFile = new File(args[1]);

        PrintWriter writer = new PrintWriter(new FileWriter(outFile));

        NexusImporter importer = new NexusImporter(new FileReader(file));

        Tree[] trees = importer.importTrees(null);

        System.out.println("Found " + trees.length + " trees.");

        for (Tree tree : trees) {

            PatristicDistanceMatrix matrix = new PatristicDistanceMatrix(tree);

            writer.println(matrix.toString());

        }
        writer.flush();

        writer.close();

    }

}
