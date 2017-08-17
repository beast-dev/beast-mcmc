package dr.evolution.tree;

/**
 * TestTreeUtils
 *
 * @author Andrew Rambaut
 */

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.treemetrics.*;
import dr.evolution.util.Taxon;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

public class TestTreeUtils {

    public static void main(String[] args) {

        try {

            //5-taxa example
            NewickImporter importer = new NewickImporter(
                    "(((('A':0.6,'B':0.6):0.1,'C':0.5):0.4,'D':0.7):0.1,'E':1.3)");
            Tree tree = importer.importNextTree();
            System.out.println("5-taxon tree: " + tree);


            System.out.println("Tree length: " + TreeUtils.getTreeLength(tree, tree.getRoot()));
            System.out.println();

            Set<Taxon> taxonSet1 = new HashSet<Taxon>();
            taxonSet1.add(tree.getTaxon(0));
            taxonSet1.add(tree.getTaxon(1));
            taxonSet1.add(tree.getTaxon(2));

            System.out.println("Subtree taxa: " + taxonSet1.toString());
            System.out.println("Subtree length: " + TreeUtils.getSubTreeLength(tree, taxonSet1));

            Set<Taxon> taxonSet2 = new HashSet<Taxon>();
            taxonSet2.add(tree.getTaxon(0));
            taxonSet2.add(tree.getTaxon(1));
            taxonSet2.add(tree.getTaxon(3));

            System.out.println("Subtree taxa: " + taxonSet2.toString());
            System.out.println("Subtree length: " + TreeUtils.getSubTreeLength(tree, taxonSet2));

            Set<Taxon> taxonSet3 = new HashSet<Taxon>();
            taxonSet3.add(tree.getTaxon(0));
            taxonSet3.add(tree.getTaxon(2));
            taxonSet3.add(tree.getTaxon(3));
            taxonSet3.add(tree.getTaxon(4));

            System.out.println("Subtree taxa: " + taxonSet3.toString());
            System.out.println("Subtree length: " + TreeUtils.getSubTreeLength(tree, taxonSet3));



        } catch (Importer.ImportException ie) {
            System.err.println(ie);
        } catch (IOException ioe) {
            System.err.println(ioe);
        }

    }

}
