package dr.app.tools;

import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.Tree;

import java.io.FileReader;
import java.io.IOException;


public class SequentialTreeReader {

    private static final Boolean DEBUG = false;

    private final String inputFileName;
    private final int burnin;
    private static BaseTreeTool.TreeProgressPrinter progressPrinter = new BaseTreeTool.TreeProgressPrinter(System.err);
    private int totalTrees;
    private int totalUsedTrees;

    private Tree currentTree;
    private int currentIndex;

    private final FileReader reader;
    private final TreeImporter importer;


    SequentialTreeReader(String inputFileName, int burnin) throws IOException {
        this.inputFileName = inputFileName;
        this.burnin = burnin;
        this.reader = new FileReader(inputFileName);
        this.importer = new NexusImporter(reader, false);

    }


    public Tree getTree(int index) throws IOException {

        if (DEBUG) {
            System.out.println(index);
        }

        if (index < burnin) {
            throw new RuntimeException("Tree " + index + " was requested, but burnin is set at " + burnin + ".");
        }

        try {
            if (currentTree == null) {
                progressPrinter.printReadingTrees();
                currentTree = importer.importNextTree();
                currentIndex = 0;
                totalTrees = 1;
                totalUsedTrees = 1;
            }

            if (index < currentIndex) {
                throw new RuntimeException("The last tree accessed was " + currentIndex + ", but " + index +
                        " was requested. Cannot go backwards. (for developers: consider using BaseTreeTools.readTrees)");
            } else if (index == currentIndex) {
                return currentTree;
            } else {
                for (int i = currentIndex; i < index; i++) {


                    if (importer.hasTree()) {
                        progressPrinter.printProgress(totalTrees);
                        currentTree = importer.importNextTree();
                        totalTrees++;
                    } else {
                        reader.close();
                        progressPrinter.printSummary(totalTrees, totalUsedTrees, burnin);
                        return null;
                    }
                }

                totalUsedTrees++;
                currentIndex = index;
                return currentTree;
            }

        } catch (Importer.ImportException e) {
            reader.close();
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            return null;
        }
    }
}
