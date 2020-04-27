package dr.app.tools;

import dr.app.util.Arguments;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

class BaseTreeTool {

    // Messages to stderr, output to stdout
    static PrintStream progressStream = System.err;

    int totalTrees = 0;
    private int totalUsedTrees = 0;

    void readTrees(List<Tree> trees, String inputFileName, int burnin) throws IOException {

        progressStream.println("Reading trees (bar assumes 10,000 trees)...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        long stepSize = 10000 / 60;

        FileReader fileReader = new FileReader(inputFileName);
        TreeImporter importer = new NexusImporter(fileReader, false);

        try {
            totalTrees = 0;
            while (importer.hasTree()) {

                Tree tree = importer.importNextTree();
                if (trees == null) {
                    trees = new ArrayList<>();
                }
                trees.add(tree);

                if (totalTrees > 0 && totalTrees % stepSize == 0) {
                    progressStream.print("*");
                    progressStream.flush();
                }
                totalTrees++;
                if (totalTrees > burnin) {
                    totalUsedTrees++;
                }
            }

        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            return;
        }

        fileReader.close();

        progressStream.println();
        progressStream.println();

        if (totalTrees < 1) {
            System.err.println("No trees");
            return;
        }
        if (totalUsedTrees < 1) {
            System.err.println("No trees past burnin (=" + burnin + ")");
            return;
        }

        progressStream.println("Total trees read: " + totalTrees);
        progressStream.println("Total trees used: " + totalUsedTrees);
    }

    PrintStream openOutputFile(String outputFileName) {
        PrintStream ps = null;

        if (outputFileName == null) {
            ps = progressStream;
        } else {
            try {
                ps = new PrintStream(new File(outputFileName));
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }

        return ps;
    }

    void addTaxonByName(Tree tree, List<Taxon> taxa, String name) {
        int taxonId = tree.getTaxonIndex(name);
        if (taxonId == -1) {
            throw new RuntimeException("Unable to find taxon '" + name + "'.");
        }
        taxa.add(tree.getTaxon(taxonId));
    }

    static void handleHelp(Arguments arguments, String[] args, Consumer<Arguments> printUsage) {

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            progressStream.println(ae);
            printUsage.accept(arguments);
            System.exit(1);
        }

        if (arguments.hasOption("help")) {
            printUsage.accept(arguments);
            System.exit(0);
        }
    }

    static void centreLine(String line, int pageWidth) {
        int n = pageWidth - line.length();
        int n1 = n / 2;
        for (int i = 0; i < n1; i++) {
            progressStream.print(" ");
        }
        progressStream.println(line);
    }

    static String[] getInputOutputFileNames(Arguments arguments, String[] args, Consumer<Arguments> usage) {

        String[] result = new String[2];

        final String[] args2 = arguments.getLeftoverArguments();

        switch (args2.length) {
            case 2:
                result[1] = args2[1];
                // fall to
            case 1:
                result[0] = args2[0];
                break;
            default: {
                System.err.println("Unknown option: " + args2[2]);
                System.err.println();
                usage.accept(arguments);
                System.exit(1);
            }
        }

        return result;
    }
}
