package dr.app.tools;

import dr.app.util.Arguments;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;

import java.io.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.function.Consumer;

class BaseTreeTool {

    static final String NAME_CONTENT = "nameContent";

    // Messages to stderr, output to stdout
    static PrintStream progressStream = System.err;

    int totalTrees = 0;
    private int totalUsedTrees = 0;

    public static class TreeProgressPrinter {
        private final PrintStream stream;
        private static long stepSize = 10000 / 60;

        public TreeProgressPrinter(PrintStream stream) {
            this.stream = stream;
        }

        public void printReadingTrees() {
            stream.println("Reading trees (bar assumes 10,000 trees)...");
            stream.println("0              25             50             75            100");
            stream.println("|--------------|--------------|--------------|--------------|");
        }

        public void printProgress(int totalTrees) {
            if (totalTrees > 0 && totalTrees % stepSize == 0) {
                progressStream.print("*");
                progressStream.flush();
            }
        }

        public void printSummary(int totalTrees, int totalUsedTrees, int burnIn) {
            stream.println();
            stream.println();

            if (totalTrees < 1) {
                System.err.println("No trees");
                return;
            }
            if (totalUsedTrees < 1) {
                System.err.println("No trees past burn-in (=" + burnIn + ")");
                return;
            }

            stream.println("Total trees read: " + totalTrees);
            stream.println("Total trees used: " + totalUsedTrees);
        }

    }

    void readTrees(List<Tree> trees, String inputFileName, int burnIn) throws IOException {

        TreeProgressPrinter progressPrinter = new TreeProgressPrinter(progressStream);
        progressPrinter.printReadingTrees();


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

                progressPrinter.printProgress(totalTrees);

                totalTrees++;
                if (totalTrees > burnIn) {
                    totalUsedTrees++;
                }
            }

        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            return;
        }

        fileReader.close();

        progressPrinter.printSummary(totalTrees, totalUsedTrees, burnIn);

    }

    protected PrintStream openOutputFile(String outputFileName) {
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

    void addTaxonByName(Tree tree, Collection<Taxon> taxa, String name) {
        int taxonId = tree.getTaxonIndex(name);
        if (taxonId == -1) {
            throw new RuntimeException("Unable to find taxon '" + name + "'.");
        }
        taxa.add(tree.getTaxon(taxonId));
    }

    private void addTaxaByNameContent(Tree tree, Collection<Taxon> taxa, String content) {
        int counter = 0;
        for (int i = 0; i < tree.getTaxonCount(); i++) {
            Taxon taxon = tree.getTaxon(i);
            String taxonName = taxon.toString();
            if (taxonName.contains(content)) {
                taxa.add(taxon);
                counter++;
            }
        }
        if (counter == 0) {
            throw new RuntimeException("Unable to find taxon with a name containing '" + content + "'.");
        }
    }

    void getTaxaToFromName(Tree tree, Collection<Taxon> taxa, String[] names, boolean basedOnContent) {
        if (names != null) {
            for (String name : names) {
                if (!basedOnContent) {
                    addTaxonByName(tree, taxa, name);
                } else {
                    addTaxaByNameContent(tree, taxa, name);
                }
            }
        }
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

    static boolean parseBasedOnNameContent(Arguments arguments) {
        String nameContentString = arguments.getStringOption(NAME_CONTENT);
        return nameContentString != null && nameContentString.compareToIgnoreCase("true") == 0;
    }

    protected static void centreLine(String line, int pageWidth) {
        centreLine(line, pageWidth, progressStream);
    }

    public static void centreLine(String line, int pageWidth, PrintStream ps) {
        int n = pageWidth - line.length();
        int n1 = n / 2;
        for (int i = 0; i < n1; i++) {
            ps.print(" ");
        }
        ps.println(line);
    }

    static String[] getInputOutputFileNames(Arguments arguments, Consumer<Arguments> usage) {

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
