package dr.app.tools;

import dr.app.beast.BeastVersion;
import dr.app.util.Arguments;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.util.Version;
import dr.xml.XMLParseException;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */

public class TreeLogToEmpiricalDistributionConverter extends BaseTreeTool {

    private final static Version version = new BeastVersion();
    private static final String BURN_IN = "burnIn";

    private TreeLogToEmpiricalDistributionConverter(String inputFileName,
                                                    String outputFileName,
                                                    int burnIn,
                                                    double startTime,
                                                    double endTime) throws IOException {

        List<Tree> trees = new ArrayList<>();

        readTrees(trees, inputFileName, burnIn);

        this.startTime = startTime;
        this.endTime = endTime;

        this.ps = openOutputFile(outputFileName);
        processTrees(trees, burnIn);
        closeOutputFile(ps);
    }

    private Tree[] getEmpiricalDistribution() throws XMLParseException {
        Tree[] trees = null;
//        NexusImporter importer = null;
//        try {
//            FileReader reader = new FileReader(file);
//            importer = new NexusImporter(reader);
//
//                trees = importer.importTrees(taxa, true); // Re-order taxon numbers to original TaxonList order
//                reader.close();
//
//        } catch (FileNotFoundException e) {
//            throw new XMLParseException(e.getMessage());
//        } catch (IOException e) {
//            throw new XMLParseException(e.getMessage());
//        } catch (Importer.ImportException e) {
//            throw new XMLParseException(e.getMessage());
//        }
        return trees;
    }

    private void processTrees(List<Tree> trees, int burnIn) {
        if (burnIn < 0) {
            burnIn = 0;
        }
        for (int i = burnIn; i < trees.size(); ++i) {
            Tree tree = trees.get(i);
            processOneTree(tree);
        }
    }

    private void processOneTree(Tree tree) {

        String treeId = tree.getId();
        if (treeId.startsWith("STATE_")) {
            treeId = treeId.replaceFirst("STATE_", "");
        }

        int count = 0;

        for (int i = 0; i < tree.getInternalNodeCount(); ++i) {
            NodeRef node = tree.getInternalNode(i);
            double nodeHeight = tree.getNodeHeight(node);
            if (nodeHeight >= startTime && nodeHeight <= endTime) {
                ++count;
            }
        }

        ps.println(treeId + "\t" + count);
    }

    private void closeOutputFile(PrintStream ps) {
        if (ps != null) {
            ps.close();
        }
    }

    private double startTime;
    private double endTime;
    private PrintStream ps;

    public static void printTitle() {
        progressStream.println();
        centreLine("CoalescentEventSummarizer " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("tool to summarize the number of coalescent events in a particular time interval", 60);
        centreLine("by", 60);
        centreLine("Marc A. Suchard and Joel Wertheim", 60);
        progressStream.println();
        progressStream.println();
    }

    public static void printUsage(Arguments arguments) {
        arguments.printUsage("CoalescentEventSummarizer", "<input-file-name> [<output-file-name>]");
        progressStream.println();
        progressStream.println("  Example: CoalescentEventSummarizer -startTime 10.0 -endTime 2.0 input.trees output.log");
        progressStream.println();
    }

    public static void main(String[] args) throws IOException {

        int burnIn = -1;
        double startTime = 0; //default start time considered to be time zero
        double endTime = Double.MAX_VALUE; //end time, backwards in time

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption(BURN_IN, "the number of states to be considered as 'burn-in' [default = 0]"),
                        new Arguments.RealOption("startTime", "The start time for to time summaries [default = 0] "),
                        new Arguments.RealOption("endTime", "The end time for to time summaries [default = Double.MAX_VALUE] "),
                        new Arguments.Option("help", "option to print this message")
                });

        handleHelp(arguments, args, TaxaMarkovJumpHistoryAnalyzer::printUsage);

        if (arguments.hasOption("startTime")) {
            startTime = arguments.getRealOption("startTime");
        }

        if (arguments.hasOption("endTime")) {
            endTime = arguments.getRealOption("endTime");
        }

        if (arguments.hasOption(BURN_IN)) {
            burnIn = arguments.getIntegerOption(BURN_IN);
            System.err.println("Ignoring a burn-in of " + burnIn + " trees.");
        }

        String[] fileNames = getInputOutputFileNames(arguments, TreeLogToEmpiricalDistributionConverter::printUsage);

        new TreeLogToEmpiricalDistributionConverter(fileNames[0], fileNames[1], burnIn,
                startTime,
                endTime
        );
        System.exit(0);
    }
}
