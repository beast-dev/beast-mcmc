package dr.app.tools;

import dr.app.util.Arguments;
import dr.app.beast.BeastVersion;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.util.Version;


import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class TreeStateTimeSummarizer extends BaseTreeTool {

    private final static Version version = new BeastVersion();
    private static final String BURN_IN = "burnIn";
    private static final String HISTORY = "history";

    private TreeStateTimeSummarizer(String inputFileName,
                                    String outputFileName,
                                    int burnIn,
                                    double startTime,
                                    double endTime,
                                    String[] states
//                                  String nodeStateAnnotation
    ) throws IOException {

        List<Tree> trees = new ArrayList<>();

        readTrees(trees, inputFileName, burnIn);

        this.startTime = startTime;
        this.endTime = endTime;

        this.ps = openOutputFile(outputFileName);
        processTrees(trees, burnIn, states);
        closeOutputFile(ps);
    }

    private void processTrees(List<Tree> trees, int burnIn, String[] states) {
        if (burnIn < 0) {
            burnIn = 0;
        }
        for (int i = burnIn; i < trees.size(); ++i) {
            Tree tree = trees.get(i);
            //TODO: your magic
            //processOneTree(tree, states);
        }
    }

    private static Object[] readMJH(NodeRef node, Tree treeTime) {
        if (treeTime.getNodeAttribute(node, HISTORY) != null) {
            return (Object[]) treeTime.getNodeAttribute(node, HISTORY);
        } else {
            return null;
        }
    }

    private class Row {
        String treeId;
        String state;
        double time;

        private static final String DELIMITER = ",";

        private Row(String treeId,
                    String state,
                    double time) {
            this.treeId = treeId;
            this.state = state;
            this.time = time;
        }

        public String toString() {
            return treeId + DELIMITER
                    + state + DELIMITER
                    + time + DELIMITER;
        }
    }

    protected PrintStream openOutputFile(String outputFileName) {

        PrintStream ps = super.openOutputFile(outputFileName);
        ps.println("treeId,state,time");
        return ps;
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
        centreLine("TreeStateTimeSummarizer " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("tool to summarize the times spent in particular discrete states according to Markov jump annotations", 60);
        centreLine("by", 60);
        centreLine("Philippe Lemey and Marc Suchard", 60);
        progressStream.println();
        progressStream.println();
    }

    public static void printUsage(Arguments arguments) {
        arguments.printUsage("TreeStateTimeSummarizer", "<input-file-name> [<output-file-name>]");
        progressStream.println();
        progressStream.println("  Example: TreeStateTimeSummarizer -taxaToProcess taxon1,taxon2 input.trees output.trees");
        progressStream.println();
    }

    //Main method
    public static void main(String[] args) throws IOException {

        int burnIn = -1;
        double startTime = 0; //default start time considered to be time zero
        double endTime = Double.MAX_VALUE; //end time, backwards in time
        String[] states = null;
        //Perhaps not needed if we rely (entirely) on the Markov jump history?
        //String nodeStateAnnotation = null;

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption(BURN_IN, "the number of states to be considered as 'burn-in' [default = 0]"),
                        new Arguments.RealOption("startTime", "The start time for to time summaries [default = 0] "),
                        new Arguments.RealOption("endTime", "The end time for to time summaries [default = Double.MAX_VALUE] "),
                        new Arguments.StringOption("states", "String", "a comma-separated list of discrete states for which times need to be summarized"),
                        //new Arguments.StringOption("nodeStateAnnotation", "String", "string used for node state annotations"),
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

        if (arguments.hasOption("states")) {
            states = Branch2dRateToGrid.parseVariableLengthStringArray(arguments.getStringOption("states"));
        } else {
            System.err.print("no states provided for time summaries... nothing to summarize");
            System.exit(-1);
        }

//        if (arguments.hasOption("nodeStateAnnotation")) {
//            nodeStateAnnotation = arguments.getStringOption("nodeStateAnnotation");
//       }

        String[] fileNames = getInputOutputFileNames(arguments, TreeStateTimeSummarizer::printUsage);

        new TreeStateTimeSummarizer(fileNames[0], fileNames[1], burnIn,
                startTime,
                endTime,
                states
//                nodeStateAnnotation
        );
        System.exit(0);
    }

}
