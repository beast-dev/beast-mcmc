package dr.app.tools;

import dr.app.beauti.types.SequenceErrorType;
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
                                    String nodeStateAnnotation) throws IOException {

        List<Tree> trees = new ArrayList<>();

        readTrees(trees, inputFileName, burnIn);

        this.startTime = startTime;
        this.endTime = endTime;
        this.nodeStateAnnotation = nodeStateAnnotation;

        this.ps = openOutputFile(outputFileName);
        processTrees(trees, burnIn);
        closeOutputFile(ps);
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

        NodeRef root = tree.getRoot();

        String rootState = (String) tree.getNodeAttribute(root, nodeStateAnnotation);
        if(rootState==null){
            System.err.println("no annotation for "+nodeStateAnnotation);
            System.exit(-1);
        }

        StateHistory history = traversePostOrder(tree, root, treeId);

        checkEqual(rootState, history.state);
        Row row = new Row(treeId, history.state, history.duration);
        ps.println(row);

    }

    private class StateHistory {
        String state;
        double duration;

        StateHistory(String state, double duration) {
            this.state = state;
            this.duration = duration;
        }
    }

    private void checkEqual(String lhs, String rhs) {
        if (!lhs.equals(rhs)) {
            System.err.println("State mismatch");
            System.exit(-1);
        }
    }

    private StateHistory traversePostOrder(Tree tree, NodeRef node, String treeId) {

        StateHistory history;

        // If the node is internal, update the partial likelihoods.
        if (tree.isExternal(node)) {
            history = new StateHistory((String) tree.getNodeAttribute(node, nodeStateAnnotation), 0.0);
        } else {
            StateHistory history0 = traversePostOrder(tree, tree.getChild(node, 0), treeId);
            StateHistory history1 = traversePostOrder(tree, tree.getChild(node, 1), treeId);

            checkEqual(history0.state, history1.state);

            history = new StateHistory(history0.state, history0.duration + history1.duration);
        }

        if (!tree.isRoot(node)) {

            double currentTime = tree.getNodeHeight(node);

            Object[] jumps = readMJH(node, tree);
            if (jumps != null) {
                for (int i = jumps.length - 1; i >= 0; --i) {
                    Object[] jump = (Object[]) jumps[i];

                    checkEqual(history.state, (String) jump[2]);

                    double jumpTime = (Double) jump[0];

                    history.duration += clippedDuration(jumpTime, currentTime); //(jumpTime - currentTime);

                    Row row = new Row(treeId, history.state, history.duration);
                    ps.println(row);

                    currentTime = jumpTime;
                    history = new StateHistory((String) jump[1], 0.0);

                }
            }

            double parentHeight = tree.getNodeHeight(tree.getParent(node));

            history.duration += clippedDuration(parentHeight, currentTime); //parentHeight - currentTime;

        }

        return history;
    }


    private double clippedDuration(double end, double start) {

        if (start > endTime || end < startTime) {
            return 0.0;
        } else {
            double clippedEnd = Math.min(endTime, end);
            double clippedStart = Math.max(startTime, start);

            return clippedEnd - clippedStart;
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
                    + time;
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
    private String nodeStateAnnotation;

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
        String nodeStateAnnotation = "states";

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption(BURN_IN, "the number of states to be considered as 'burn-in' [default = 0]"),
                        new Arguments.RealOption("startTime", "The start time for to time summaries [default = 0] "),
                        new Arguments.RealOption("endTime", "The end time for to time summaries [default = Double.MAX_VALUE] "),
                        new Arguments.StringOption("nodeStateAnnotation", "String", "string used for node state annotations [default = states] "),
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

        if (arguments.hasOption("nodeStateAnnotation")) {
            nodeStateAnnotation = arguments.getStringOption("nodeStateAnnotation");
        } else {
            System.err.print("no state annotation string provided... trying with 'states...\n");
        }

        String[] fileNames = getInputOutputFileNames(arguments, TreeStateTimeSummarizer::printUsage);

        new TreeStateTimeSummarizer(fileNames[0], fileNames[1], burnIn,
                startTime,
                endTime,
                nodeStateAnnotation
        );
        System.exit(0);
    }

}
