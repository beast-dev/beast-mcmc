package dr.app.tools;

        import dr.app.beast.BeastVersion;
        import dr.app.util.Arguments;
        import dr.evolution.tree.NodeRef;
        import dr.evolution.tree.Tree;
        import dr.evolution.tree.TreeUtils;
        import dr.util.Version;

        import java.io.*;
        import java.util.*;

/**
 * @author Philippe Lemey
 * @author Marc Suchard
 */
public class PersistenceSummarizer extends BaseTreeTool {

    private final static Version version = new BeastVersion();

    private static final String HISTORY = "history";
    private static final String BURN_IN = "burnIn";
    private static final String[] falseTrue = {"false", "true"};

//    private static final boolean NEW_OUTPUT = true;

    private PersistenceSummarizer(String inputFileName,
                                  String outputFileName,
                                  int burnIn,
                                  double[] evaluationTimes,
                                  double[] ancestralTimes,
//                                         double mrsd,
                                  String nodeStateAnnotation
    ) throws IOException {

//        List<Tree> trees = new ArrayList<>();

//        readTrees(trees, inputFileName, burnIn);

        SequentialTreeReader treeReader = new SequentialTreeReader(inputFileName, burnIn);
        //       this.mrsd = mrsd;
        //        this.evaluationTime = evaluationTime;
        //        this.ancestryTime = ancestryTime;

        this.ps = openOutputFile(outputFileName);
        //processTrees(trees, burnIn, evaluationTime, ancestryTime, nodeStateAnnotation);
        processTrees(treeReader, burnIn, evaluationTimes, ancestralTimes, nodeStateAnnotation);
        closeOutputFile(ps);
    }

    private void processTrees(SequentialTreeReader treeReader, int burnIn, double[] evaluationTimes, double[] ancestralTimes, String nodeStateAnnotation) throws IOException {
        if (burnIn < 0) {
            burnIn = 0;
        }

        int index = burnIn;
        Tree tree;

        while (treeReader.getTree(index) != null) {
            tree = treeReader.getTree(index);
            for (int i = 0; i < evaluationTimes.length; i++) {
                processOneTree(tree, evaluationTimes[i], ancestralTimes[i], nodeStateAnnotation);
            }
            index++;
        }
    }

    private void processOneTree(Tree tree, double evaluationTime, double ancestralTime, String nodeStateAnnotation) {

        String treeId = tree.getId();
        if (treeId.startsWith("STATE_")) {
            treeId = treeId.replaceFirst("STATE_", "");
        }

        Set nodes = new HashSet();

        for (int i = 0; i < tree.getNodeCount(); ++i) {
            NodeRef node = tree.getNode(i);

            if (!tree.isRoot(node)){
                double nodeHeight = tree.getNodeHeight(node);
                NodeRef parentNode = tree.getParent(node);
                double parentNodeHeight = tree.getNodeHeight(parentNode);
                if ((nodeHeight <= evaluationTime) && (evaluationTime < parentNodeHeight)) {
//                    System.out.println("node "+i);

                    Set nodeDescendants = TreeUtils.getExternalNodes(tree, node);
                    Set totalEventDescendents = null;
                    NodeRef originalNode = tree.getNode(i);

                    String nodeState = (String) tree.getNodeAttribute(node, nodeStateAnnotation);
                    if (nodeState == null) {
                        throw new RuntimeException("Could not locate node state annotation '" + nodeStateAnnotation +
                                "' for node " + node.getNumber());
                    }

                    String currentState = nodeState;
                    String ancestralState = nodeState;
                    double currentStateTime = 0;
                    boolean jumpOccured = false;
                    Object[] jumps = readCJH(node, tree);

                    //when no jumps are found, we need to go back it's ancestry to find a branch with a jump
                    if (jumps == null) {
                        currentStateTime += parentNodeHeight - evaluationTime;
//                        System.out.println("1\t"+currentStateTime);
                        node = tree.getParent(node);
                        while (!jumpOccured && nodeHeight<ancestralTime) {
                        //while (!jumpOccured) {
                            if (!tree.isRoot(node)){
                                nodeHeight = tree.getNodeHeight(node);
                                parentNodeHeight = tree.getNodeHeight(tree.getParent(node));
                                jumps = readCJH(node, tree);
                                if (jumps == null) {
                                    currentStateTime += parentNodeHeight - nodeHeight;
//                                    System.out.println("2\t"+currentStateTime);
                                    totalEventDescendents = TreeUtils.getExternalNodes(tree, node);
                                } else {
                                    Object[] jump = getMostRecentJump(jumps);
                                    ancestralState = (String) jump[1];
                                    currentStateTime += (Double) jump[0] - nodeHeight;
//                                    System.out.println("3\t"+currentStateTime);
                                    jumpOccured = true;
                                    totalEventDescendents = TreeUtils.getExternalNodes(tree, node);
                                }
                                node = tree.getParent(node);
                            } else {
                                currentStateTime += parentNodeHeight - nodeHeight;
//                                System.out.println("4\t"+currentStateTime);
                                totalEventDescendents = TreeUtils.getExternalNodes(tree, node);
                                jumpOccured = true;
                            }
                        }
                    //jumps
                    } else {

                        double mostRecentJumpTimeLargerThanEvaluationTime = Double.MAX_VALUE;
                        double oldestJumpTimeSmallerThanEvaluationTime = Double.MIN_VALUE;
                        Object[] recentJump = new Object[3];
                        boolean jumpBeforeEvaluationTime = false;
                        for (int j = jumps.length - 1; j >= 0; j--) {
                            Object[] currentJump = (Object[]) jumps[j];
                            double jumpTime = (Double) currentJump[0];
                            if (jumpTime < evaluationTime){
                                if (jumpTime > oldestJumpTimeSmallerThanEvaluationTime){
                                    currentState = (String) currentJump[1];
                                }
                            } else if (jumpTime < mostRecentJumpTimeLargerThanEvaluationTime) {
                                mostRecentJumpTimeLargerThanEvaluationTime = jumpTime;
                                recentJump = currentJump;
                                jumpBeforeEvaluationTime = true;
                            }
                        }
                        if (jumpBeforeEvaluationTime){
                            currentStateTime += (Double) recentJump[0] - evaluationTime;
//                            System.out.println("5\t"+currentStateTime);
                            ancestralState = (String) recentJump[1];
                            totalEventDescendents = TreeUtils.getExternalNodes(tree, node);
                        } else {
                            currentStateTime += parentNodeHeight - evaluationTime;
//                            System.out.println("6\t"+currentStateTime);
                            currentState = (String) tree.getNodeAttribute(tree.getParent(node), nodeStateAnnotation);
                            ancestralState = (String) tree.getNodeAttribute(tree.getParent(node), nodeStateAnnotation);
                            jumpOccured = false;
                            node = tree.getParent(node);
                            while (!jumpOccured && nodeHeight<ancestralTime) {
                            //while (!jumpOccured) {
                                if (!tree.isRoot(node)) {
                                    nodeHeight = tree.getNodeHeight(node);
                                    parentNodeHeight = tree.getNodeHeight(tree.getParent(node));
                                    jumps = readCJH(node, tree);
                                    if (jumps == null) {
                                        currentStateTime += parentNodeHeight - nodeHeight;
//                                        System.out.println("7\t"+currentStateTime);
                                        totalEventDescendents = TreeUtils.getExternalNodes(tree, node);
                                    } else {
                                         //iterate over jumps: get the most recent one
                                        Object[] jump = getMostRecentJump(jumps);
                                        ancestralState = (String) jump[1];
                                        currentStateTime += (Double) jump[0] - nodeHeight;
//                                        System.out.println("8\t"+currentStateTime);
                                        jumpOccured = true;
                                        totalEventDescendents = TreeUtils.getExternalNodes(tree, node);
                                    }
                                    node = tree.getParent(node);
                                } else {
                                    currentStateTime += parentNodeHeight - nodeHeight;
//                                    System.out.println("9\t"+currentStateTime);
                                    totalEventDescendents = TreeUtils.getExternalNodes(tree, node);
                                    jumpOccured = true;
                                }
                            }
                        }
                    }

                    double independenceTime = getIndependenceTime(tree,nodes,originalNode,evaluationTime,currentState,nodeStateAnnotation);
                    nodes.add(originalNode);

                    Row row = new Row(treeId, evaluationTime, ancestralTime, originalNode.getNumber(), node.getNumber(), currentState, ancestralState, currentStateTime, independenceTime,
                            nodeDescendants.size(),getSameStateDescendants(nodeDescendants,tree,currentState,nodeStateAnnotation, 0),
                            totalEventDescendents.size(),getSameStateDescendants(totalEventDescendents,tree,currentState,nodeStateAnnotation, 0),
                            nodesAfterEvalTime(totalEventDescendents, tree, evaluationTime),getSameStateDescendants(totalEventDescendents,tree,currentState,nodeStateAnnotation, evaluationTime));
                    ps.println(row);
                }
            }

        }
   }

    private Object[] getMostRecentJump(Object[] jumps){
            double mostRecentJumpTime = Double.MAX_VALUE;
            Object[] recentJump = new Object[3];
            //get most recent jump
            for (int j = jumps.length - 1; j >= 0; j--) {
                Object[] currentJump = (Object[]) jumps[j];
                double jumpTime = (Double) currentJump[0];
                if (jumpTime < mostRecentJumpTime) {
                    mostRecentJumpTime = jumpTime;
                    recentJump = currentJump;
                }
            }
        return recentJump;
    }

    private double getIndependenceTime(Tree tree, Set nodes, NodeRef node, double evaluationTime, String currentState, String nodeStateAnnotation){
        double independenceTime =  tree.getNodeHeight(tree.getRoot());
        Iterator iter = nodes.iterator();
        while (iter.hasNext()) {
            NodeRef currentNode = (NodeRef) iter.next();
            NodeRef mrcaNode = TreeUtils.getCommonAncestorNode(tree,node,currentNode);
            double mrcaNodeHeight = tree.getNodeHeight(mrcaNode);

            if ((mrcaNodeHeight-evaluationTime)<independenceTime){
                String currentNodeState = (String) tree.getNodeAttribute(currentNode, nodeStateAnnotation);
                //we change the node state if it has changed in between evaluationTime and the node, so to what it should be at evaluationTime
                Object[] jumps = readCJH(currentNode, tree);
                if(jumps!=null){
                    double oldestJumpTimeSmallerThanEvaluationTime = Double.MIN_VALUE;
                    for (int j = jumps.length - 1; j >= 0; j--) {
                        Object[] currentJump = (Object[]) jumps[j];
                        double jumpTime = (Double) currentJump[0];
                        if (jumpTime < evaluationTime){
                            if (jumpTime > oldestJumpTimeSmallerThanEvaluationTime){
                                currentNodeState = (String) currentJump[1];
                            }
                        }
                    }
                }

                if(currentState.equalsIgnoreCase(currentNodeState)){
                    if (stateMaintainedBetweenNodes(currentState, tree, node, mrcaNode, nodeStateAnnotation)){
                        if (stateMaintainedBetweenNodes(currentNodeState, tree, currentNode, mrcaNode, nodeStateAnnotation)){
                            independenceTime = (mrcaNodeHeight-evaluationTime);
                        }
                    }
                }
            }
        }
        return independenceTime;
    }


    private int getSameStateDescendants(Set leafs, Tree tree, String state, String nodeStateAnnotation, double evaluationTime){
        int descendents = 0;
        Iterator iter = leafs.iterator();
        while (iter.hasNext()) {
            NodeRef currentNode = (NodeRef)iter.next();
            String nodeState = (String) tree.getNodeAttribute(currentNode, nodeStateAnnotation);
            if (nodeState.equalsIgnoreCase(state)){
                if(evaluationTime > 0){
                    if (tree.getNodeHeight(currentNode) < evaluationTime){
                        descendents ++;
                    }
                } else{
                    descendents ++;
                }
            }
        }
        return descendents;
    }

    private int nodesAfterEvalTime(Set<NodeRef> tips, Tree tree, double evaluationTime) {
        int count = 0;
        Iterator iter = tips.iterator();
        while (iter.hasNext()) {
            NodeRef currentNode = (NodeRef) iter.next();
            if (tree.getNodeHeight(currentNode) < evaluationTime) {
                count++;
            }
        }
        return count;
    }

    private boolean stateMaintainedBetweenNodes(String state, Tree tree, NodeRef node, NodeRef mrcaNode, String nodeStateAnnotation){
        boolean stateMaintained = true;
        String nodeState = state;
        while (!node.equals(mrcaNode)){
            NodeRef parentNode = tree.getParent(node);
            String parentNodeState = (String) tree.getNodeAttribute(parentNode, nodeStateAnnotation);
            if (!nodeState.equalsIgnoreCase(parentNodeState)){
                stateMaintained = false;
            }
            node = parentNode;
            nodeState = parentNodeState;
        }
        return stateMaintained;
    }

    private class Row {
        String treeId;
        double evaluationTime;
        double ancestralTime;
        int startNodeID;
        int endNodeID;
        String startLocation;
        String endLocation;
        double time;
        double independenceTime;
        int numberOfDescendants;
        int numberOfDescendantsOfSameState;
        int totalNumberOfDescendantsFromUniqueEvent;
        int totalnumberOfDescendantsFromUniqueEventAndSameState;
        int numberOfDescendantsFromUniqueEventAfterEvalTime;
        int numberOfDescendantsFromUniqueEventAndSameStateAfterEvalTime;

        private static final String DELIMITER = ",";

        private Row(String treeId, double evaluationTime, double ancestralTime,
                    int startNodeID, int endNodeID,
                    String startLocation, String endLocation,
                    double time, double independenceTime, int numberOfDescendants, int numberOfDescendantsOfSameState,
                    int totalNumberOfDescendantsFromUniqueEvent, int totalnumberOfDescendantsFromUniqueEventAndSameState,
                    int numberOfDescendantsFromUniqueEventAfterEvalTime, int numberOfDescendantsFromUniqueEventAndSameStateAfterEvalTime
        ) {
            this.treeId = treeId;
            this.evaluationTime = evaluationTime;
            this.startNodeID = startNodeID;
            this.endNodeID = endNodeID;
            this.ancestralTime = ancestralTime;
            this.startLocation = startLocation;
            this.endLocation = endLocation;
            this.time = time;
            this.independenceTime = independenceTime;
            this.numberOfDescendants = numberOfDescendants;
            this.numberOfDescendantsOfSameState = numberOfDescendantsOfSameState;
            this.totalNumberOfDescendantsFromUniqueEvent = totalNumberOfDescendantsFromUniqueEvent;
            this.totalnumberOfDescendantsFromUniqueEventAndSameState = totalnumberOfDescendantsFromUniqueEventAndSameState;
            this.numberOfDescendantsFromUniqueEventAfterEvalTime = numberOfDescendantsFromUniqueEventAfterEvalTime;
            this.numberOfDescendantsFromUniqueEventAndSameStateAfterEvalTime = numberOfDescendantsFromUniqueEventAndSameStateAfterEvalTime;
        }

        public String toString() {
            return treeId + DELIMITER + evaluationTime + DELIMITER + ancestralTime + DELIMITER
                    + startNodeID + DELIMITER + endNodeID + DELIMITER
                    + startLocation + DELIMITER + endLocation + DELIMITER
                    + time + DELIMITER + independenceTime + DELIMITER + numberOfDescendants + DELIMITER + numberOfDescendantsOfSameState
                    + DELIMITER + totalNumberOfDescendantsFromUniqueEvent + DELIMITER + totalnumberOfDescendantsFromUniqueEventAndSameState
                    + DELIMITER + numberOfDescendantsFromUniqueEventAfterEvalTime + DELIMITER + numberOfDescendantsFromUniqueEventAndSameStateAfterEvalTime
                    ;
        }
    }


    private double adjust(double time) {
        if (mrsd < Double.MAX_VALUE) {
            time = mrsd - time;
        }
        return time;
    }


    private static Object[] readCJH(NodeRef node, Tree treeTime) {
        if (treeTime.getNodeAttribute(node, HISTORY) != null) {
            return (Object[]) treeTime.getNodeAttribute(node, HISTORY);
        } else {
            return null;
        }
    }

    protected PrintStream openOutputFile(String outputFileName) {
        PrintStream ps = super.openOutputFile(outputFileName);
        ps.println("treeId,evaluationTime,ancestralTime,evaluationNodeID,ancestralNodeID,stateAtEvaluationTime,ancestralState,persistenceTime,independenceTime,descendants,descendantsOfSameState,totalDescendantsFromUnique,totalDescendantsFromUniqueOfSameState,descendantsFromUniqueAfterEvalTime,descendantsFromUniqueOfSameStateAfterEvalTime");
        return ps;
    }

    private void closeOutputFile(PrintStream ps) {
        if (ps != null) {
            ps.close();
        }
    }

    private double mrsd;
//    private double ancestryTime;

    private PrintStream ps;

    public static void printTitle() {
        progressStream.println();
        centreLine("PersistVsIntroSummarizer " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("tool to summarize state persistence at a particular timepoint", 60);
        centreLine("by", 60);
        centreLine("Philippe Lemey and Marc Suchard", 60);
        progressStream.println();
        progressStream.println();
    }

    public static void printUsage(Arguments arguments) {

        arguments.printUsage("PersistVsIntroSummarizer", "<input-file-name> [<output-file-name>]");
        progressStream.println();
        progressStream.println("  Example: PersistVsIntroSummarizer -burnIn 1000 -evaluationTime 0.5 -ancestralTime 1.0 -nodeStateAnnotation states input.trees output.trees");
        progressStream.println();
    }

    private static final String EVALUATION_TIME = "evaluationTime";
    private static final String ANCESTRAL_TIME = "ancestralTime";
    private static final String NODE_STATE_ANNOTATION = "nodeStateAnnotation";

    //Main method
    public static void main(String[] args) throws IOException, Arguments.ArgumentException {

        int burnIn = -1;
        double[] evaluationTimes = new double[1];
        double[] ancestralTimes = new double[]{Double.MAX_VALUE};
        String nodeStateAnnotation = null;

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption(BURN_IN, "the number of states to be considered as 'burn-in' [default = 0]"),
//                        new Arguments.RealOption("mrsd", "The most recent sampling time to convert heights to times [default=MAX_VALUE]"),
                        new Arguments.RealArrayOption(EVALUATION_TIME, -1, "The time(s) at which the ancestral persistence of lineages is evaluated"),
//                        new Arguments.RealOption("independenceTime", "The time for which a lineage should not share a common ancestor with another lineage to be called a unique persistence/introduction  [default=MAX_VALUE]"),
                        new Arguments.RealArrayOption(ANCESTRAL_TIME, -1, "The time(s) in the past until which the the ancestral persistence of lineages is evaluated [default=MAX_VALUE]"),
                        new Arguments.StringOption(NODE_STATE_ANNOTATION, "String", "use node state annotations as poor proxy to MJs based on a annotation string for the discrete trait"),
                        new Arguments.Option("help", "option to print this message"),
                });


        handleHelp(arguments, args, TaxaMarkovJumpHistoryAnalyzer::printUsage);

//        if (arguments.hasOption("mrsd")) {
//            mrsd = arguments.getRealOption("mrsd");
//        }

        if (arguments.hasOption(EVALUATION_TIME)) {
            evaluationTimes = arguments.getRealArrayOption(EVALUATION_TIME);
        }

        if (arguments.hasOption(ANCESTRAL_TIME)) {
            ancestralTimes = arguments.getRealArrayOption(ANCESTRAL_TIME);
        }

        if (evaluationTimes.length != ancestralTimes.length) {
            throw new Arguments.ArgumentException("The number of " + EVALUATION_TIME + " arguments (" +
                    evaluationTimes.length + ") must equal the number of " + ANCESTRAL_TIME +
                    " arguments (" + ancestralTimes.length + ")");
        }

        if (arguments.hasOption(BURN_IN)) {
            burnIn = arguments.getIntegerOption(BURN_IN);
            System.err.println("Ignoring a burn-in of " + burnIn + " trees.");
        }

        if (arguments.hasOption(NODE_STATE_ANNOTATION)) {
            nodeStateAnnotation = arguments.getStringOption(NODE_STATE_ANNOTATION);
        } else {
            throw new Arguments.ArgumentException("Must include state annotation name via argument " +
                    NODE_STATE_ANNOTATION + " (e.g. -" + NODE_STATE_ANNOTATION + " yourStateName)");
        }

        String[] fileNames = getInputOutputFileNames(arguments, PersistenceSummarizer::printUsage);

        new PersistenceSummarizer(fileNames[0], fileNames[1], burnIn,
 //               mrsd,
                evaluationTimes,
                ancestralTimes,
                nodeStateAnnotation
        );
        System.exit(0);
    }
}