package dr.app.tools;

import dr.app.beast.BeastVersion;
import dr.app.util.Arguments;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.util.Version;

import javax.json.Json;
import javax.json.JsonObject;
import javax.json.JsonObjectBuilder;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;
import java.util.function.Consumer;

public class TransmissionChainSummarizer extends BaseTreeTool {
    private final static Version version = new BeastVersion();

    private static final String HISTORY = "history";
    private static final String BURN_IN = "burnIn";
    private static final String NODE_STATE_ANNOTATION = "nodeStateAnnotation";
    private static final String ANNOTATION_STATES = "annotationStates";
    private PrintStream ps;

    public TransmissionChainSummarizer(
            String inputFileName,
            String outputFileName,
            int burnIn,
            String nodeStateAnnotation,
            String[] annotationStates
    ) throws IOException {
        SequentialTreeReader treeReader = new SequentialTreeReader(inputFileName, burnIn);

        this.ps = openOutputFile(outputFileName);
        processTrees(treeReader, burnIn, nodeStateAnnotation, annotationStates);
        closeOutputFile(ps);
    }

    protected PrintStream openOutputFile(String outputFileName) {
        PrintStream ps = super.openOutputFile(outputFileName);
        ps.println(Row.header);
        return ps;
    }

    private void closeOutputFile(PrintStream ps) {
        if (ps != null) {
            ps.close();
        }
    }

    private void processTrees(SequentialTreeReader treeReader, int burnIn, String nodeStateAnnotation, String[] annotationStates) throws IOException {
        if (burnIn < 0) {
            burnIn = 0;
        }

        int index = burnIn;
        Tree tree;

        while (treeReader.getTree(index) != null) {
            tree = treeReader.getTree(index);
            processOneTree(tree, nodeStateAnnotation, Arrays.asList(annotationStates));
            index++;
        }
    }

    private class Row {
        String treeId;
        int currentNodeID;
        int ancestralNodeID;
        String sourceLocation;
        String currentLocation;
        double introductionTime;
        int numberOfDescendants;
        int numberOfDescendantsOfSameState;
        int numberOfPersistentDescendantsOfSameState;
        double lengthOfTransmissionChain;
        double rootHeight;
        boolean isExternal;
        String persistentDescendantStateCounts;

        private static final String DELIMITER = "\t";

        private static final String header = "treeId\t" +
                "startNodeID\t"+
                "endNodeID\t"+
                "sourceLocation\t"+
                "currentLocation\t"+
                "introductionTime\t"+
                "numberOfDescendants\t"+
                "numberOfDescendantsOfSameState\t"+
                "numberOfPersistentDescendantsOfSameState\t"+
                "lengthOfTransmissionChain\t"+
                "rootHeight\t"+
                "isExternal\t"+
                "persistentDescendantStateCounts";

        private Row(
                String treeId,
                int currentNodeID,
                int ancestralNodeID,
                String sourceLocation,
                String currentLocation,
                double introductionTime,
                int numberOfDescendants,
                int numberOfDescendantsOfSameState,
                int numberOfPersistentDescendantsOfSameState,
                double lengthOfTransmissionChain,
                double rootHeight,
                boolean isExternal,
                String persistentDescendantStateCounts
        ) {
            this.treeId = treeId;
            this.currentNodeID = currentNodeID;
            this.ancestralNodeID = ancestralNodeID;
            this.sourceLocation = sourceLocation;
            this.currentLocation = currentLocation;
            this.introductionTime = introductionTime;
            this.numberOfDescendants = numberOfDescendants;
            this.numberOfDescendantsOfSameState = numberOfDescendantsOfSameState;
            this.numberOfPersistentDescendantsOfSameState = numberOfPersistentDescendantsOfSameState;
            this.lengthOfTransmissionChain = lengthOfTransmissionChain;
            this.rootHeight = rootHeight;
            this.isExternal = isExternal;
            this.persistentDescendantStateCounts = persistentDescendantStateCounts;
        }

        public String toString() {
            return String.join(
                    DELIMITER,
                    treeId,
                    String.valueOf(currentNodeID),
                    String.valueOf(ancestralNodeID),
                    sourceLocation,
                    currentLocation,
                    String.valueOf(introductionTime),
                    String.valueOf(numberOfDescendants),
                    String.valueOf(numberOfDescendantsOfSameState),
                    String.valueOf(numberOfPersistentDescendantsOfSameState),
                    String.valueOf(lengthOfTransmissionChain),
                    String.valueOf(rootHeight),
                    String.valueOf(isExternal),
                    persistentDescendantStateCounts
            );
        }
    }

    private String convertToJson(HashMap<String, Integer> map){
        JsonObjectBuilder builder = Json.createObjectBuilder();
        map.forEach(builder::add);
        JsonObject obj = builder.build();
        return obj.toString();
    }

    private void processOneTree(Tree tree, String nodeStateAnnotation, List<String> annotationStates){
        String treeId = tree.getId();
        if (treeId.startsWith("STATE_")) {
            treeId = treeId.replaceFirst("STATE_", "");
        }

        for (int i = 0; i < tree.getNodeCount(); ++i) {
            NodeRef node = tree.getNode(i);
            if(tree.isRoot(node))
                continue;

            String nodeState = (String) tree.getNodeAttribute(node, nodeStateAnnotation);

            if(!annotationStates.contains(nodeState))
                continue;

            double nodeHeight = tree.getNodeHeight(node);
            NodeRef parentNode = tree.getParent(node);
            double parentNodeHeight = tree.getNodeHeight(parentNode);

            Set nodeDescendants = TreeUtils.getExternalNodes(tree, node);
            Set totalEventDescendants = null;

            if (nodeState == null) {
                throw new RuntimeException("Could not locate node state annotation '" + nodeStateAnnotation +
                        "' for node " + node.getNumber());
            }

            String currentState = nodeState;
            String ancestralState =  (String) tree.getNodeAttribute(parentNode, nodeStateAnnotation);;

            // If state of parentNode is the same as the current node skip
            if(currentState.equalsIgnoreCase(ancestralState))
                continue;

            double introductionTime = 0;
            Object[] jumps = readCJH(node, tree);

            //When no jumps are found, record height of parent node
            if (jumps == null) {
                // Unclear if this is even required here
                totalEventDescendants = TreeUtils.getExternalNodes(tree, node);
                introductionTime = parentNodeHeight;
            } else {
                // If jumps are found, record the most recent jump
                Object[] recentJump = getMostRecentJump(jumps);
                ancestralState = (String) recentJump[1];
                currentState = (String) recentJump[2];
                introductionTime = (Double) recentJump[0];
                totalEventDescendants = TreeUtils.getExternalNodes(tree, node);
            }

            Set<NodeRef> persistentDescendants = new HashSet<NodeRef>();
            getPersistentDescendants(tree, node, nodeStateAnnotation, nodeState, persistentDescendants);

            HashMap<String, Integer> map = new HashMap<String, Integer>();
            getPersistentDescendantStateCounts(tree, node, nodeStateAnnotation,  nodeState, map);


            Row row = new Row(
                    treeId,
                    parentNode.getNumber(),
                    node.getNumber(),
                    ancestralState,
                    currentState,
                    introductionTime,
                    totalEventDescendants.size(),
                    getSameStateDescendants(nodeDescendants,tree,currentState,nodeStateAnnotation),
                    persistentDescendants.size(),
                    getLengthOfTransmissionChain(tree, node, nodeState, persistentDescendants),
                    tree.getNodeHeight(tree.getRoot()),
                    tree.isExternal(node),
                    convertToJson(map)
                );
            ps.println(row);

        }
    }

    private Object[] getEarliestJump(Object[] jumps){
        double earliestJumpTime = Double.MIN_VALUE;
        Object[] earliestJump = new Object[3];
        //get earliest jump
        for (int j = jumps.length - 1; j >= 0; j--) {
            Object[] currentJump = (Object[]) jumps[j];
            double jumpTime = (Double) currentJump[0];
            if (jumpTime > earliestJumpTime) {
                earliestJumpTime = jumpTime;
                earliestJump = currentJump;
            }
        }
        return earliestJump;
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

    private int getSameStateDescendants(Set leafs, Tree tree, String state, String nodeStateAnnotation){
        int descendents = 0;
        Iterator iter = leafs.iterator();
        while (iter.hasNext()) {
            NodeRef currentNode = (NodeRef)iter.next();
            String nodeState = (String) tree.getNodeAttribute(currentNode, nodeStateAnnotation);
            if (nodeState.equalsIgnoreCase(state)){
                descendents++;
            }
        }
        return descendents;
    }

    // Get all nodes in a transmission chain that persists within a given state
    private void getPersistentDescendants(Tree tree, NodeRef node, String nodeStateAnnotation, String annotationState, Set<NodeRef> set){
        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef childNode = tree.getChild(node, i);
            String childNodeState = (String) tree.getNodeAttribute(childNode, nodeStateAnnotation);
            if(!childNodeState.equalsIgnoreCase(annotationState))
                continue;
            if(tree.isExternal(childNode))
                set.add(childNode);
            getPersistentDescendants(tree, childNode, nodeStateAnnotation, annotationState, set);
        }
    }

    // Get count of states of each immediate jump from a transmission chain in a different state
    private void getPersistentDescendantStateCounts(Tree tree, NodeRef node, String nodeStateAnnotation, String annotationState, HashMap<String, Integer> map){
        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef childNode = tree.getChild(node, i);
            String childNodeState = (String) tree.getNodeAttribute(childNode, nodeStateAnnotation);
            if(tree.isExternal(childNode)) {
                Integer count = map.containsKey(childNodeState) ? map.get(childNodeState) : 0;
                count += 1;
                map.put(childNodeState, count);
            } else if(!childNodeState.equalsIgnoreCase(annotationState)) {
                Integer count = map.containsKey(childNodeState) ? map.get(childNodeState) : 0;
                count += 1;
                map.put(childNodeState, count);
                continue;
            }
            // Only traverse path if childNode is in same state as annotationState
            getPersistentDescendantStateCounts(tree, childNode, nodeStateAnnotation, annotationState, map);
        }
    }

    // Length of the longest persisting path in clade as a proportion of rootHeight
    private double getLengthOfTransmissionChain(Tree tree, NodeRef node, String nodeState, Set<NodeRef> persistentDescendants){
        if(persistentDescendants.size() == 0)
            return 0;
        double rootHeight = tree.getNodeHeight(tree.getRoot());
        double minHeight = Double.MAX_VALUE;
        double nodeHeight = tree.getNodeHeight(node);
        Iterator<NodeRef> itr = persistentDescendants.iterator();
        NodeRef currNode;
        double currHeight;
        while(itr.hasNext()){
            currNode = itr.next();
            currHeight = tree.getNodeHeight(currNode);
            if(tree.isExternal(currNode)){
                currHeight = tree.getNodeHeight(currNode);
                if(currHeight < minHeight){
                    minHeight = currHeight;
                }
            } else {
                // Iterate over MJs of child node to get earliest jump time from source, nodeState
                for(int i =0;i<tree.getChildCount(currNode);i++){
                    Object[] jumps = readCJH(tree.getChild(currNode, i), tree);
                    Object[] earliestJump;
                    if(jumps != null){
                        earliestJump = getEarliestJump(jumps);
                        if(!nodeState.equalsIgnoreCase((String) earliestJump[1])){
                            throw new RuntimeException("Node state, "+nodeState+" and source of earliest jump"+(String) earliestJump[1]+" do not match!");
                        }
                        currHeight = (Double) earliestJump[0];
                    }
                    if(currHeight < minHeight){
                        minHeight = currHeight;
                    }
                }
            }
        }
        return (nodeHeight - minHeight)/rootHeight;
    }

    private static Object[] readCJH(NodeRef node, Tree treeTime) {
        if (treeTime.getNodeAttribute(node, HISTORY) != null) {
            return (Object[]) treeTime.getNodeAttribute(node, HISTORY);
        } else {
            return null;
        }
    }

    public static void printTitle() {
        progressStream.println();
        centreLine("TransmissionChainSummarizer " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("tool to summarize transmission chains", 60);
        centreLine("by", 60);
        centreLine("Karthik G and Marc Suchard", 60);
        progressStream.println();
        progressStream.println();
    }

    public static void printUsage(Arguments arguments) {

        arguments.printUsage("TransmissionChainSummarizer", "<input-file-name> [<output-file-name>]");
        progressStream.println();
        progressStream.println("  Example: TransmissionChainSummarizer -burnIn 1000 -nodeStateAnnotation states -annotateStates location1,location2 input.trees output.tsv");
        progressStream.println();
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

    //Main method
    public static void main(String[] args) throws IOException, Arguments.ArgumentException {

        int burnIn = -1;
        String nodeStateAnnotation = null;
        String[] annotationStates = new String[1];

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption(BURN_IN, "the number of states to be considered as 'burn-in' [default = 0]"),
                        new Arguments.StringOption(NODE_STATE_ANNOTATION, "String", "use node state annotations as poor proxy to MJs based on a annotation string for the discrete trait"),
                        new Arguments.StringOption(ANNOTATION_STATES, "String", "States to consider"),
                        new Arguments.Option("help", "option to print this message"),
                });


        handleHelp(arguments, args, TransmissionChainSummarizer::printUsage);


        if (arguments.hasOption(ANNOTATION_STATES)) {
            annotationStates = arguments.getStringOption(ANNOTATION_STATES).split(",");
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

        new TransmissionChainSummarizer(
                fileNames[0],
                fileNames[1],
                burnIn,
                nodeStateAnnotation,
                annotationStates
        );
        System.exit(0);
    }
}
