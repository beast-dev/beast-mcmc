package dr.app.tools;

import dr.app.beast.BeastVersion;
import dr.app.util.Arguments;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.util.DataTable;
import dr.util.Version;

import javax.json.Json;
import javax.json.JsonObject;
import javax.json.JsonObjectBuilder;
import java.io.FileReader;
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
    private static final String MERGE_STATES = "mergeStates";
    private static final String TIME_SLICES = "timeSlice";
    private static final String MOST_RECENT_SAMPLING_DATE = "mrsd";
    private PrintStream ps;
    private HashMap<String, String> mergedStates;
    private double timeSlice = Double.MAX_VALUE;
    public static final String SPECIAL_CHARACTERS_REGEX = ".*[\\s\\.;,\"\'].*";
    private HashMap<String, Integer> idMap = new HashMap<String, Integer>();

    public TransmissionChainSummarizer(
            String inputFileName,
            String outputFileName,
            int burnIn,
            String nodeStateAnnotation,
            String[] annotationStates,
            HashMap<String, String> mergedStates,
            double timeSlice
    ) throws IOException {
        SequentialTreeReader treeReader = new SequentialTreeReader(inputFileName, burnIn);

        this.ps = openOutputFile(outputFileName);
        this.mergedStates = mergedStates;
        this.timeSlice = timeSlice;
        processTrees(treeReader, burnIn, nodeStateAnnotation, annotationStates);
        closeOutputFile(ps);
    }

    protected PrintStream openOutputFile(String outputFileName) {
        PrintStream ps = super.openOutputFile(outputFileName);
        return ps;
    }

    private void closeOutputFile(PrintStream ps) {
        if (ps != null) {
            ps.close();
        }
    }

    private HashMap<String, Integer> createAndPrintIdMap(Tree tree){
        int k = 1;
        int taxonCount = tree.getTaxonCount();
        List<String> names = new ArrayList<String>();
        for (int i = 0; i < tree.getTaxonCount(); i++) {
            names.add(tree.getTaxonId(i));
        }
        HashMap<String, Integer> idMap = new HashMap<String, Integer>();
        for (String name : names) {
            idMap.put(name, k);
            if (name.matches(SPECIAL_CHARACTERS_REGEX)) {
                name = "'" + name + "'";
            }
            if (k < names.size()) {
                this.ps.println("#"+k + "\t" + name);
            }
            k += 1;
        }
        return idMap;
    }

    private void processTrees(SequentialTreeReader treeReader, int burnIn, String nodeStateAnnotation, String[] annotationStates) throws IOException {
        if (burnIn < 0) {
            burnIn = 0;
        }

        int index = burnIn;
        Tree tree;

        while (treeReader.getTree(index) != null) {
            tree = treeReader.getTree(index);
            if(index == burnIn){
                this.idMap = createAndPrintIdMap(tree);
                this.ps.println(Row.header);
            }
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
        int numberOfBranches;
        int numberOfDescendants;
        int numberOfDescendantsOfSameState;
        int numberOfPersistentDescendantsOfSameState;
        String persistentDescendantIds;
        double lengthOfTransmissionChain;
        double rootHeight;
        boolean isExternal;
        String jumpsToDifferentState;

        private static final String DELIMITER = "\t";

        private static final String header = "treeId\t" +
                "startNodeID\t"+
                "endNodeID\t"+
                "sourceLocation\t"+
                "currentLocation\t"+
                "introductionTime\t"+
                "numberOfBranches\t"+
                "numberOfDescendants\t"+
                "numberOfDescendantsOfSameState\t"+
                "numberOfPersistentDescendantsOfSameState\t"+
                "persistentDescendantsOfSameState\t"+
                "lengthOfTransmissionChain\t"+
                "rootHeight\t"+
                "isExternal\t"+
                "jumpsToDifferentState";

        private Row(
                String treeId,
                int currentNodeID,
                int ancestralNodeID,
                String sourceLocation,
                String currentLocation,
                double introductionTime,
                int numberOfBranches,
                int numberOfDescendants,
                int numberOfDescendantsOfSameState,
                int numberOfPersistentDescendantsOfSameState,
                String persistentDescendantIds,
                double lengthOfTransmissionChain,
                double rootHeight,
                boolean isExternal,
                String jumpsToDifferentState
        ) {
            this.treeId = treeId;
            this.currentNodeID = currentNodeID;
            this.ancestralNodeID = ancestralNodeID;
            this.sourceLocation = sourceLocation;
            this.currentLocation = currentLocation;
            this.introductionTime = introductionTime;
            this.numberOfBranches = numberOfBranches;
            this.numberOfDescendants = numberOfDescendants;
            this.numberOfDescendantsOfSameState = numberOfDescendantsOfSameState;
            this.numberOfPersistentDescendantsOfSameState = numberOfPersistentDescendantsOfSameState;
            this.persistentDescendantIds = persistentDescendantIds;
            this.lengthOfTransmissionChain = lengthOfTransmissionChain;
            this.rootHeight = rootHeight;
            this.isExternal = isExternal;
            this.jumpsToDifferentState = jumpsToDifferentState;
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
                    String.valueOf(numberOfBranches),
                    String.valueOf(numberOfDescendants),
                    String.valueOf(numberOfDescendantsOfSameState),
                    String.valueOf(numberOfPersistentDescendantsOfSameState),
                    persistentDescendantIds,
                    String.valueOf(lengthOfTransmissionChain),
                    String.valueOf(rootHeight),
                    String.valueOf(isExternal),
                    jumpsToDifferentState
            );
        }
    }

    private String idstoString(Tree tree, Set<NodeRef> persistentDescendants){
        Iterator iter = persistentDescendants.iterator();
        List<String> ids = new ArrayList<String>();
        int ctr = 0;
        NodeRef currentNode;
        while (iter.hasNext()) {
            currentNode = (NodeRef)iter.next();
            ids.add(String.valueOf(this.idMap.get(tree.getNodeTaxon(currentNode).getId())));
        }
        return String.join(",", ids);
    }

    private String convertToJson(HashMap<String, Integer> map){
        JsonObjectBuilder builder = Json.createObjectBuilder();
        map.forEach(builder::add);
        JsonObject obj = builder.build();
        return obj.toString();
    }

    // Get node state after merging states
    private String getMergedState(String nodeState){
        String key = nodeState.toLowerCase();
        if(this.mergedStates.containsKey(key)){
            return this.mergedStates.get(key);
        }
        return nodeState;
    }

    private String getMergedState(Tree tree, NodeRef node, String nodeStateAnnotation){
        String nodeState = (String) tree.getNodeAttribute(node, nodeStateAnnotation);
        String key = nodeState.toLowerCase();
        if(this.mergedStates.containsKey(key)){
            return this.mergedStates.get(key);
        }
        return nodeState;
    }

    private void processOneTree(Tree tree, String nodeStateAnnotation, List<String> annotationStates){
        String treeId = tree.getId();
        if (treeId.startsWith("STATE_")) {
            treeId = treeId.replaceFirst("STATE_", "");
        }

        for (int i = 0; i < tree.getNodeCount(); ++i) {
            NodeRef node = tree.getNode(i);
            double nodeHeight = tree.getNodeHeight(node);
//            if(tree.isRoot(node))
//                continue;

            String nodeState = getMergedState(tree, node, nodeStateAnnotation);

            if(!annotationStates.contains(nodeState) || nodeHeight < this.timeSlice)
                continue;

            NodeRef parentNode = tree.isRoot(node) ? null : tree.getParent(node);
            double parentNodeHeight = tree.isRoot(node) ? tree.getNodeHeight(node) : tree.getNodeHeight(parentNode);

            Set<NodeRef> nodeDescendants = getDescendants(tree, node);
            if (nodeState == null) {
                throw new RuntimeException("Could not locate node state annotation '" + nodeStateAnnotation +
                        "' for node " + node.getNumber());
            }

            String currentState = nodeState;
            String ancestralState = tree.isRoot(node) ? null : getMergedState(tree, parentNode, nodeStateAnnotation);

            // If state of parentNode is the same as the current node skip
            if(currentState.equalsIgnoreCase(ancestralState))
                continue;

            double introductionTime = 0;
            Object[] jumps = readCJH(node, tree);

            //When no jumps are found, record height of parent node
            if (jumps == null) {
                // Unclear if this is even required here
                introductionTime = parentNodeHeight;
            } else {
                // If jumps are found, record the most recent jump
                Object[] recentJump = getMostRecentJump(jumps);
                ancestralState = (String) recentJump[1];
                currentState = (String) recentJump[2];
                introductionTime = (Double) recentJump[0];
            }

            // External nodes in a transmission chain that persist within a given state
            Set<NodeRef> persistentDescendants = new HashSet<NodeRef>();
            // Count of each immediate jump from a transmission chain into a different state
            HashMap<String, Integer> map = new HashMap<String, Integer>();
            // Length of persistence of transmission chain in same state
            double heightOfTransmissionChain = nodeHeight;

            heightOfTransmissionChain = traversePersistentChain(
                    tree,
                    node,
                    nodeStateAnnotation,
                    nodeState,
                    persistentDescendants,
                    map,
                    heightOfTransmissionChain
            );

            double rootHeight = tree.getNodeHeight(tree.getRoot());
            Row row = new Row(
                    treeId,
                    parentNode == null ? -1 : parentNode.getNumber(),
                    node.getNumber(),
                    ancestralState,
                    currentState,
                    introductionTime,
                    getNumberOfBranches(tree, node, currentState, nodeStateAnnotation),
                    nodeDescendants.size(),
                    getSameStateDescendants(nodeDescendants, tree, currentState, nodeStateAnnotation),
                    persistentDescendants.size(),
                    idstoString(tree, persistentDescendants),
                    (introductionTime - heightOfTransmissionChain)/rootHeight,
                    rootHeight,
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

    private Set<NodeRef> getDescendants(Tree tree, NodeRef node){
        Set<NodeRef> descendants = new HashSet<NodeRef>();
        Set<NodeRef> nodeDescendants = TreeUtils.getExternalNodes(tree, node);
        Iterator iter = nodeDescendants.iterator();
        while (iter.hasNext()) {
            NodeRef currentNode = (NodeRef)iter.next();
            if(tree.getNodeHeight(currentNode) < this.timeSlice)
                continue;
            descendants.add(currentNode);
        }
        return descendants;
    }

    // Get all descendants in the same state regardless of persistence
    private int getSameStateDescendants(Set leafs, Tree tree, String state, String nodeStateAnnotation){
        int descendents = 0;
        Iterator iter = leafs.iterator();
        while (iter.hasNext()) {
            NodeRef currentNode = (NodeRef)iter.next();
            String nodeState = getMergedState(tree, currentNode, nodeStateAnnotation);
            if (nodeState.equalsIgnoreCase(state)){
                descendents++;
            }
        }
        return descendents;
    }

    // Get number of persistent branches i.e., number of branches of the transmission chain in the same state
    private int getNumberOfBranches(Tree tree, NodeRef node, String state, String nodeStateAnnotation){
        int childCount = tree.getChildCount(node);
        int numberOfBranches = childCount;
        for (int i = 0; i < childCount; i++) {
            NodeRef childNode = tree.getChild(node, i);
            if(tree.getNodeHeight(childNode) < this.timeSlice)
                continue;
            String nodeState = getMergedState(tree, childNode, nodeStateAnnotation);
            if (!nodeState.equalsIgnoreCase(state)){
                continue;
            }
            numberOfBranches += getNumberOfBranches(tree, childNode, state, nodeStateAnnotation);
        }
        return numberOfBranches;
    }


    private double traversePersistentChain(
            Tree tree,
            NodeRef node,
            String nodeStateAnnotation,
            String annotationState,
            Set<NodeRef> persistentDescendants,
            HashMap<String, Integer> jumpsToDifferentState,
            double minHeight) {
        Integer count = 0;
        double currHeight = 0;
        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef childNode = tree.getChild(node, i);
            if(tree.getNodeHeight(childNode) < this.timeSlice) {
                if(this.timeSlice < minHeight)
                    minHeight = this.timeSlice;
                continue;
            }
            String childNodeState = getMergedState(tree, childNode, nodeStateAnnotation);
            if(tree.isExternal(childNode)) {
                if (childNodeState.equalsIgnoreCase(annotationState)){
                    // If node is external and is in same state as chain, add to persistent descendants
                    persistentDescendants.add(childNode);
                    // If node is external and is in same state as chain, update minimum height as required
                    currHeight = tree.getNodeHeight(childNode);
                    if(currHeight < minHeight){
                        minHeight = currHeight;
                    }
                } else {
                    // If a leaf has different state, increment count of jumps out of transmission chain to a different state
                    count = jumpsToDifferentState.containsKey(childNodeState) ? jumpsToDifferentState.get(childNodeState) : 0;
                    count += 1;
                    jumpsToDifferentState.put(childNodeState, count);
                    // If node is external and is in different state as chain, check MJs to update minimum height as required
                    Object[] jumps = readCJH(childNode, tree);
                    Object[] earliestJump;
                    if(jumps != null){
                        earliestJump = getEarliestJump(jumps);
                        if(!annotationState.equalsIgnoreCase((String) earliestJump[1])){
                            throw new RuntimeException("Persistent chain state, "+annotationState+" and source of earliest jump"+(String) earliestJump[1]+" do not match!");
                        }
                        currHeight = (Double) earliestJump[0];
                    }
                    if(currHeight < minHeight){
                        minHeight = currHeight;
                    }
                }
            } else {
                if(!childNodeState.equalsIgnoreCase(annotationState)){
                    // If state of internal node doesn't match state of transmission chain, increment count of jumps out of transmission chain to a different state
                    count = jumpsToDifferentState.containsKey(childNodeState) ? jumpsToDifferentState.get(childNodeState) : 0;
                    count += 1;
                    jumpsToDifferentState.put(childNodeState, count);

                    // If node is internal and is in different state as chain, check MJs to update minimum height as required
                    Object[] jumps = readCJH(childNode, tree);
                    Object[] earliestJump;
                    if(jumps != null){
                        earliestJump = getEarliestJump(jumps);
                        if(!annotationState.equalsIgnoreCase((String) earliestJump[1])){
                            throw new RuntimeException("Persistent chain state, "+annotationState+" and source of earliest jump, "+(String) earliestJump[1]+" do not match!");
                        }
                        currHeight = (Double) earliestJump[0];
                    } else {
                        // If Mjs are null use nodeHeight
                        currHeight = tree.getNodeHeight(childNode);
                    }
                    if(currHeight < minHeight){
                        minHeight = currHeight;
                    }
                    continue;
                }
            }
            minHeight = traversePersistentChain(tree, childNode, nodeStateAnnotation, annotationState, persistentDescendants, jumpsToDifferentState, minHeight);
        }
        return minHeight;
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

    private Object[] readCJH(NodeRef node, Tree treeTime) {
        if (treeTime.getNodeAttribute(node, HISTORY) != null) {
            Object[] jumps = (Object[]) treeTime.getNodeAttribute(node, HISTORY);
            for (int j = jumps.length - 1; j >= 0; j--) {
                Object[] currentJump = (Object[]) jumps[j];
                currentJump[1] = getMergedState((String) currentJump[1]);
                currentJump[2] = getMergedState((String) currentJump[2]);
            }
            return jumps;
        }
        return null;
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
        String mergedStatesFileName = new String();
        HashMap<String, String> mergedStates = new HashMap<String, String>();
        double timeSlice = 0;
        double mrsd = 0;

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption(BURN_IN, "the number of states to be considered as 'burn-in' [default = 0]"),
                        new Arguments.StringOption(NODE_STATE_ANNOTATION, "String", "use node state annotations as poor proxy to MJs based on a annotation string for the discrete trait"),
                        new Arguments.StringOption(ANNOTATION_STATES, "String", "States to consider"),
                        new Arguments.StringOption(MERGE_STATES, "String", "TSV file with states to be merged into a new state"),
                        new Arguments.StringOption(TIME_SLICES, "String", "Get number of descendants at given node height or time (if mrsd is supplied)"),
                        new Arguments.StringOption(MOST_RECENT_SAMPLING_DATE, "String", "Most recent sampling date"),
                        new Arguments.Option("help", "option to print this message"),
                });


        handleHelp(arguments, args, TransmissionChainSummarizer::printUsage);


        if (arguments.hasOption(ANNOTATION_STATES)) {
            annotationStates = arguments.getStringOption(ANNOTATION_STATES).split(",");
        }

        if (arguments.hasOption(MERGE_STATES)) {
            mergedStatesFileName = arguments.getStringOption(MERGE_STATES);
            DataTable<String[]> dataTable = DataTable.Text.parse(new FileReader(mergedStatesFileName));
            if(dataTable.getColumnCount() != 1){
                throw new Arguments.ArgumentException(MERGE_STATES + " file should have two columns: State and New_State");
            }
            String[] rowLabels = dataTable.getRowLabels();
            for (int i = 0; i <dataTable.getRowCount(); i++) {
                String[] states = dataTable.getRow(i);
                mergedStates.put(rowLabels[i].toLowerCase(), states[0]);
            }
        }

        if(arguments.hasOption(TIME_SLICES)){
            timeSlice =  Double.parseDouble(arguments.getStringOption(TIME_SLICES));
        }

        if(arguments.hasOption(MOST_RECENT_SAMPLING_DATE)){
            mrsd = Double.parseDouble(arguments.getStringOption(MOST_RECENT_SAMPLING_DATE));
            timeSlice = mrsd - timeSlice;
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
                annotationStates,
                mergedStates,
                timeSlice
        );
        System.exit(0);
    }
}
