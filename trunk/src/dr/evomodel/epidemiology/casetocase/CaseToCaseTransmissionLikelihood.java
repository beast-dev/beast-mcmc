package dr.evomodel.epidemiology.casetocase;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.MathUtils;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.*;

import java.io.*;
import java.util.*;

/**
 * A likelihood function for transmission between identified epidemiological cases
 *
 * Currently works only for fixed trees and estimates the network and epidemiological parameters. Timescale
 * must be in days. Python scripts to write XML for it and analyse the posterior set of networks exist; contact MH.
 *
 * @author Matthew Hall
 * @author Matthew Hall
 * @author Andrew Rambaut
 * @version $Id: $
 */
public class CaseToCaseTransmissionLikelihood extends AbstractModelLikelihood implements Loggable, Citable {

    /* The phylogenetic tree. */

    private TreeModel virusTree;

    /* Mapping of cases to branches on the tree; old version is stored before operators are applied */

    private AbstractCase[] branchMap;
    private AbstractCase[] storedBranchMap;

    /* The log likelihood of the subtree from the parent branch of the referred-to node downwards; old version is stored
    before operators are applied.*/

    private double[] nodeLogLikelihoods;
    private double[] storedNodeLogLikelihoods;

    /* Whether operations have required the recalculation of the log likelihoods of the subtree from the parent branch
    of the referred-to node downwards.
     */

    private boolean[] nodeRecalculationNeeded;
    private boolean[] storedRecalculationArray;

    /* Matches cases to external nodes */

    private HashMap<AbstractCase, NodeRef> tipMap;
    private double estimatedLastSampleTime;
    boolean verbose;

    /**
     * The set of cases
     */
    private AbstractOutbreak cases;
    private boolean likelihoodKnown = false;
    private double logLikelihood;

    // for extended version

    private boolean extended;
    private boolean[] switchLocks;
    private boolean[] creepLocks;

    // occasionally needed for debugging purposes

/*    private HashMap<AbstractCase, Double> treeInfectionDates;
    private HashMap<AbstractCase, Integer> treeInfectionDays;
    private HashMap<AbstractCase, Double> treeExamDates;
    private HashMap<AbstractCase, Integer> treeExamDays;
    private HashMap<AbstractCase, Double> dataInfectiousnessDates;
    private HashMap<AbstractCase, Double> dataExamDates;*/

    // PUBLIC STUFF

    // Name

    public static final String CASE_TO_CASE_TRANSMISSION_LIKELIHOOD = "caseToCaseTransmissionLikelihood";

    // Basic constructor.

    public CaseToCaseTransmissionLikelihood(TreeModel virusTree, AbstractOutbreak caseData,
                                            String startingNetworkFileName, boolean extended)
            throws TaxonList.MissingTaxonException {
        this(CASE_TO_CASE_TRANSMISSION_LIKELIHOOD, virusTree, caseData, startingNetworkFileName, extended);
    }

    // Legacy constructor

    public CaseToCaseTransmissionLikelihood(String name, TreeModel virusTree, AbstractOutbreak caseData, String
            startingNetworkFileName){
        this(name, virusTree, caseData, startingNetworkFileName, false);
    }


    // Constructor for an instance with a non-default name

    public CaseToCaseTransmissionLikelihood(String name, TreeModel virusTree, AbstractOutbreak caseData, String
            startingNetworkFileName, boolean extended) {

        super(name);
        this.virusTree = virusTree;
        addModel(virusTree);
        Date lastSampleDate = getLatestTaxonDate(virusTree);

        /* We assume samples were taken at the end of the date of the last tip */

        estimatedLastSampleTime = lastSampleDate.getTimeValue();
        cases = caseData;
        this.extended = extended;
        addModel(cases);
        verbose = false;

        tipMap = new HashMap<AbstractCase, NodeRef>();

        nodeRecalculationNeeded = new boolean[virusTree.getNodeCount()];
        Arrays.fill(nodeRecalculationNeeded, true);
        nodeLogLikelihoods = new double[virusTree.getNodeCount()];

        if(startingNetworkFileName==null){
            branchMap = paintRandomNetwork();
        } else {
            branchMap = paintSpecificNetwork(startingNetworkFileName);
        }

        for(int i=0; i<virusTree.getExternalNodeCount(); i++){
            TreeModel.Node currentExternalNode = (TreeModel.Node)virusTree.getExternalNode(i);
            Taxon currentTaxon = currentExternalNode.taxon;
            for(AbstractCase thisCase : cases.getCases()){
                for(Taxon caseTaxon: thisCase.getAssociatedTaxa()){
                    if(caseTaxon.equals(currentTaxon)){
                        tipMap.put(thisCase, currentExternalNode);
                    }
                }
            }
        }
        if(extended){
            switchLocks = new boolean[virusTree.getInternalNodeCount()];
            creepLocks = new boolean[virusTree.getInternalNodeCount()];
            for(int i=0; i<virusTree.getInternalNodeCount(); i++){
                creepLocks[virusTree.getInternalNode(i).getNumber()] = isCreepLocked(virusTree.getInternalNode(i));
            }
        }

    }

    /* Get the date of the last tip */

    private static Date getLatestTaxonDate(TreeModel tree){
        Date latestDate = new Date(new java.util.Date(-1000000), Units.Type.DAYS);
        for(int i=0;i<tree.getTaxonCount();i++){
            if(tree.getTaxon(i).getDate().after(latestDate)){
                latestDate = tree.getTaxon(i).getDate();
            }
        }
        return latestDate;
    }

    public boolean tipLinked(NodeRef node){
        NodeRef tip = tipMap.get(branchMap[node.getNumber()]);
        if(tip==node){
            return true;
        }
        NodeRef parent = tip;
        while(parent!=virusTree.getRoot()){
            parent = virusTree.getParent(parent);
            if(parent==node){
                return true;
            }
        }
        return false;
    }

    public boolean[] getSwitchLocks(){
        return switchLocks;
    }

    public boolean[] getCreepLocks(){
        return creepLocks;
    }

    public boolean isSwitchLocked(NodeRef node){
        if(virusTree.isExternal(node)){
            throw new RuntimeException("Checking the lock on an external node");
        }
        return !tipLinked(node) || countChildrenWithSamePainting(node) == 2;
    }

    public boolean isCreepLocked(NodeRef node){
        if(virusTree.isExternal(node)){
            throw new RuntimeException("Checking the lock on an external node");
        }
        if(tipLinked(node) && samePaintingUpTree(node, false).contains(virusTree.getRoot().getNumber())){
            return true;
        } else if(tipLinked(node)){
            return doubleChildAncestor(node);
        } else {
            return doubleChildDescendant(node);
        }
    }

    // An ancestor of this node, or this node has the same painting and two children with that same painting

    public boolean doubleChildAncestor(NodeRef node){
        AbstractCase nodePainting = branchMap[node.getNumber()];
        NodeRef ancestor = node;
        while(ancestor!=null && branchMap[ancestor.getNumber()]==nodePainting){
            if(countChildrenWithSamePainting(ancestor)==2){
                return true;
            }
            ancestor=virusTree.getParent(ancestor);
        }
        return false;
    }

    //  A descendant of this node, or this node has the same painting and two children with that same painting

    public boolean doubleChildDescendant(NodeRef node){
        AbstractCase nodePainting = branchMap[node.getNumber()];
        return countChildrenWithSamePainting(node) == 2
                || doubleChildDescendant(virusTree.getChild(node, 0))
                || doubleChildDescendant(virusTree.getChild(node, 1));
    }

    public boolean isExtended(){
        return extended;
    }

    public void changeSwitchLock(int index, boolean value){
        switchLocks[index]=value;
    }

    public void recalculateLocks(){
        for(int i=0; i<virusTree.getInternalNodeCount(); i++){
            switchLocks[virusTree.getInternalNode(i).getNumber()]=isSwitchLocked(virusTree.getInternalNode(i));
            creepLocks[virusTree.getInternalNode(i).getNumber()]=isCreepLocked(virusTree.getInternalNode(i));
        }
    }

    //Counts the children of the current node which have the same painting as itself under the current map.
    //This will always be 1 if extended==false.

    public int countChildrenWithSamePainting(NodeRef node){
        if(virusTree.isExternal(node)){
            return -1;
        } else {
            int count = 0;
            AbstractCase parentCase = branchMap[node.getNumber()];
            for(int i=0; i<virusTree.getChildCount(node); i++){
                if(branchMap[virusTree.getChild(node,i).getNumber()]==parentCase){
                    count++;
                }
            }
            return count;
        }
    }


    private static NodeRef sibling(TreeModel tree, NodeRef node){
        if(tree.isRoot(node)){
            return null;
        } else {
            NodeRef parent = tree.getParent(node);
            for(int i=0; i<tree.getChildCount(parent); i++){
                if(tree.getChild(parent,i)!=node){
                    return tree.getChild(parent,i);
                }
            }
        }
        return null;
    }

    //Return a set of nodes that are not descendants of (or equal to) the current node and have the same painting as it.
    //The node recalculation is going to need reworking once the time comes to test the extended version,

    public HashSet<Integer> samePaintingUpTree(NodeRef node, boolean flagForRecalc){
        if(!nodeRecalculationNeeded[node.getNumber()] && flagForRecalc){
            extendedflagForRecalculation(virusTree, virusTree.getChild(node,0), nodeRecalculationNeeded);
            extendedflagForRecalculation(virusTree, virusTree.getChild(node,1), nodeRecalculationNeeded);
        }
        HashSet<Integer> out = new HashSet<Integer>();
        AbstractCase painting = branchMap[node.getNumber()];
        NodeRef ancestorNode = virusTree.getParent(node);
        while(branchMap[ancestorNode.getNumber()]==painting){
            out.add(ancestorNode.getNumber());
            if(extended){
                if(countChildrenWithSamePainting(ancestorNode)==2){
                    NodeRef otherChild = sibling(virusTree, ancestorNode);
                    out.add(otherChild.getNumber());
                    out.addAll(samePaintingDownTree(otherChild, flagForRecalc));
                }
            }
            ancestorNode = virusTree.getParent(ancestorNode);
        }
        return out;
    }

    //Return a set of nodes that are descendants (and not equal to) the current node and have the same painting as it.


    public HashSet<Integer> samePaintingDownTree(NodeRef node, boolean flagForRecalc){
        HashSet<Integer> out = new HashSet<Integer>();
        AbstractCase painting = branchMap[node.getNumber()];
        boolean creepsFurther = false;
        for(int i=0; i<virusTree.getChildCount(node); i++){
            if(branchMap[virusTree.getChild(node,i).getNumber()]==painting){
                creepsFurther = true;
                out.add(virusTree.getChild(node,i).getNumber());
                out.addAll(samePaintingDownTree(virusTree.getChild(node,i), flagForRecalc));
            }
        }
        if(flagForRecalc && !creepsFurther){
            extendedflagForRecalculation(virusTree, virusTree.getChild(node,0), nodeRecalculationNeeded);
            extendedflagForRecalculation(virusTree, virusTree.getChild(node,1), nodeRecalculationNeeded);
        }
        return out;
    }

    // Return the node numbers of the entire subtree with the same painting as this node (including itself)

    public HashSet<Integer> samePainting(NodeRef node, boolean flagForRecalc){
        HashSet<Integer> out = new HashSet<Integer>();
        out.add(node.getNumber());
        out.addAll(samePaintingDownTree(node, flagForRecalc));
        out.addAll(samePaintingUpTree(node, flagForRecalc));
        return out;
    }

    // The flags indicate if a node's painting has changed. Only tree operators can use this shortcut; if parameters of
    // the epidemiological model have changed then the whole tree's likelihood needs recalculating. Somewhat out of
    // date at this point.

    public static void extendedflagForRecalculation(TreeModel tree, NodeRef node, boolean[] flags){
        flags[node.getNumber()]=true;
        for(int i=0; i<tree.getChildCount(node); i++){
            flags[tree.getChild(node,i).getNumber()]=true;
        }
        NodeRef currentNode=node;
        while(!tree.isRoot(currentNode) && !flags[currentNode.getNumber()]){
            currentNode = tree.getParent(currentNode);
            flags[currentNode.getNumber()]=true;
        }
    }

    public static void flagForRecalculation(NodeRef node, boolean[] flags){
        flags[node.getNumber()]=true;
    }

    public void changeCreepLock(int index, boolean value){
        creepLocks[index]=value;
    }

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************


    protected final void handleModelChangedEvent(Model model, Object object, int index) {
        Arrays.fill(nodeRecalculationNeeded, true);
        likelihoodKnown = false;
    }


    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************


    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        Arrays.fill(nodeRecalculationNeeded, true);
        likelihoodKnown = false;
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    /**
     * Stores the precalculated state (in this case the node labels and subtree likelihoods)
     */
    protected final void storeState() {
        storedBranchMap = Arrays.copyOf(branchMap, branchMap.length);
        storedNodeLogLikelihoods = Arrays.copyOf(nodeLogLikelihoods, nodeLogLikelihoods.length);
        storedRecalculationArray = Arrays.copyOf(nodeRecalculationNeeded, nodeRecalculationNeeded.length);
    }

    /**
     * Restores the precalculated state.
     */
    protected final void restoreState() {
        branchMap = storedBranchMap;
        nodeLogLikelihoods = storedNodeLogLikelihoods;
        nodeRecalculationNeeded = storedRecalculationArray;
        likelihoodKnown = false;
    }

    protected final void acceptState() {
    } // nothing to do

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    public final Model getModel() {
        return this;
    }

    public final AbstractCase[] getBranchMap(){
        return branchMap;
    }

    public final void setBranchMap(AbstractCase[] map){
        branchMap = map;
    }

    public final TreeModel getTree(){
        return virusTree;
    }

    public final double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    public final void makeDirty() {
        // so every switch operator is forcing recalculation of the whole tree? Can't be necessary.
        Arrays.fill(nodeRecalculationNeeded, true);
        if(extended){
            recalculateLocks();
        }
        likelihoodKnown = false;
    }

    public void makeDirty(boolean cleanTree){
        if(cleanTree){
            makeDirty();
        } else {
            if(extended){
                recalculateLocks();
            }
            likelihoodKnown = false;
        }
    }

    /**
     * Calculates the log likelihood of this set of node labels given the tree.
     */
    public double calculateLogLikelihood() {
        return calculateLogLikelihood(branchMap);
    }

    public double calculateLogLikelihood(AbstractCase[] map){
        // temporary stuff

/*
        treeInfectionDates = new HashMap<AbstractCase,Double>();
        treeInfectionDays = new HashMap<AbstractCase,Integer>();
        treeExamDates = new HashMap<AbstractCase, Double>();
        treeExamDays = new HashMap<AbstractCase, Integer>();
        dataExamDates = new HashMap<AbstractCase, Double>();
        dataInfectiousnessDates = new HashMap<AbstractCase, Double>();
*/


        NodeRef root = virusTree.getRoot();

/*        buildDataTable(root,map,true,true,"test_table.csv");*/

        return calculateNodeTransmissionLogLikelihood(root, Integer.MIN_VALUE, map);
    }

    /* Return the integer day on which the given node occurred */

    public int getNodeDay(NodeRef node){
        Double nodeHeight = getHeight(node);
        Double estimatedNodeTime = estimatedLastSampleTime-nodeHeight;
        Date nodeDate = new Date(estimatedNodeTime, Units.Type.DAYS, false);
        /*Since day t begins at time t-1, to find the day of a given decimal time value, take the ceiling.*/
        return (int)Math.ceil(nodeDate.getTimeValue());
    }

    /**
     * Given a node, calculates the log likelihood of its parent branch and then goes down the tree and calculates the
     * log likelihoods of lower branches.
     */

    private double calculateNodeTransmissionLogLikelihood(NodeRef node, int parentDay, AbstractCase[]
            currentBranchMap) {
        double logLikelihood=0;
        AbstractCase nodeCase = currentBranchMap[node.getNumber()];
        Integer nodeDay = getNodeDay(node);
        /* Likelihood of the branch leading to the node.*/
        if(!nodeRecalculationNeeded[node.getNumber()]){
            logLikelihood = logLikelihood + nodeLogLikelihoods[node.getNumber()];
        } else {
            double nodeLogLikelihood = 0;
            if(virusTree.isRoot(node)){
                /*Likelihood of root node case being infectious by the time of the root.*/
                nodeLogLikelihood =  Math.log(cases.noTransmissionBranchLikelihood(nodeCase, nodeDay));
            } else {
                /*Likelihood of this case being infected at the time of the parent node and infectious by the time
                of the child node.*/
                NodeRef parent = virusTree.getParent(node);
                nodeLogLikelihood = Math.log(cases.transmissionBranchLikelihood(currentBranchMap[parent.getNumber()],
                        nodeCase, parentDay, nodeDay));
            }
            nodeLogLikelihoods[node.getNumber()]=nodeLogLikelihood;
            nodeRecalculationNeeded[node.getNumber()]=false;
            logLikelihood += nodeLogLikelihood;
        }
        if(!virusTree.isExternal(node)){
            for(int i=0; i<virusTree.getChildCount(node); i++){
                logLikelihood = logLikelihood + calculateNodeTransmissionLogLikelihood(virusTree.getChild(node,i),
                        nodeDay, currentBranchMap);
            }
        }
        return logLikelihood;
    }

    public boolean[] getRecalculationArray(){
        return nodeRecalculationNeeded;
    }

    private double getHeight(NodeRef node){
        return virusTree.getNodeHeight(node);
    }

    /* Writes a table of farm infection, infectiousness and exam dates based on the current painting */

/*
    private void buildDataTable(NodeRef node, AbstractCase[] currentBranchMap, boolean firstLoop,
                                boolean doOutput, String fileName){
        AbstractCase nodeCase = currentBranchMap[node.getNumber()];
        if(virusTree.isRoot(node)){
            treeInfectionDates.put(nodeCase, -1.0);
            treeInfectionDays.put(nodeCase, -1);
        } else {
            NodeRef parent = virusTree.getParent(node);
            AbstractCase parentCase = currentBranchMap[parent.getNumber()];
            if(!parentCase.equals(nodeCase)){
                treeInfectionDates.put(nodeCase, getHeight(parent));
                treeInfectionDays.put(nodeCase, getNodeDay(parent));
            }
        }
        if(virusTree.isExternal(node)){
            treeExamDates.put(nodeCase, getHeight(node));
            treeExamDays.put(nodeCase, getNodeDay(node));
        } else {
            for(int i=0;i<2;i++){
                buildDataTable(virusTree.getChild(node,i),currentBranchMap,false,false,null);
            }
        }
        if(firstLoop){
            for(AbstractCase currentCase: cases.getCases()){
                dataExamDates.put(currentCase,
                        currentCase.getEndOfInfectiousDateModeHeight(getLatestTaxonDate(virusTree)));
                dataInfectiousnessDates.put(currentCase,
                        currentCase.getInfectiousDateModeHeight(getLatestTaxonDate(virusTree)));
            }
            if(doOutput){
                try{
                    BufferedWriter writer = new BufferedWriter(new FileWriter(fileName));
                    String header = "Name,tree_infection_date,tree_infection_day,tree_exam_date,tree_exam_day,"+
                            "data_exam_date,data_infectious_date";
                    writer.write(header);
                    writer.newLine();
                    for(AbstractCase currentCase: cases.getCases()){
                        String line = currentCase.getName() + ","
                                + Double.toString(treeInfectionDates.get(currentCase)) + "," +
                                Integer.toString(treeInfectionDays.get(currentCase)) + "," +
                                Double.toString(treeExamDates.get(currentCase)) + "," +
                                Integer.toString(treeExamDays.get(currentCase)) + "," +
                                Double.toString(dataExamDates.get(currentCase)) + "," +
                                Double.toString(dataInfectiousnessDates.get(currentCase));
                        writer.write(line);
                        writer.newLine();
                    }
                    writer.close();
                } catch(IOException e){
                    System.out.println("Write to file failed (IOException");
                }
            }
        }
    }
*/

        /* Return the case whose infection resulted from this node (if internal) */

    private AbstractCase getInfected(NodeRef node){
        if(virusTree.isExternal(node)){
            throw new RuntimeException("Node is external");
        } else {
            for(int i=0;i<2;i++){
                NodeRef child = virusTree.getChild(node, i);
                if(branchMap[child.getNumber()]!=branchMap[node.getNumber()]){
                    return branchMap[child.getNumber()];
                }
            }
        }
        throw new RuntimeException("Node appears to have two children with its infector");
    }

    /* Return the case which infected this case */

    public AbstractCase getInfector(AbstractCase thisCase){
        NodeRef tip = tipMap.get(thisCase);
        return getInfector(tip);
    }

    /* Return the case which was the infector in the infection event represented by this node */

    public AbstractCase getInfector(NodeRef node){
        if(virusTree.isRoot(node) || node.getNumber()==virusTree.getRoot().getNumber()){
            return null;
        } else {
            AbstractCase nodeCase = branchMap[node.getNumber()];
            if(branchMap[virusTree.getParent(node).getNumber()]!=nodeCase){
                return branchMap[virusTree.getParent(node).getNumber()];
            } else {
                return getInfector(virusTree.getParent(node));
            }
        }
    }

    public boolean checkPaintingIntegrity(boolean verbose){
        boolean out=true;
        for(int i=0; i<virusTree.getInternalNodeCount(); i++){
            NodeRef node = virusTree.getInternalNode(i);
            NodeRef firstChild = virusTree.getChild(node,0);
            NodeRef secondChild = virusTree.getChild(node,1);
            NodeRef parent = virusTree.getParent(node);
            if(branchMap[node.getNumber()]!=branchMap[firstChild.getNumber()] &&
                    branchMap[node.getNumber()]!=branchMap[secondChild.getNumber()] &&
                    (!extended || branchMap[node.getNumber()]!=branchMap[parent.getNumber()])){
                out = false;
                if(!verbose){
                    break;
                } else {
                    System.out.println("Node "+node.getNumber()+" failed painting integrity check:");
                    System.out.println("Node painting: "+branchMap[node.getNumber()].getName());
                    System.out.println("Parent painting: "+branchMap[parent.getNumber()].getName());
                    System.out.println("Child 1 painting: "+branchMap[firstChild.getNumber()].getName());
                    System.out.println("Child 2 painting: "+branchMap[secondChild.getNumber()].getName());
                    System.out.println();
                }
            }
        }
        return out;
    }


    /* Return the case (painting) of the parent of this node */

    public AbstractCase getParentCase(NodeRef node){
        return branchMap[virusTree.getParent(node).getNumber()];
    }

    /* Populates the branch map for external nodes */

    private AbstractCase[] prepareExternalNodeMap(AbstractCase[] map){
        for(int i=0; i<virusTree.getExternalNodeCount(); i++){
            TreeModel.Node currentExternalNode = (TreeModel.Node)virusTree.getExternalNode(i);
            Taxon currentTaxon = currentExternalNode.taxon;
            for(AbstractCase thisCase : cases.getCases()){
                for(Taxon caseTaxon: thisCase.getAssociatedTaxa()){
                    if(caseTaxon.equals(currentTaxon)){
                        map[currentExternalNode.getNumber()]=thisCase;
                    }
                }
            }
        }
        return map;
    }

/*  The CSV file should have a header line, and then a line for each node with its id in the first column and the id
    of its parent in the second. The node with no parent has "Start" in the second column.*/

    private AbstractCase[] paintSpecificNetwork(String networkFileName){
        System.out.println("Using specified starting network.");
        try{
            BufferedReader reader = new BufferedReader (new FileReader(networkFileName));
            HashMap<AbstractCase, AbstractCase> specificParentMap = new HashMap<AbstractCase, AbstractCase>();
            // skip header line
            reader.readLine();
            String currentLine = reader.readLine();
            while(currentLine!=null){
                currentLine = currentLine.replace("\"", "");
                String[] splitLine = currentLine.split("\\,");
                if(!splitLine[1].equals("Start")){
                    specificParentMap.put(cases.getCase(splitLine[0]), cases.getCase(splitLine[1]));
                } else {
                    specificParentMap.put(cases.getCase(splitLine[0]),null);
                }
                currentLine = reader.readLine();
            }
            reader.close();
            return paintSpecificNetwork(specificParentMap);
        } catch(IOException e){
            System.out.println("Cannot read file: " + networkFileName );
            return null;
        }
    }


    /* Takes a HashMap referring each case to its parent, and tries to paint the tree with it. This only works on
     * the non-extended version right now, watch it. */

    private AbstractCase[] paintSpecificNetwork(HashMap<AbstractCase, AbstractCase> map){
        Arrays.fill(nodeRecalculationNeeded,true);
        AbstractCase[] branchArray = new AbstractCase[virusTree.getNodeCount()];
        branchArray = prepareExternalNodeMap(branchArray);
        TreeModel.Node root = (TreeModel.Node)virusTree.getRoot();
        specificallyPaintNode(root, branchArray, map);
        return branchArray;
    }


    /* Paints a phylogenetic tree node and its children according to a specified map of child to parent cases */

    private AbstractCase specificallyPaintNode(TreeModel.Node node, AbstractCase[] map,
                                               HashMap<AbstractCase,AbstractCase> parents){
        if(node.isExternal()){
            return map[node.getNumber()];
        } else {
            AbstractCase[] childPaintings = new AbstractCase[2];
            for(int i=0; i<node.getChildCount(); i++){
                childPaintings[i] = specificallyPaintNode(node.getChild(i), map, parents);
            }
            if(parents.get(childPaintings[1])==childPaintings[0]){
                map[node.getNumber()]=childPaintings[0];
            } else if(parents.get(childPaintings[0])==childPaintings[1]){
                map[node.getNumber()]=childPaintings[1];
            } else {
                throw new RuntimeException("This network does not appear to be compatible with the tree");
            }
            return map[node.getNumber()];
        }
    }

    public void writeNetworkToFile(String fileName){
        try{
            BufferedWriter writer = new BufferedWriter(new FileWriter(fileName));
            writer.write("Case,Parent");
            writer.newLine();
            for(int i=0; i<virusTree.getExternalNodeCount(); i++){
                TreeModel.Node extNode = (TreeModel.Node)virusTree.getExternalNode(i);
                String tipName = extNode.taxon.toString();
                String infector;
                try{
                    infector = getInfector(extNode).getName();
                } catch(NullPointerException e){
                    infector = "Start";
                }
                writer.write(tipName + "," + infector);
                writer.newLine();
            }
            writer.close();
        } catch(IOException e) {
            System.out.println("Failed to write to file");
        }

    }

    /*Given a new tree with no labels, associates each of the terminal branches with the relevant case and then
  * generates a random painting of the rest of the tree to start off with. If checkNonZero is true in
  * randomlyPaintNode then the network will be checked to prohibit links with zero (or rounded to zero)
  * likelihood first.*/

    private AbstractCase[] paintRandomNetwork(){
        boolean gotOne = false;
        AbstractCase[] map = null;
        int tries = 1;
        System.out.println("Generating a random starting painting of the tree (checking nonzero likelihood for all " +
                "branches and repeating up to 100 times until a start with nonzero likelihood is found)");
        System.out.print("Attempt: ");
        while(!gotOne){
            System.out.print(tries + "...");
            map = prepareExternalNodeMap(new AbstractCase[virusTree.getNodeCount()]);
            //Warning - if the BadPaintingException in randomlyPaintNode might be caused by a bug rather than both
            //likelihoods rounding to zero, you want to stop catching this to investigate.
            try{
                paintRandomNetwork(map,true);
                if(calculateLogLikelihood(map)!=Double.NEGATIVE_INFINITY){
                    gotOne = true;
                    System.out.print("found.");
                }
            } catch(BadPaintingException ignored){}
            tries++;
            if(tries==101){
                System.out.println("giving " +
                        "up.");
                break;
            }
        }
        System.out.println();
        return map;
    }



    /* Paints a phylogenetic tree with a random compatible painting; if checkNonZero is true, make sure all branch
   likelihoods are nonzero in the process (this sometimes still results in a zero likelihood for the whole tree, but
   is much less likely to).
    */

    private AbstractCase[] paintRandomNetwork(AbstractCase[] map, boolean checkNonZero){
        Arrays.fill(nodeRecalculationNeeded,true);
        TreeModel.Node root = (TreeModel.Node)virusTree.getRoot();
        randomlyPaintNode(root, map, checkNonZero);
        return map;
    }

    private AbstractCase randomlyPaintNode(TreeModel.Node node, AbstractCase[] map, boolean checkNonZero){
        if(node.isExternal()){
            return map[node.getNumber()];
        } else {
            AbstractCase[] choices = new AbstractCase[2];
            for(int i=0; i<node.getChildCount(); i++){
                if((map[node.getChild(i).getNumber()]==null)){
                    choices[i] = randomlyPaintNode(node.getChild(i),map,checkNonZero);
                } else {
                    choices[i] = map[node.getChild(i).getNumber()];
                }
            }
            int randomSelection = MathUtils.nextInt(2);
            AbstractCase decision;
            if(checkNonZero){
                Double[] branchLogLs = new Double[2];
                for(int i=0; i<2; i++){
                    branchLogLs[i]= cases.transmissionBranchLogLikelihood(choices[i], choices[1 - i], getNodeDay(node),
                            getNodeDay(node.getChild(1 - i)));
                }
                if(branchLogLs[0]==Double.NEGATIVE_INFINITY && branchLogLs[1]==Double.NEGATIVE_INFINITY){
                    throw new BadPaintingException("Both branch possibilities have zero likelihood: "
                            +node.toString()+", cases " + choices[0].getName() + " and " + choices[1].getName() + ".");
                } else if(branchLogLs[0]==Double.NEGATIVE_INFINITY || branchLogLs[1]==Double.NEGATIVE_INFINITY){
                    if(branchLogLs[0]==Double.NEGATIVE_INFINITY){
                        decision = choices[1];
                    } else {
                        decision = choices[0];
                    }
                } else {
                    decision = choices[randomSelection];
                }
            } else {
                decision = choices[randomSelection];
            }
            map[node.getNumber()]=decision;
            return decision;
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
    public static final String VIRUS_TREE = "virusTree";
    public static final String STARTING_NETWORK = "startingNetwork";

    public String getParserName() {
        return CASE_TO_CASE_TRANSMISSION_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        FlexibleTree flexTree = (FlexibleTree) xo.getElementFirstChild(VIRUS_TREE);

        String startingNetworkFileName=null;

        if(xo.hasChildNamed(STARTING_NETWORK)){
            startingNetworkFileName = (String) xo.getElementFirstChild(STARTING_NETWORK);
        }

        TreeModel virusTree = new TreeModel(flexTree);

        AbstractOutbreak caseSet = (AbstractOutbreak) xo.getChild(AbstractOutbreak.class);

        CaseToCaseTransmissionLikelihood likelihood;

        final boolean extended = xo.getBooleanAttribute("extended");

        try {
            likelihood = new CaseToCaseTransmissionLikelihood(virusTree, caseSet, startingNetworkFileName,
                    extended);
        } catch (TaxonList.MissingTaxonException e) {
            throw new XMLParseException(e.toString());
        }

        return likelihood;
    }

    public String getParserDescription() {
        return "This element represents a likelihood function for case to case transmission.";
    }

    public Class getReturnType() {
        return CaseToCaseTransmissionLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule("extended"),
            new ElementRule("virusTree", Tree.class, "The tree"),
            new ElementRule(AbstractOutbreak.class, "The set of cases"),
            new ElementRule("startingNetwork", String.class, "A CSV file containing a specified starting network",
                    true)
    };
};


    //************************************************************************
    // Loggable implementation
    //************************************************************************

    public LogColumn[] getColumns(){

        LogColumn[] columns = new LogColumn[cases.size()];

        for(int i=0; i< cases.size(); i++){
            final AbstractCase infected = cases.getCase(i);
            columns[i] = new LogColumn.Abstract(infected.toString()){
                protected String getFormattedValue() {
                    if(getInfector(infected)==null){
                        return "Start";
                    } else {
                        return getInfector(infected).toString();
                    }
                }
            };
        }

        return columns;

    }

    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(new Citation(new Author[]{new Author("M", "Hall"), new Author("A", "Rambaut")},
                Citation.Status.IN_PREPARATION));
        return citations;
    }

public class BadPaintingException extends RuntimeException{

    public BadPaintingException(String s){
        super(s);
    }

}

}
