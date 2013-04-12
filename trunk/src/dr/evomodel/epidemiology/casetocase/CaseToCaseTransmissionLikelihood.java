package dr.evomodel.epidemiology.casetocase;
import dr.evolution.tree.*;
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
 * @author Andrew Rambaut
 * @version $Id: $
 */
public class CaseToCaseTransmissionLikelihood extends AbstractModelLikelihood implements Loggable, Citable,
        TreeTraitProvider {

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
    protected TreeTraitProvider.Helper treeTraits = new Helper();

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

    // for Felsenstein version:

    private boolean felsenstein;
    private double[] tripletLikelihoods;
    private double[] storedTripletLikelihoods;
    private boolean currentReconstructionExists;

    // PUBLIC STUFF

    // Name

    public static final String CASE_TO_CASE_TRANSMISSION_LIKELIHOOD = "caseToCaseTransmissionLikelihood";
    public static final String PAINTINGS_KEY = "paintings";

    // Basic constructor.

    public CaseToCaseTransmissionLikelihood(TreeModel virusTree, AbstractOutbreak caseData,
                                            String startingNetworkFileName, boolean extended, boolean felsenstein)
            throws TaxonList.MissingTaxonException {
        this(CASE_TO_CASE_TRANSMISSION_LIKELIHOOD, virusTree, caseData, startingNetworkFileName, extended, felsenstein);
    }

    // Legacy constructor

    public CaseToCaseTransmissionLikelihood(String name, TreeModel virusTree, AbstractOutbreak caseData, String
            startingNetworkFileName){
        this(name, virusTree, caseData, startingNetworkFileName, false, false);
    }


    // Constructor for an instance with a non-default name

    public CaseToCaseTransmissionLikelihood(String name, TreeModel virusTree, AbstractOutbreak caseData, String
            startingNetworkFileName, boolean extended, boolean felsenstein) {

        super(name);
        this.virusTree = virusTree;
        addModel(virusTree);
        Date lastSampleDate = getLatestTaxonDate(virusTree);

        /* We assume samples were taken at the end of the date of the last tip */

        estimatedLastSampleTime = lastSampleDate.getTimeValue();
        cases = caseData;
        this.extended = extended;
        this.felsenstein = felsenstein;
        addModel(cases);
        verbose = false;

        tipMap = new HashMap<AbstractCase, NodeRef>();

        //map the cases to the external nodes
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

        if(felsenstein){
            tripletLikelihoods = new double[virusTree.getInternalNodeCount()
                    *virusTree.getExternalNodeCount()
                    *virusTree.getExternalNodeCount()
                    *virusTree.getExternalNodeCount()];
            currentReconstructionExists = false;
            branchMap = new AbstractCase[virusTree.getNodeCount()];
        } else {
            nodeRecalculationNeeded = new boolean[virusTree.getNodeCount()];
            Arrays.fill(nodeRecalculationNeeded, true);
            nodeLogLikelihoods = new double[virusTree.getNodeCount()];

            //paint the starting network onto the tree

            if(startingNetworkFileName==null){
                branchMap = paintRandomNetwork();
            } else {
                branchMap = paintSpecificNetwork(startingNetworkFileName);
            }

        }
        treeTraits.addTrait(PAINTINGS_KEY, new TreeTrait.S() {
            public String getTraitName() {
                return "host_case";
            }

            public Intent getIntent() {
                return Intent.NODE;
            }

            public String getTrait(Tree tree, NodeRef node) {
                return getPaintingForNode(tree, node);
            }


        });

        //@todo deal with felsenstein==true and extended==true (probably by prohibiting it)
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

    private NodeRef[] getChildren(NodeRef node){
        NodeRef[] children = new NodeRef[virusTree.getChildCount(node)];
        for(int i=0; i<virusTree.getChildCount(node); i++){
            children[i] = virusTree.getChild(node,i);
        }
        return children;
    }

    private int[] getChildDays(NodeRef node){
        int[] childDays = new int[virusTree.getChildCount(node)];
        for(int i=0; i<virusTree.getChildCount(node); i++){
            childDays[i] = getNodeDay(virusTree.getChild(node,i));
        }
        return childDays;
    }


    // ************************************************************************************
    // EXTENDED VERSION METHODS
    // ************************************************************************************

    /* check if the given node is tip-linked under the current painting (the tip corresponding to its painting is
    a descendant of it
     */

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
        if(!felsenstein){
            Arrays.fill(nodeRecalculationNeeded, true);
        } else {
            currentReconstructionExists = false;
        }
        likelihoodKnown = false;
    }


    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************


    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if(!felsenstein){
            Arrays.fill(nodeRecalculationNeeded, true);
        } else {
            currentReconstructionExists = false;
        }
        likelihoodKnown = false;
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    /**
     * Stores the precalculated state (in this case the node labels and subtree likelihoods)
     */
    protected final void storeState() {
        if(!felsenstein){
            storedBranchMap = Arrays.copyOf(branchMap, branchMap.length);
            storedNodeLogLikelihoods = Arrays.copyOf(nodeLogLikelihoods, nodeLogLikelihoods.length);
            storedRecalculationArray = Arrays.copyOf(nodeRecalculationNeeded, nodeRecalculationNeeded.length);
        } else {
            storedTripletLikelihoods = Arrays.copyOf(tripletLikelihoods, tripletLikelihoods.length);

        }
    }

    /**
     * Restores the precalculated state.
     */
    protected final void restoreState() {
        if(!felsenstein){
            branchMap = storedBranchMap;
            nodeLogLikelihoods = storedNodeLogLikelihoods;
            nodeRecalculationNeeded = storedRecalculationArray;
        } else {
            tripletLikelihoods = storedTripletLikelihoods;
        }
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
        if(!felsenstein){
            Arrays.fill(nodeRecalculationNeeded, true);
        } else {
            currentReconstructionExists = false;
        }
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
        if(felsenstein){
            NodeRef root = virusTree.getRoot();
            double totalLikelihood = 0;
            for(AbstractCase thisCase: cases.getCases()){
                totalLikelihood += L(root,thisCase);
            }
            return Math.log(totalLikelihood);
        } else {
            if(!checkPaintingIntegrity(true)){
                throw new RuntimeException("Not a painting");
            }
            return calculateLogLikelihood(branchMap);
        }
    }

    public double calculateLogLikelihood(AbstractCase[] map){
        NodeRef root = virusTree.getRoot();
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
        return getInfector(thisCase, branchMap);
    }

    /* Return the case which was the infector in the infection event represented by this node */

    public AbstractCase getInfector(NodeRef node){
        return getInfector(node, branchMap);
    }

    public AbstractCase getInfector(AbstractCase thisCase, AbstractCase[] branchMap){
        NodeRef tip = tipMap.get(thisCase);
        return getInfector(tip, branchMap);
    }

    public AbstractCase getInfector(NodeRef node, AbstractCase[] branchMap){
        if(virusTree.isRoot(node) || node.getNumber()==virusTree.getRoot().getNumber()){
            return null;
        } else {
            AbstractCase nodeCase = branchMap[node.getNumber()];
            if(branchMap[virusTree.getParent(node).getNumber()]!=nodeCase){
                return branchMap[virusTree.getParent(node).getNumber()];
            } else {
                return getInfector(virusTree.getParent(node), branchMap);
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
                throw new RuntimeException("Failed to find a starting transmission network with nonzero likelihood");
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
    // Felsenstein version methods
    //************************************************************************


    /* Check whether the tree will admit a particular painting for a particular node. This is true if a
     * descendant external node is the one corresponding to the case */

    private boolean treeWillAdmitCaseHere (NodeRef node, AbstractCase thisCase){
        if(virusTree.isExternal(node)){
            return tipMap.get(thisCase).toString().equals(node.toString());
        } else {
            boolean check1 = treeWillAdmitCaseHere(virusTree.getChild(node,0), thisCase);
            boolean check2 = treeWillAdmitCaseHere(virusTree.getChild(node,1), thisCase);
            return check1 || check2;
        }
    }

    /* The likelihood of observing a case i at node alpha*/

    private double L(NodeRef alpha, AbstractCase i){
        if(virusTree.isExternal(alpha)){
            return tipMap.get(i).toString().equals(alpha.toString()) ? 1 : 0;
        } else {
            double totalLikelihood = 0;
            NodeRef[] children = getChildren(alpha);
            for(AbstractCase case_x: cases.getCases()){
                for(AbstractCase case_y: cases.getCases()){
                    double P_term = P(i, case_x, case_y, alpha);
                    double L_term = 0;
                    int noTips = virusTree.getExternalNodeCount();
                    if(P_term!=0){
                        L_term = L(children[0],case_x)*L(children[1],case_y);
                    }
                    tripletLikelihoods[(alpha.getNumber()-noTips)*noTips*noTips*noTips
                            + cases.getCases().indexOf(i)*noTips*noTips
                            + cases.getCases().indexOf(case_x)*noTips
                            + cases.getCases().indexOf(case_y)] = L_term;

                    totalLikelihood = totalLikelihood + L_term*P_term;
                }
            }
            return totalLikelihood;
        }
    }

    /* The probability that you observe an a at alpha given a b at one child and a c at the other*/

    private double P(AbstractCase a, AbstractCase b, AbstractCase c, NodeRef alpha){
        // Internal nodes only
        if(virusTree.isExternal(alpha)){
            throw new RuntimeException("alpha must be external");
        }
        // this is a standard painting, so a must be either b or c but not both
        if((a!=b && a!=c) || (b==c)){
            return 0;
        }
        NodeRef[] children = getChildren(alpha);
        //Tree must admit the cases in these positions
        if(!treeWillAdmitCaseHere(alpha, a)
                || !treeWillAdmitCaseHere(children[0],b)
                || !treeWillAdmitCaseHere(children[1],c)) {
            return 0;
        }
        int alphaDay = getNodeDay(alpha);
        int[] childDays = getChildDays(alpha);
        if (a==b){
            double c_term = cases.transmissionBranchLikelihood(a, c, alphaDay, childDays[1]);
            double a_term = cases.noTransmissionBranchLikelihood(a, alphaDay)
                    /cases.noTransmissionBranchLikelihood(a, childDays[0]);
            if(!Double.isNaN(a_term)){
                return c_term*a_term;
            } else {
                return 0;
            }
        } else {
            double b_term = cases.transmissionBranchLikelihood(a, b, alphaDay, childDays[0]);
            double a_term = cases.noTransmissionBranchLikelihood(a, alphaDay)
                    /cases.noTransmissionBranchLikelihood(a, childDays[1]);
            if(!Double.isNaN(a_term)){
                return b_term*a_term;
            } else {
                return 0;
            }
        }
    }

    /* Sample a random transmission tree from the likelihoods for the paintings of each triplet of nodes */

    private AbstractCase[] sampleTransmissionTree(){
        AbstractCase[] samplePainting = new AbstractCase[virusTree.getNodeCount()];
        for(int i=0; i<virusTree.getExternalNodeCount(); i++){
            TreeModel.Node tip = (TreeModel.Node)virusTree.getExternalNode(i);
            samplePainting[tip.getNumber()]=cases.getCase(tip.taxon.toString());
        }
        NodeRef root = virusTree.getRoot();
        int noTips = virusTree.getExternalNodeCount();
        double[] rootLikelihoods = Arrays.copyOfRange(tripletLikelihoods,
                (root.getNumber()-noTips)*noTips*noTips*noTips,
                (root.getNumber()-noTips+1)*noTips*noTips*noTips);
        int dim = virusTree.getExternalNodeCount();
        int choice = MathUtils.randomChoicePDF(rootLikelihoods);
        int[] choices = new int[3];
        choices[0] = (choice - (choice % (dim*dim)))/(dim*dim);
        choices[1] = ((choice % (dim*dim)) - ((choice % dim*dim) % dim)) / dim;
        choices[2] = (choice % (dim*dim)) % dim;
        samplePainting[root.getNumber()]=cases.getCase(choices[0]);
        samplePainting[virusTree.getChild(root,0).getNumber()]=cases.getCase(choices[1]);
        samplePainting[virusTree.getChild(root,1).getNumber()]=cases.getCase(choices[2]);

        for(int i=0; i<2; i++){
            fillUp(virusTree.getChild(root,i),samplePainting);
        }

        return samplePainting;
    }


    private void fillUp(NodeRef node, AbstractCase[] painting){

        if(!virusTree.isExternal(node)){
            int dim = virusTree.getExternalNodeCount();
            int nodePaintingNumber = cases.getCases().indexOf(painting[node.getNumber()]);
            double[] childLikelihoods = Arrays.copyOfRange(tripletLikelihoods,
                    (node.getNumber()-dim)*dim*dim*dim + nodePaintingNumber*dim*dim,
                    (node.getNumber()-dim)*dim*dim*dim + (nodePaintingNumber+1)*dim*dim);
            int choice = MathUtils.randomChoicePDF(childLikelihoods);
            int[] choices = new int[2];
            choices[0] = (choice - (choice % dim))/dim;
            choices[1] = choice % dim;
            for(int i=0; i<2; i++){
                NodeRef child = virusTree.getChild(node,i);
                painting[child.getNumber()]=cases.getCase(choices[i]);
                if(!virusTree.isExternal(child)){
                    fillUp(child, painting);
                }
            }
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

            TreeModel virusTree = (TreeModel) xo.getElementFirstChild(VIRUS_TREE);

            String startingNetworkFileName=null;

            if(xo.hasChildNamed(STARTING_NETWORK)){
                startingNetworkFileName = (String) xo.getElementFirstChild(STARTING_NETWORK);
            }


            AbstractOutbreak caseSet = (AbstractOutbreak) xo.getChild(AbstractOutbreak.class);

            CaseToCaseTransmissionLikelihood likelihood;

            final boolean extended = xo.getBooleanAttribute("extended");
            final boolean felsenstein = xo.getBooleanAttribute("felsenstein");

            try {
                likelihood = new CaseToCaseTransmissionLikelihood(virusTree, caseSet, startingNetworkFileName,
                        extended, felsenstein);
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
                AttributeRule.newBooleanRule("felsenstein"),
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

        if(felsenstein){

            for(int i=0; i<cases.size(); i++){
                final AbstractCase infected = cases.getCase(i);
                columns[i] = new LogColumn.Abstract(infected.toString()){
                    @Override
                    protected String getFormattedValue() {
                        if(!currentReconstructionExists){
                            branchMap = sampleTransmissionTree();
                        }
                        if(getInfector(infected, branchMap)==null){
                            return "Start";
                        } else {
                            return getInfector(infected, branchMap).toString();
                        }
                    }
                };
            }

            return columns;

        } else {

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

    }

    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(new Citation(new Author[]{new Author("M", "Hall"), new Author("A", "Rambaut")},
                Citation.Status.IN_PREPARATION));
        return citations;
    }

    // **************************************************************
    // TreeTraitProvider IMPLEMENTATION
    // **************************************************************

    public TreeTrait[] getTreeTraits() {
        return treeTraits.getTreeTraits();
    }

    public TreeTrait getTreeTrait(String key) {
        return treeTraits.getTreeTrait(key);
    }

    public String getPaintingForNode(Tree tree, NodeRef node) {
        if (tree != virusTree) {
            throw new RuntimeException("Can only reconstruct states on treeModel given to constructor");
        }

        if (!likelihoodKnown) {
            calculateLogLikelihood();
            likelihoodKnown = true;
        }

        if(felsenstein && !currentReconstructionExists){
            branchMap = sampleTransmissionTree();
        }

        return branchMap[node.getNumber()].toString();

    }



    public class BadPaintingException extends RuntimeException{

        public BadPaintingException(String s){
            super(s);
        }
    }
}




