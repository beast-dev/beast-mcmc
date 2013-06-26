package dr.evomodel.epidemiology.casetocase;
import dr.app.tools.NexusExporter;
import dr.evolution.tree.*;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.LikelihoodCore;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.LogTricks;
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
 * Timescale must be in days. Python scripts to write XML for it and analyse the posterior set of networks exist; contact MH.
 *
 * @author Matthew Hall
 * @author Andrew Rambaut
 * @version $Id: $
 */
public class CaseToCaseTransmissionLikelihood extends AbstractModelLikelihood implements Loggable, Citable,
        TreeTraitProvider {

    /* The phylogenetic tree. */

    private TreeModel virusTree;
    private int noTips;

    /* Mapping of cases to branches on the tree; old version is stored before operators are applied */

    private AbstractCase[] branchMap;
    private AbstractCase[] storedBranchMap;

    /* The log likelihood of the subtree from the parent branch of the referred-to node downwards; old version is stored
    before TT operators are applied.*/

    private double[] nodeLogLikelihoods;
    private double[] storedNodeLogLikelihoods;

    /* Whether operations have required the recalculation of the log likelihoods of the subtree from the parent branch
    of the referred-to node upwards.
     */

    private boolean[] nodeRecalculationNeeded;
    private boolean[] storedRecalculationArray;

    /* Matches cases to external nodes */

    private HashMap<AbstractCase, Integer> tipMap;
    private double estimatedLastSampleTime;
    boolean verbose;
    protected TreeTraitProvider.Helper treeTraits = new Helper();

    /* possible paintings for each node */

    private HashMap<Integer, HashSet<AbstractCase>> nodePaintingPossibilities;
    private HashMap<Integer, HashSet<AbstractCase>> storedNodePaintingPossibilities;

    /**
     * The set of cases
     */
    private AbstractOutbreak cases;
    private HashSet<AbstractCase> caseSet;
    private boolean likelihoodKnown = false;
    private double logLikelihood;
    private double storedLogLikelihood;

    // for extended version (not implemented at present)

    private boolean extended;
    private boolean[] switchLocks;
    private boolean[] creepLocks;

    // for reconstructing version version:

    private boolean sampleTTs;
    private double[] rootLikelihoods;
    private double[] storedRootLikelihoods;
    private double[] subLikelihoods;
    private double[] storedSubLikelihoods;
    private boolean currentReconstructionExists;



    // PUBLIC STUFF

    // Name

    public static final String CASE_TO_CASE_TRANSMISSION_LIKELIHOOD = "caseToCaseTransmissionLikelihood";
    public static final String PAINTINGS_KEY = "paintings";

    // Basic constructor.

    public CaseToCaseTransmissionLikelihood(TreeModel virusTree, AbstractOutbreak caseData,
                                            String startingNetworkFileName, boolean extended, boolean sampleTTs)
            throws TaxonList.MissingTaxonException {
        this(CASE_TO_CASE_TRANSMISSION_LIKELIHOOD, virusTree, caseData, startingNetworkFileName, extended, sampleTTs);
    }

    // Legacy constructor

    public CaseToCaseTransmissionLikelihood(String name, TreeModel virusTree, AbstractOutbreak caseData, String
            startingNetworkFileName){
        this(name, virusTree, caseData, startingNetworkFileName, false, false);
    }


    // Constructor for an instance with a non-default name

    public CaseToCaseTransmissionLikelihood(String name, TreeModel virusTree, AbstractOutbreak caseData, String
            startingNetworkFileName, boolean extended, boolean sampleTTs) {

        super(name);
        this.virusTree = virusTree;
        noTips = virusTree.getExternalNodeCount();
        addModel(virusTree);
        Date lastSampleDate = getLatestTaxonDate(virusTree);

        /* We assume samples were taken at the end of the date of the last tip */

        estimatedLastSampleTime = lastSampleDate.getTimeValue();
        cases = caseData;
        caseSet = new HashSet<AbstractCase>(cases.getCases());
        this.extended = extended;
        this.sampleTTs = sampleTTs;
        addModel(cases);
        verbose = false;

        //map cases to tips (needed for both versions)

        branchMap = new AbstractCase[virusTree.getNodeCount()];
        storedBranchMap = new AbstractCase[virusTree.getNodeCount()];
        prepareExternalNodeMap(branchMap);

        tipMap = new HashMap<AbstractCase, Integer>();

        //map the cases to the external nodes
        for(int i=0; i<virusTree.getExternalNodeCount(); i++){
            TreeModel.Node currentExternalNode = (TreeModel.Node)virusTree.getExternalNode(i);
            Taxon currentTaxon = currentExternalNode.taxon;
            for(AbstractCase thisCase : cases.getCases()){
                for(Taxon caseTaxon: thisCase.getAssociatedTaxa()){
                    if(caseTaxon.equals(currentTaxon)){
                        tipMap.put(thisCase, currentExternalNode.getNumber());
                    }
                }
            }
        }

        if(!sampleTTs){
            rootLikelihoods = new double[virusTree.getExternalNodeCount()];
            storedRootLikelihoods = new double[virusTree.getExternalNodeCount()];
            subLikelihoods = new double[virusTree.getInternalNodeCount()
                    *virusTree.getExternalNodeCount()*virusTree.getExternalNodeCount()];
            storedSubLikelihoods = new double[virusTree.getInternalNodeCount()
                    *virusTree.getExternalNodeCount()*virusTree.getExternalNodeCount()];
            currentReconstructionExists = false;
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

        nodePaintingPossibilities = possiblePaintingsMap(virusTree.getRoot(),
                new HashMap<Integer, HashSet<AbstractCase>>());

        treeTraits.addTrait(PAINTINGS_KEY, new TreeTrait.S() {
            public String getTraitName() {
                return PAINTINGS_KEY;
            }

            public Intent getIntent() {
                return Intent.NODE;
            }

            public String getTrait(Tree tree, NodeRef node) {
                return getPaintingForNode(tree, node);
            }


        });

/*      @todo deal with sampleTTs==true and extended==true (probably by prohibiting it)
        not considering the extended case just yet
        if(extended){
            switchLocks = new boolean[virusTree.getInternalNodeCount()];
            creepLocks = new boolean[virusTree.getInternalNodeCount()];
            for(int i=0; i<virusTree.getInternalNodeCount(); i++){
                creepLocks[virusTree.getInternalNode(i).getNumber()] = isCreepLocked(virusTree.getInternalNode(i));
            }
        }*/

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
    // EXTENDED VERSION METHODS - NOT CURRENTLY WORKING
    // ************************************************************************************

    /* check if the given node is tip-linked under the current painting (the tip corresponding to its painting is
    a descendant of it
     */

    public boolean tipLinked(NodeRef node){
        NodeRef tip = virusTree.getNode(tipMap.get(branchMap[node.getNumber()]));
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

    public void recalculatePaintingSets(){
        nodePaintingPossibilities.clear();
        nodePaintingPossibilities = possiblePaintingsMap(virusTree.getRoot(),nodePaintingPossibilities);
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

    // find all valid paintings of the current node

    public HashSet<AbstractCase> possiblePaintings(NodeRef node){
        HashSet<AbstractCase> out = new HashSet<AbstractCase>();
        if(virusTree.isExternal(node)){
            out.add(branchMap[node.getNumber()]);
            return out;
        } else {
            for(int i=0; i<virusTree.getChildCount(node); i++){
                out.addAll(possiblePaintings(virusTree.getChild(node,i)));
            }
            return out;
        }
    }

    // find all paintings of the current node and build an array of this information for this node and its children

    public HashMap<Integer,HashSet<AbstractCase>> possiblePaintingsMap(NodeRef node,
                                                                       HashMap<Integer,HashSet<AbstractCase>> map){
        HashSet<AbstractCase> out = new HashSet<AbstractCase>();
        if(virusTree.isExternal(node)){
            out.add(branchMap[node.getNumber()]);
            map.put(node.getNumber(),out);
        } else {
            for(int i=0; i<virusTree.getChildCount(node); i++){
                out.addAll(possiblePaintingsMap(virusTree.getChild(node,i),map)
                        .get(virusTree.getChild(node,i).getNumber()));
            }
            map.put(node.getNumber(),out);
        }
        return map;
    }


    //Return a set of nodes that are not descendants of (or equal to) the current node and have the same painting as it.
    //The node recalculation is going to need reworking once the time comes to test the extended version,

    public HashSet<Integer> samePaintingDownTree(NodeRef node, boolean flagForRecalc){
        if(!nodeRecalculationNeeded[node.getNumber()] && flagForRecalc){
            extendedflagForRecalculation(virusTree, node, nodeRecalculationNeeded);
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
                    out.addAll(samePaintingUpTree(otherChild, flagForRecalc));
                }
            }
            ancestorNode = virusTree.getParent(ancestorNode);
        }
        return out;
    }

    //Return a set of nodes that are descendants (and not equal to) the current node and have the same painting as it.

    public HashSet<Integer> samePaintingUpTree(NodeRef node, boolean flagForRecalc){
        HashSet<Integer> out = new HashSet<Integer>();
        AbstractCase painting = branchMap[node.getNumber()];
        boolean creepsFurther = false;
        for(int i=0; i<virusTree.getChildCount(node); i++){
            if(branchMap[virusTree.getChild(node,i).getNumber()]==painting){
                creepsFurther = true;
                out.add(virusTree.getChild(node,i).getNumber());
                out.addAll(samePaintingUpTree(virusTree.getChild(node,i), flagForRecalc));
            }
        }
        if(flagForRecalc && !creepsFurther){
            extendedflagForRecalculation(virusTree, node, nodeRecalculationNeeded);
        }
        return out;
    }

    // Return the node numbers of the entire subtree with the same painting as this node (including itself)

    public Integer[] samePainting(NodeRef node, boolean flagForRecalc){
        HashSet<Integer> out = new HashSet<Integer>();
        out.add(node.getNumber());
        out.addAll(samePaintingDownTree(node, flagForRecalc));
        out.addAll(samePaintingUpTree(node, flagForRecalc));
        return out.toArray(new Integer[out.size()]);
    }

    // change flags to indicate that something needs recalculation further up the tree

    private static void extendedflagForRecalculation(TreeModel tree, NodeRef node, boolean[] flags){
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

    public void extendedflagForRecalculation(TreeModel tree, NodeRef node){
        extendedflagForRecalculation(tree, node, nodeRecalculationNeeded);
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
        if(sampleTTs){
//            Arrays.fill(nodeRecalculationNeeded, true);
        } else {
            currentReconstructionExists = false;
        }
        if(model instanceof Tree){
            recalculatePaintingSets();
        }
        likelihoodKnown = false;
    }


    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************


    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if(sampleTTs){
//            Arrays.fill(nodeRecalculationNeeded, true);
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
        if(sampleTTs){
            System.arraycopy(branchMap, 0, storedBranchMap, 0, branchMap.length);
//            storedNodeLogLikelihoods = Arrays.copyOf(nodeLogLikelihoods, nodeLogLikelihoods.length);
//            storedRecalculationArray = Arrays.copyOf(nodeRecalculationNeeded, nodeRecalculationNeeded.length);
        } else {
            System.arraycopy(subLikelihoods, 0, storedSubLikelihoods, 0, subLikelihoods.length);
            System.arraycopy(rootLikelihoods, 0, storedRootLikelihoods, 0, rootLikelihoods.length);
        }
        storedNodePaintingPossibilities = new HashMap<Integer, HashSet<AbstractCase>>(nodePaintingPossibilities);
        storedLogLikelihood = logLikelihood;
    }

    /**
     * Restores the precalculated state.
     */

    protected final void restoreState() {
        if(sampleTTs){
            branchMap = storedBranchMap;
//            nodeLogLikelihoods = storedNodeLogLikelihoods;
//            nodeRecalculationNeeded = storedRecalculationArray;
        } else {
            subLikelihoods = storedSubLikelihoods;
            rootLikelihoods = storedRootLikelihoods;
        }
        nodePaintingPossibilities = storedNodePaintingPossibilities;
        logLikelihood = storedLogLikelihood;
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

    public final TreeModel getVirusTree(){
        return virusTree;
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
        if(sampleTTs){
//            Arrays.fill(nodeRecalculationNeeded, true);
        } else {
            currentReconstructionExists = false;
        }
        if(extended & sampleTTs){
            recalculateLocks();
        }
        likelihoodKnown = false;
        recalculatePaintingSets();
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
    private double calculateLogLikelihood() {
        if(!sampleTTs){
            return totalTreeLogLikelihood(false);
        } else {
            if(!checkPaintingIntegrity(branchMap, true)){
                throw new RuntimeException("Not a painting");
            }
            return calculateLogLikelihood(branchMap);
        }
    }

    private double calculateLogLikelihood(AbstractCase[] map){
        double paintingLogLikelihood = paintingLogLikelihood(map);
        double logNormalisationValue = totalTreeLogLikelihood(false);
        if(paintingLogLikelihood>logNormalisationValue){
            debugOutputTree();
            throw new RuntimeException("Painting likelihood larger than normalisation value. Investigate.");
        }
        return paintingLogLikelihood - logNormalisationValue;
    }

    private double totalTreeLogLikelihood(boolean recordForReconstruction){
        double[] rootLikelihoods = prune(virusTree.getRoot(), recordForReconstruction);
        double normalisationValue = Double.NEGATIVE_INFINITY;
        for (double rootLikelihood : rootLikelihoods) {
            normalisationValue = LogTricks.logSum(normalisationValue, rootLikelihood);
        }
        return normalisationValue;
    }

    private double[] prune(NodeRef node, boolean recordForReconstruction){
        if(!extended){
            if(virusTree.isExternal(node)){
                AbstractCase tipPainting = branchMap[node.getNumber()];
                NodeRef parent = virusTree.getParent(node);
                double[] out = new double[noTips];
                Arrays.fill(out, Double.NEGATIVE_INFINITY);
                HashSet<AbstractCase> possiblePaintings;
                possiblePaintings = nodePaintingPossibilities.get(parent.getNumber());
                for(AbstractCase parentPainting: possiblePaintings){
                    double value = cases.logP(parentPainting, tipPainting, getNodeTime(parent), getNodeTime(node), false);
                    out[cases.getCases().indexOf(parentPainting)]=value;
                }
                return out;
            } else {
                NodeRef parent = virusTree.getParent(node);
                double[][] subtreeLikelihoods = new double[2][noTips];
                double[] out = new double[noTips];
                Arrays.fill(out, Double.NEGATIVE_INFINITY);
                for(int i=0; i<2; i++){
                    System.arraycopy(prune(virusTree.getChild(node,i), recordForReconstruction), 0,
                            subtreeLikelihoods[i], 0, noTips);
                }
                HashSet<AbstractCase> possibleNodePaintings = nodePaintingPossibilities.get(node.getNumber());
                if(parent!=null){
                    // not the root - the likelihood of this subtree given that the parent node has a given painting
                    HashSet<AbstractCase> possibleParentPaintings = nodePaintingPossibilities.get(parent.getNumber());
                    for(AbstractCase parentPainting: possibleParentPaintings){
                        int parentIndex=cases.getCases().indexOf(parentPainting);
                        double sum = Double.NEGATIVE_INFINITY;
                        for(AbstractCase nodePainting: possibleNodePaintings){
                            // is the tip with the parent painting a descendant of this node?
                            boolean paintingForcedByParent = possibleNodePaintings.contains(parentPainting);
                            // If paintingForcedByParent is true, then the likelihood is nonzero only if the node
                            // has the same painting as its parent
                            boolean treeCompatibilityCheck = !paintingForcedByParent || parentPainting==nodePainting;
                            if(treeCompatibilityCheck){
                                double term = 0;
                                int childIndex = cases.getCases().indexOf(nodePainting);
                                for(int i=0; i<2; i++){
                                    term = term + subtreeLikelihoods[i][childIndex];
                                }
                                term = term + cases.logP(parentPainting,nodePainting,getNodeTime(parent),
                                        getNodeTime(node), false);
                                if(recordForReconstruction){
                                    subLikelihoods[(node.getNumber()-noTips)*noTips*noTips
                                            +parentIndex*noTips
                                            +childIndex] = term;
                                }
                                sum = LogTricks.logSum(sum, term);
                            }
                        }
                        out[parentIndex]=sum;
                    }
                    return out;
                } else {
                    // root - these are the likelihoods of the root paintings
                    for(AbstractCase nodePainting: possibleNodePaintings){
                        double term = 0;
                        int index = cases.getCases().indexOf(nodePainting);
                        for(int i=0; i<2; i++){
                            term+=subtreeLikelihoods[i][index];
                        }
                        term+=cases.logProbInfectiousBy(nodePainting,getNodeTime(node), false);
                        if(recordForReconstruction){
                            rootLikelihoods[index]=term;
                        }
                        out[index]=term;
                    }
                    return out;
                }
            }
        } else {
            // extended version
            if(virusTree.isExternal(node)){
                AbstractCase tipPainting = branchMap[node.getNumber()];
                NodeRef parent = virusTree.getParent(node);
                double[] out = new double[noTips];
                Arrays.fill(out, Double.NEGATIVE_INFINITY);
                for(AbstractCase parentPainting: caseSet){
                    double value = cases.logP(parentPainting, tipPainting, getNodeTime(parent), getNodeTime(node),
                            true);
                    out[cases.getCases().indexOf(parentPainting)]=value;
                }
                return out;
            } else {
                NodeRef parent = virusTree.getParent(node);
                double[][] subtreeLikelihoods = new double[2][noTips];
                double[] out = new double[noTips];
                Arrays.fill(out, Double.NEGATIVE_INFINITY);
                for(int i=0; i<2; i++){
                    System.arraycopy(prune(virusTree.getChild(node,i), recordForReconstruction), 0,
                            subtreeLikelihoods[i], 0, noTips);
                }
                HashSet<AbstractCase> nodeDescendantTips = nodePaintingPossibilities.get(node.getNumber());
                if(parent!=null){
                    // not the root - the likelihood of this subtree given that the parent node has a given painting
                    HashSet<AbstractCase> parentDescendantTips = nodePaintingPossibilities.get(parent.getNumber());
                    for(AbstractCase parentPainting: caseSet){
                        int parentIndex=cases.getCases().indexOf(parentPainting);
                        double sum = Double.NEGATIVE_INFINITY;
                        for(AbstractCase nodePainting: caseSet){
                            // is the tip with the parent painting a descendant of this node?
                            boolean paintingForcedByParent = nodeDescendantTips.contains(parentPainting);
                            // are the paintings NOT paintings of a tip descended from the nodes?
                            boolean nodeCreep = !nodeDescendantTips.contains(nodePainting);
                            boolean parentCreep = !parentDescendantTips.contains(parentPainting);
                            // valid combinations:
                            // 1) no creep in either case, paintingForcedByParent = false
                            // 2) no creep in either case, paintingForcedByParent = true and both nodes have same
                            // painting
                            // 3) parent creep but no node creep
                            // 4) node creep, parent and child have same painting
                            boolean treeCompatibilityCheck = (!nodeCreep && !parentCreep && (!paintingForcedByParent
                                    || parentPainting == nodePainting)) || (parentCreep && !nodeCreep) ||
                                    (nodeCreep && (parentPainting == nodePainting));
                            if(treeCompatibilityCheck){
                                double term = 0;
                                int childIndex = cases.getCases().indexOf(nodePainting);
                                for(int i=0; i<2; i++){
                                    term = term + subtreeLikelihoods[i][childIndex];
                                }
                                term = term + cases.logP(parentPainting,nodePainting,getNodeTime(parent),
                                        getNodeTime(node), true);
                                if(recordForReconstruction){
                                    subLikelihoods[(node.getNumber()-noTips)*noTips*noTips
                                            +parentIndex*noTips
                                            +childIndex] = term;
                                }
                                sum = LogTricks.logSum(sum, term);
                            }
                        }
                        out[parentIndex]=sum;
                    }
                    return out;
                } else {
                    // root - these are the likelihoods of the root paintings
                    for(AbstractCase nodePainting: caseSet){
                        double term = 0;
                        int index = cases.getCases().indexOf(nodePainting);
                        for(int i=0; i<2; i++){
                            term+=subtreeLikelihoods[i][index];
                        }
                        term+=cases.logProbInfectiousBy(nodePainting,getNodeTime(node), true);
                        if(recordForReconstruction){
                            rootLikelihoods[index]=term;
                        }
                        out[index]=term;
                    }
                    return out;

                }
            }
        }
    }

    private double paintingLogLikelihood(AbstractCase[] map){
        return prune(map, virusTree.getRoot());
    }

    private double prune(AbstractCase[] map, NodeRef node){
        if(virusTree.isExternal(node)){
            NodeRef parent = virusTree.getParent(node);
            AbstractCase tipPainting = map[node.getNumber()];
            AbstractCase parentPainting = map[parent.getNumber()];
            return cases.logP(parentPainting, tipPainting, getNodeTime(parent), getNodeTime(node), extended);
        } else {
            NodeRef parent = virusTree.getParent(node);
            double[] subtreeLikelihoods = new double[2];
            for(int i=0; i<2; i++){
                subtreeLikelihoods[i] = prune(map, virusTree.getChild(node,i));
            }
            AbstractCase nodePainting = map[node.getNumber()];
            if(parent!=null){
                // not the root - the likelihood of this subtree given that the parent node has a given painting
                AbstractCase parentPainting = map[parent.getNumber()];
                double term = 0;
                for(int i=0; i<2; i++){
                    term+=subtreeLikelihoods[i];
                }
                term += cases.logP(parentPainting, nodePainting, getNodeTime(parent), getNodeTime(node), extended);
                return term;

            } else {
                // root - these are the likelihoods of the root paintings
                double term = 0;
                for(int i=0; i<2; i++){
                    term+=subtreeLikelihoods[i];
                }
                term+=cases.logProbInfectiousBy(nodePainting, getNodeTime(node), extended);
                return term;
            }
        }
    }



    /* Return the integer day on which the given node occurred */

    public int getNodeDay(NodeRef node){
        Double nodeHeight = getHeight(node);
        Double estimatedNodeTime = estimatedLastSampleTime-nodeHeight;
        Date nodeDate = new Date(estimatedNodeTime, Units.Type.DAYS, false);
        /*Since day t begins at time t-1, to find the day of a given decimal time value, take the ceiling.*/
        return (int)Math.ceil(nodeDate.getTimeValue());
    }

    /* Return the double time at which the given node occurred */

    public double getNodeTime(NodeRef node){
        double nodeHeight = getHeight(node);
        return estimatedLastSampleTime-nodeHeight;
    }

    /**
     * Given a node, calculates the log likelihood of its parent branch and then goes down the tree and calculates the
     * log likelihoods of lower branches.
     */

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
        NodeRef tip = virusTree.getNode(tipMap.get(thisCase));
        return getInfector(tip, branchMap);
    }

    public AbstractCase getInfector(NodeRef node, AbstractCase[] branchMap){

        try{
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
        } catch (NullPointerException e){
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
    }

    public boolean checkPaintingIntegrity(AbstractCase[] map, boolean verbose){
        boolean out=true;
        for(int i=0; i<virusTree.getInternalNodeCount(); i++){
            NodeRef node = virusTree.getInternalNode(i);
            NodeRef firstChild = virusTree.getChild(node,0);
            NodeRef secondChild = virusTree.getChild(node,1);
            NodeRef parent = virusTree.getParent(node);
            if(map[node.getNumber()]!=map[firstChild.getNumber()] &&
                    map[node.getNumber()]!=map[secondChild.getNumber()] &&
                    (!extended || map[node.getNumber()]!=map[parent.getNumber()])){
                out = false;
                if(!verbose){
                    break;
                } else {
                    System.out.println("Node "+node.getNumber()+" failed painting integrity check:");
                    System.out.println("Node painting: "+map[node.getNumber()].getName());
                    System.out.println("Parent painting: "+map[parent.getNumber()].getName());
                    System.out.println("Child 1 painting: "+map[firstChild.getNumber()].getName());
                    System.out.println("Child 2 painting: "+map[secondChild.getNumber()].getName());
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
            System.arraycopy(map, 0, branchMap, 0, map.length);
            //Warning - if the BadPaintingException in randomlyPaintNode might be caused by a bug rather than both
            //likelihoods rounding to zero, you want to stop catching this to investigate.
            paintRandomNetwork(map,true);
            if(paintingLogLikelihood(map)!=Double.NEGATIVE_INFINITY){
                gotOne = true;
                System.out.println("found.");
            }
            tries++;
            if(tries==101){
                System.out.println("giving " +
                        "up.");
                throw new RuntimeException("Failed to find a starting transmission network with nonzero likelihood");
            }
        }
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
                    branchLogLs[i]= cases.logP(choices[i], choices[1-i], getNodeTime(node),
                            getNodeTime(virusTree.getChild(node,1-i)), extended);
                }
                if(branchLogLs[0]==Double.NEGATIVE_INFINITY && branchLogLs[1]==Double.NEGATIVE_INFINITY){

                    for(int i=0; i<2; i++){
                        branchLogLs[i]= cases.logP(choices[i], choices[1-i], getNodeTime(node),
                                getNodeTime(virusTree.getChild(node,1-i)), extended);
                    }
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

    /* Sample a random transmission tree from the likelihoods for the paintings of each triplet of nodes */

    private AbstractCase[] sampleTransmissionTree(){

        Arrays.fill(rootLikelihoods, Double.NEGATIVE_INFINITY);
        Arrays.fill(subLikelihoods, Double.NEGATIVE_INFINITY);

        totalTreeLogLikelihood(true);

        AbstractCase[] samplePainting = new AbstractCase[virusTree.getNodeCount()];
        samplePainting = prepareExternalNodeMap(samplePainting);
        NodeRef root = virusTree.getRoot();
        try{
            int choice = MathUtils.randomChoicePDF(convertAndNormalise(rootLikelihoods));
            samplePainting[root.getNumber()]=cases.getCase(choice);
        } catch (Error e){
            debugOutputTree();
            throw new RuntimeException("Problem with painting reconstruction");
        }

        for(int i=0; i<2; i++){
            fillUp(virusTree.getChild(root,i),samplePainting);
        }

        if(!checkPaintingIntegrity(samplePainting, true)){
            throw new RuntimeException("The reconstruction method isn't outputting proper paintings");
        }
        currentReconstructionExists = true;
        return samplePainting;

    }

    //randomChoicePDF is failing because sometimes all the probabilities round to zero on the real scale. This takes
    // an array of probabilities on the log scale, normalises, and returns the normalised probabilities in real space

    private double[] convertAndNormalise(double[] array){
        double sum = Double.NEGATIVE_INFINITY;
        double[] out = new double[array.length];
        for(double value: array){
            sum = LogTricks.logSum(sum,value);
        }
        for(int i=0; i<array.length; i++){
            out[i] = Math.exp(array[i] - sum);
        }
        return out;
    }


    private void fillUp(NodeRef node, AbstractCase[] painting){

        if(!virusTree.isExternal(node)){
            int dim = virusTree.getExternalNodeCount();
            AbstractCase parentPainting = painting[virusTree.getParent(node).getNumber()];
            int parentPaintingNumber = cases.getCases().indexOf(painting[virusTree.getParent(node).getNumber()]);
            double[] childLikelihoods = Arrays.copyOfRange(subLikelihoods,
                    (node.getNumber()-dim)*dim*dim + parentPaintingNumber*dim,
                    (node.getNumber()-dim)*dim*dim + (parentPaintingNumber+1)*dim);
            try{
                int choice = MathUtils.randomChoicePDF(convertAndNormalise(childLikelihoods));
                painting[node.getNumber()]=cases.getCase(choice);
            } catch(Error e){
                debugOutputTree();
                throw new RuntimeException("Problem with painting reconstruction");
            }
            for(int i=0; i<2; i++){
                NodeRef child = virusTree.getChild(node,i);
                if(!virusTree.isExternal(child)){
                    fillUp(child, painting);
                }
            }
        }
    }

    private void debugOutputTree(){
        try{
            FlexibleTree treeCopy = new FlexibleTree(virusTree);
            for(int j=0; j<treeCopy.getNodeCount(); j++){
                FlexibleNode node = (FlexibleNode)treeCopy.getNode(j);
                node.setAttribute("Number", node.getNumber());
            }
            NexusExporter testTreesOut = new NexusExporter(new PrintStream("testTrees.nex"));
            testTreesOut.exportTree(treeCopy);
        } catch (IOException ignored) {System.out.println("IOException");}
        if(!sampleTTs){
            Arrays.fill(subLikelihoods, Double.NEGATIVE_INFINITY);
            Arrays.fill(rootLikelihoods, Double.NEGATIVE_INFINITY);
            totalTreeLogLikelihood(true);
        } else {
            totalTreeLogLikelihood(false);
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
            final boolean sampleTTs = xo.getBooleanAttribute("sampleTTs");

            try {
                likelihood = new CaseToCaseTransmissionLikelihood(virusTree, caseSet, startingNetworkFileName,
                        extended, sampleTTs);
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
                AttributeRule.newBooleanRule("extended"),  // don't do this!!
                AttributeRule.newBooleanRule("sampleTTs"),
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

        if(!sampleTTs){

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

        if(!sampleTTs && !currentReconstructionExists){
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




