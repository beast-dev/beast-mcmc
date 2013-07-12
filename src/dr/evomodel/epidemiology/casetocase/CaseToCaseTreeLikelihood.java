package dr.evomodel.epidemiology.casetocase;
import dr.app.tools.NexusExporter;
import dr.evolution.tree.*;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.AbstractTreeLikelihood;
import dr.evomodel.treelikelihood.GeneralLikelihoodCore;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
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
import java.util.logging.Logger;

/**
 * Handles manipulation of the tree partition, and likelihood of the infection times.
 *
 * This is not an actual implementation of AbstractTreeLikelihood since it doesn't
 *
 * @author Matthew Hall
 * @author Andrew Rambaut
 * @version $Id: $
 */
public class CaseToCaseTreeLikelihood extends AbstractTreeLikelihood implements Loggable, Citable,
        TreeTraitProvider {

    /* The phylogenetic tree. */

    private int noTips;

    /* Mapping of cases to branches on the tree; old version is stored before operators are applied */

    private AbstractCase[] branchMap;
    private AbstractCase[] storedBranchMap;

/* The log likelihood of the subtree from the parent branch of the referred-to node downwards; old version is stored
    before TT operators are applied. @todo work out if these are still necessary*/


    private double[] nodeLogLikelihoods;
    private double[] storedNodeLogLikelihoods;

/* Whether operations have required the recalculation of the log likelihoods of the subtree from the parent branch
    of the referred-to node upwards. @todo work out if these are still necessary
     */

    private boolean[] nodeRecalculationNeeded;
    private boolean[] storedRecalculationArray;


    /* Matches cases to external nodes */

    private HashMap<AbstractCase, Integer> tipMap;
    private double estimatedLastSampleTime;
    boolean verbose;
    protected TreeTraitProvider.Helper treeTraits = new Helper();

    /* the partitions that all descendant tips of a current node belong to */

    private HashMap<Integer, HashSet<AbstractCase>> descendantTipPartitions;
    private HashMap<Integer, HashSet<AbstractCase>> storedDescendantTipPartitions;

    /**
     * The set of cases
     */
    private AbstractOutbreak cases;
    private HashSet<AbstractCase> caseSet;
    private boolean likelihoodKnown = false;
    private double logLikelihood;
    private double storedLogLikelihood;
    private CaseToCaseLikelihoodCore core;

    // for individual partition calculations

    private double[] branchLogProbs;
    private double totalLogProb;

    // for normalisation

    private double[] normRootPartials;
    private double[] normProbs;
    private double normTotalProb;


    // where along the relevant branches the infections happen

    private Parameter infectionTimes;

    //ugly, ugly hack. Maybe not necessary. For heaven's sake don't estimate it, at any rate. Set it to be the
    //earliest reasonable time of first infection

    private Parameter earliestPossibleFirstInfection;

    // for extended version

    private boolean extended;

    // PUBLIC STUFF

    // Name

    public static final String CASE_TO_CASE_TREE_LIKELIHOOD = "caseToCaseTreeLikelihood";
    public static final String PARTITIONS_KEY = "partition";

    // Basic constructor.

    public CaseToCaseTreeLikelihood(TreeModel virusTree, AbstractOutbreak caseData,
                                    String startingNetworkFileName, Parameter infectionTimes,
                                    Parameter earliestPossibleFirstInfection, boolean extended)
            throws TaxonList.MissingTaxonException {
        this(CASE_TO_CASE_TREE_LIKELIHOOD, virusTree, caseData, startingNetworkFileName, infectionTimes,
                earliestPossibleFirstInfection, extended);
    }

    // Constructor for an instance with a non-default name

    public CaseToCaseTreeLikelihood(String name, TreeModel virusTree, AbstractOutbreak caseData, String
            startingNetworkFileName, Parameter infectionTimes, Parameter earliestPossibleFirstInfection,
                                    boolean extended) {
        super(name, caseData, virusTree);
        if(stateCount!=nodeCount){
            throw new RuntimeException("There are duplicate tip cases.");
        }

        noTips = virusTree.getExternalNodeCount();

        Date lastSampleDate = getLatestTaxonDate(virusTree);
        estimatedLastSampleTime = lastSampleDate.getTimeValue();

        cases = caseData;
        caseSet = new HashSet<AbstractCase>(cases.getCases());
        this.extended = extended;
        addModel(cases);
        verbose = false;

        //map cases to tips

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

        this.infectionTimes = infectionTimes;
        this.earliestPossibleFirstInfection = earliestPossibleFirstInfection;


        // @todo is the following still necessary?

        nodeRecalculationNeeded = new boolean[virusTree.getNodeCount()];
        Arrays.fill(nodeRecalculationNeeded, true);
        nodeLogLikelihoods = new double[virusTree.getNodeCount()];

        core = new CaseToCaseLikelihoodCore();
        core.initialize(noTips, 1, 1, false);

        branchLogProbs = new double[nodeCount];
        totalLogProb = 0;

        normRootPartials = new double[stateCount];
        normProbs = new double[stateCount*stateCount];

        //paint the starting network onto the tree

        if(startingNetworkFileName==null){
            partitionAccordingToRandomTT();
        } else {
            partitionAccordingToSpecificTT(startingNetworkFileName);
        }

        descendantTipPartitions = new HashMap<Integer, HashSet<AbstractCase>>();
        descendantTipPartitions(virusTree.getRoot(), descendantTipPartitions);

        treeTraits.addTrait(PARTITIONS_KEY, new TreeTrait.S() {
            public String getTraitName() {
                return PARTITIONS_KEY;
            }

            public Intent getIntent() {
                return Intent.NODE;
            }

            public String getTrait(Tree tree, NodeRef node) {
                return getNodePartition(tree, node);
            }


        });

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
        NodeRef[] children = new NodeRef[treeModel.getChildCount(node)];
        for(int i=0; i< treeModel.getChildCount(node); i++){
            children[i] = treeModel.getChild(node,i);
        }
        return children;
    }


    // ************************************************************************************
    // EXTENDED VERSION METHODS
    // ************************************************************************************

    /* check if the given node is tip-linked under the current painting (the tip corresponding to its painting is
    a descendant of it
     */

    public boolean tipLinked(NodeRef node){
        NodeRef tip = treeModel.getNode(tipMap.get(branchMap[node.getNumber()]));
        if(tip==node){
            return true;
        }
        NodeRef parent = tip;
        while(parent!= treeModel.getRoot()){
            parent = treeModel.getParent(parent);
            if(parent==node){
                return true;
            }
        }
        return false;
    }

    public boolean isExtended(){
        return extended;
    }

    public void recalculateDescTipPartitions(){
        descendantTipPartitions.clear();
        descendantTipPartitions(treeModel.getRoot(), descendantTipPartitions);
    }

    //Counts the children of the current node which have the same painting as itself under the current map.
    //This will always be 1 if extended==false.

    public int countChildrenInSamePartition(NodeRef node){
        if(treeModel.isExternal(node)){
            return -1;
        } else {
            int count = 0;
            AbstractCase parentCase = branchMap[node.getNumber()];
            for(int i=0; i< treeModel.getChildCount(node); i++){
                if(branchMap[treeModel.getChild(node,i).getNumber()]==parentCase){
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

    // find all paritions of the descendent tips of the current node. If map is specified then it makes a map of node
    // number to possible partitions; map can be null.

    public HashSet<AbstractCase> descendantTipPartitions(NodeRef node, HashMap<Integer, HashSet<AbstractCase>> map){
        HashSet<AbstractCase> out = new HashSet<AbstractCase>();
        if(treeModel.isExternal(node)){
            out.add(branchMap[node.getNumber()]);
            if(map!=null){
                map.put(node.getNumber(), out);
            }
            return out;
        } else {
            for(int i=0; i< treeModel.getChildCount(node); i++){
                out.addAll(descendantTipPartitions(treeModel.getChild(node, i), map));
            }
            if(map!=null){
                map.put(node.getNumber(), out);
            }
            return out;
        }
    }

    //Return a set of nodes that are not descendants of (or equal to) the current node and are in the same partition as
    // it. If flagForRecalc is true, then this also sets the flags for likelihood recalculation for all these nodes
    // to true

    public HashSet<Integer> samePartitionDownTree(NodeRef node, boolean flagForRecalc){
        if(!nodeRecalculationNeeded[node.getNumber()] && flagForRecalc){
            flagForDescendantRecalculation(treeModel, node, nodeRecalculationNeeded);
        }
        HashSet<Integer> out = new HashSet<Integer>();
        AbstractCase painting = branchMap[node.getNumber()];
        NodeRef ancestorNode = treeModel.getParent(node);
        while(branchMap[ancestorNode.getNumber()]==painting){
            out.add(ancestorNode.getNumber());
            if(extended){
                if(countChildrenInSamePartition(ancestorNode)==2){
                    NodeRef otherChild = sibling(treeModel, ancestorNode);
                    out.add(otherChild.getNumber());
                    out.addAll(samePartitionUpTree(otherChild, flagForRecalc));
                }
            }
            ancestorNode = treeModel.getParent(ancestorNode);
        }
        return out;
    }

    //Return a set of nodes that are descendants (and not equal to) the current node and are in the same parition as it.

    public HashSet<Integer> samePartitionUpTree(NodeRef node, boolean flagForRecalc){
        HashSet<Integer> out = new HashSet<Integer>();
        AbstractCase painting = branchMap[node.getNumber()];
        boolean creepsFurther = false;
        for(int i=0; i< treeModel.getChildCount(node); i++){
            if(branchMap[treeModel.getChild(node,i).getNumber()]==painting){
                creepsFurther = true;
                out.add(treeModel.getChild(node,i).getNumber());
                out.addAll(samePartitionUpTree(treeModel.getChild(node, i), flagForRecalc));
            }
        }
        if(flagForRecalc && !creepsFurther){
            flagForDescendantRecalculation(treeModel, node, nodeRecalculationNeeded);
        }
        return out;
    }

    // Return the node numbers of the entire subtree in the same partition as this node (including itself)

    public Integer[] samePartition(NodeRef node, boolean flagForRecalc){
        HashSet<Integer> out = new HashSet<Integer>();
        out.add(node.getNumber());
        out.addAll(samePartitionDownTree(node, flagForRecalc));
        out.addAll(samePartitionUpTree(node, flagForRecalc));
        return out.toArray(new Integer[out.size()]);
    }

    // change flags to indicate that something needs recalculation further up the tree

    private static void flagForDescendantRecalculation(TreeModel tree, NodeRef node, boolean[] flags){
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

    public void flagForDescendantRecalculation(TreeModel tree, NodeRef node){
        flagForDescendantRecalculation(tree, node, nodeRecalculationNeeded);
    }

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************


    protected final void handleModelChangedEvent(Model model, Object object, int index) {
        // @todo check if this is doing what it's supposed to.
        Arrays.fill(nodeRecalculationNeeded, true);
        if(model instanceof Tree && !extended){
            recalculateDescTipPartitions();
        }
        likelihoodKnown = false;
    }


    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************


    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // @todo check if this is doing what it's supposed to
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
        System.arraycopy(branchMap, 0, storedBranchMap, 0, branchMap.length);
        storedNodeLogLikelihoods = Arrays.copyOf(nodeLogLikelihoods, nodeLogLikelihoods.length);
        storedRecalculationArray = Arrays.copyOf(nodeRecalculationNeeded, nodeRecalculationNeeded.length);
        storedDescendantTipPartitions = new HashMap<Integer, HashSet<AbstractCase>>(descendantTipPartitions);
        storedLogLikelihood = logLikelihood;
    }

    /**
     * Restores the precalculated state.
     */

    protected final void restoreState() {
        branchMap = storedBranchMap;
        nodeLogLikelihoods = storedNodeLogLikelihoods;
        nodeRecalculationNeeded = storedRecalculationArray;
        descendantTipPartitions = storedDescendantTipPartitions;
        logLikelihood = storedLogLikelihood;
    }

    protected final void acceptState() {
    } // nothing to do

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    public final AbstractCase[] getBranchMap(){
        return branchMap;
    }

    public final TreeModel getTreeModel(){
        return treeModel;
    }

    public final void setBranchMap(AbstractCase[] map){
        branchMap = map;
    }

    public final TreeModel getTree(){
        return treeModel;
    }

    public final void makeDirty() {
        //@todo I don't think you need to force recalculation of the whole tree here, but check
//            Arrays.fill(nodeRecalculationNeeded, true);

        likelihoodKnown = false;
        if(!extended){
            recalculateDescTipPartitions();
        }
    }

    public void makeDirty(boolean cleanTree){
        if(!cleanTree){
            makeDirty();
        } else {
            Arrays.fill(nodeRecalculationNeeded, true);
            likelihoodKnown = false;
            if(!extended){
                recalculateDescTipPartitions();
            }
        }
    }

    /**
     * Calculates the log likelihood of this set of node labels given the tree.
     */
    protected double calculateLogLikelihood() {
        //@todo a consistent system of checking painting integrity
        if(!checkPartitions(branchMap, true)){
            throw new RuntimeException("Partition rules are violated");
        }


        final NodeRef root = treeModel.getRoot();

        // unnormalised probability
        traverse(treeModel, root, branchMap,
                infectionTimes.getParameterValues(), earliestPossibleFirstInfection.getParameterValue(0));
        // normalisation value
        traverse(treeModel, root);

        double logL = (normTotalProb) * patternWeights[0];

        if (logL == Double.NEGATIVE_INFINITY) {
            Logger.getLogger("dr.evomodel").info("TreeLikelihood, " + this.getId() + ", turning on partial likelihood scaling to avoid precision loss");

            // We probably had an underflow... turn on scaling
            core.setUseScaling(true);

            // and try again...
            updateAllNodes();
            updateAllPatterns();
            traverse(treeModel, root);

            logL = (normTotalProb) * patternWeights[0];

        }

        for (int i = 0; i < nodeCount; i++) {
            updateNode[i] = false;
        }

        return logL;

    }

    // normalisation value

    private boolean traverse(Tree tree, NodeRef node){
        boolean update = false;

        int nodeNum = node.getNumber();

        NodeRef parent = tree.getParent(node);

        if (updateNode[nodeNum] && parent!=null) {
            HashSet<AbstractCase> nodeDescTips = descendantTipPartitions.get(node.getNumber());
            HashSet<AbstractCase> parentDescTips = descendantTipPartitions.get(parent.getNumber());
            for(int i=0; i<stateCount; i++){
                AbstractCase origin = cases.getCase(i);
                for(int j=0; j<stateCount; j++){
                    AbstractCase destination = cases.getCase(j);
                    // this can probably be sped up because a lot of entries are going to be zero.

                    // is the tip in parent's partition a descendant of this node? If so, the node must be in
                    // this partition also
                    boolean paintingForcedByParent = nodeDescTips.contains(origin);
                    boolean treeCompatibilityCheck;
                    if(!extended){
                        // valid combinations:
                        // 1) paintingForcedByParent = false
                        // 2) paintingForcedByParent = true and both nodes in same partition
                        treeCompatibilityCheck = nodeDescTips.contains(destination) && parentDescTips.contains(origin)
                                && (!paintingForcedByParent || origin==destination);
                    } else {
                        boolean nodeCreep = !nodeDescTips.contains(destination);
                        boolean parentCreep = !parentDescTips.contains(origin);
                        // valid combinations:
                        // 1) no creep in either case, paintingForcedByParent = false
                        // 2) no creep in either case, paintingForcedByParent = true and both nodes in same
                        // partition
                        // 3) parent creep but no node creep
                        // 4) node creep, parent and child in same partition
                        treeCompatibilityCheck =  nodeDescTips.contains(destination) &&
                                parentDescTips.contains(origin) &&
                                ((!nodeCreep && !parentCreep && (!paintingForcedByParent || origin == destination))
                                        || (parentCreep && !nodeCreep) || (nodeCreep && (origin == destination)));
                    }
                    if(!treeCompatibilityCheck){
                        normProbs[stateCount*i + j]=0;
                    } else if(origin==destination) {
                        normProbs[stateCount*i + j]=1;
                    } else {
                        normProbs[stateCount*i + j]=cases.probXInfectedByYBetweenTandU(destination, origin,
                                getNodeTime(parent), getNodeTime(node));
                    }
                }
            }
            // bothered by the possibility of rounding to zero in this
            core.setNodeMatrix(nodeNum, 1, normProbs);
            update = true;
        }
        if (!tree.isExternal(node)) {

            // Traverse down the two child nodes
            NodeRef child1 = tree.getChild(node, 0);
            final boolean update1 = traverse(tree, child1);

            NodeRef child2 = tree.getChild(node, 1);
            final boolean update2 = traverse(tree, child2);

            if (update1 || update2) {

                final int childNum1 = child1.getNumber();
                final int childNum2 = child2.getNumber();

                core.setNodePartialsForUpdate(nodeNum);

                core.calculatePartials(childNum1, childNum2, nodeNum, null);

                if (parent == null) {
                    // No parent this is the root of the tree -
                    // calculate the pattern likelihoods

                    //

                    double[] partials = getRootPartials();

                    core.calculateLogLikelihoods(partials, normTotalProb);
                }
                update = true;
            }
        }
        return update;
    }

    public final double[] getRootPartials() {
        if (normRootPartials == null) {
            normRootPartials = new double[stateCount];
        }

        int nodeNum = treeModel.getRoot().getNumber();

        core.getPartials(nodeNum, normRootPartials);

        // now need to deal with the root branch

        for(int i=0; i<noTips; i++){
            normRootPartials[i] *= cases.probYInfectedBetweenTandU
                    (cases.getCase(i),heightToTime(earliestPossibleFirstInfection.getParameterValue(0)),
                            getNodeTime(treeModel.getRoot()));
        }


        return normRootPartials;
    }






    // unnormalised probability of this branchMap and corresponding infTimes

    private boolean traverse(Tree tree, NodeRef node, AbstractCase[] branchMap, double[] infTimes,
                             double earliestFirstInfection){

        // @todo probably better with two updateNodes - one for the normalisation procedure and one for this
        // @todo because changes to infTimes will matter in the latter case but not the former

        boolean update = false;

        int nodeNum = node.getNumber();

        NodeRef parent = tree.getParent(node);

        if (updateNode[nodeNum]) {
            if(parent!=null){
                AbstractCase nodeCase = branchMap[node.getNumber()];
                AbstractCase parentCase = branchMap[parent.getNumber()];
                if(nodeCase!=parentCase){
                    final double branchLength = tree.getNodeHeight(parent) - tree.getNodeHeight(node);
                    final double infectionTime = heightToTime(tree.getNodeHeight(node)
                            - branchLength*infTimes[cases.getCaseIndex(nodeCase)]);
                    branchLogProbs[node.getNumber()]
                            = cases.logProbXInfectedByYAtTimeT(nodeCase, parentCase, infectionTime);
                } else {
                    branchLogProbs[node.getNumber()] = 0;
                }
            } else {
                AbstractCase nodeCase = branchMap[node.getNumber()];
                final double branchLength = earliestFirstInfection - tree.getNodeHeight(node);
                final double infectionTime = heightToTime(tree.getNodeHeight(node)
                        - branchLength*infTimes[cases.getCaseIndex(nodeCase)]);
                branchLogProbs[node.getNumber()]
                        = cases.logProbYInfectedAtTimeT(nodeCase, infectionTime);
            }
        }
        if (!tree.isExternal(node)) {
            NodeRef child1 = tree.getChild(node, 0);
            final boolean update1 = traverse(tree, child1, branchMap, infTimes, earliestFirstInfection);
            NodeRef child2 = tree.getChild(node, 1);
            final boolean update2 = traverse(tree, child2, branchMap, infTimes, earliestFirstInfection);
            if (update1 || update2) {
                if (parent == null) {
                    totalLogProb = 0;
                    for(int i=0; i<nodeCount; i++){
                        totalLogProb += branchLogProbs[i];
                    }
                    if(totalLogProb == Double.NEGATIVE_INFINITY){
                        throw new RuntimeException("Total log probability of partition is zero");
                    }
                }
                update = true;
            }
        }
        return(update);
    }


    // start of deprecated stuff

/*
    private double calculateLogLikelihood(AbstractCase[] map){
        double partitionLogLikelihood = partitionLogLikelihood(map);
        double logNormalisationValue = totalTreeLogLikelihood();
        if(partitionLogLikelihood>logNormalisationValue){
            debugOutputTree();
            throw new RuntimeException("Painting likelihood larger than normalisation value. Investigate.");
        }
        return partitionLogLikelihood - logNormalisationValue;
    }

    private double totalTreeLogLikelihood(){
        double[] rootLikelihoods = prune(treeModel.getRoot());
        double normalisationValue = Double.NEGATIVE_INFINITY;
        for (double rootLikelihood : rootLikelihoods) {
            normalisationValue = LogTricks.logSum(normalisationValue, rootLikelihood);
        }
        return normalisationValue;
    }

    private double[] prune(NodeRef node){
        if(!extended){
            if(treeModel.isExternal(node)){
                AbstractCase tipPainting = branchMap[node.getNumber()];
                NodeRef parent = treeModel.getParent(node);
                double[] out = new double[noTips];
                Arrays.fill(out, Double.NEGATIVE_INFINITY);
                HashSet<AbstractCase> possiblePaintings;
                possiblePaintings = descendantTipPartitions.get(parent.getNumber());
                for(AbstractCase parentPainting: possiblePaintings){
                    double value = cases.logProbXInfectedByYBetweenTandU(parentPainting, tipPainting, getNodeTime(parent), getNodeTime(node), false);
                    out[cases.getCases().indexOf(parentPainting)]=value;
                }
                return out;
            } else {
                NodeRef parent = treeModel.getParent(node);
                double[][] subtreeLikelihoods = new double[2][noTips];
                double[] out = new double[noTips];
                Arrays.fill(out, Double.NEGATIVE_INFINITY);
                for(int i=0; i<2; i++){
                    System.arraycopy(prune(treeModel.getChild(node,i)), 0,
                            subtreeLikelihoods[i], 0, noTips);
                }
                HashSet<AbstractCase> possibleNodePaintings = descendantTipPartitions.get(node.getNumber());
                if(parent!=null){
                    // not the root - the likelihood of this subtree given that the parent node has a given painting
                    HashSet<AbstractCase> possibleParentPaintings = descendantTipPartitions.get(parent.getNumber());
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
                        out[index]=term;
                    }
                    return out;
                }
            }
        } else {
            // extended version
            if(treeModel.isExternal(node)){
                AbstractCase tipPainting = branchMap[node.getNumber()];
                NodeRef parent = treeModel.getParent(node);
                double[] out = new double[noTips];
                Arrays.fill(out, Double.NEGATIVE_INFINITY);
                for(AbstractCase parentPainting: caseSet){
                    double value = cases.logP(parentPainting, tipPainting, getNodeTime(parent), getNodeTime(node),
                            true);
                    out[cases.getCases().indexOf(parentPainting)]=value;
                }
                return out;
            } else {
                NodeRef parent = treeModel.getParent(node);
                double[][] subtreeLikelihoods = new double[2][noTips];
                double[] out = new double[noTips];
                Arrays.fill(out, Double.NEGATIVE_INFINITY);
                for(int i=0; i<2; i++){
                    System.arraycopy(prune(treeModel.getChild(node,i)), 0,
                            subtreeLikelihoods[i], 0, noTips);
                }
                HashSet<AbstractCase> nodeDescendantTips = descendantTipPartitions.get(node.getNumber());
                if(parent!=null){
                    // not the root - the likelihood of this subtree given that the parent node has a given painting
                    HashSet<AbstractCase> parentDescendantTips = descendantTipPartitions.get(parent.getNumber());
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
                        out[index]=term;
                    }
                    return out;

                }
            }
        }
    }

    private double partitionLogLikelihood(AbstractCase[] map){
        return prune(map, treeModel.getRoot());
    }

    private double prune(AbstractCase[] map, NodeRef node){
        if(treeModel.isExternal(node)){
            NodeRef parent = treeModel.getParent(node);
            AbstractCase tipPainting = map[node.getNumber()];
            AbstractCase parentPainting = map[parent.getNumber()];
            return cases.logP(parentPainting, tipPainting, getNodeTime(parent), getNodeTime(node), extended);
        } else {
            NodeRef parent = treeModel.getParent(node);
            double[] subtreeLikelihoods = new double[2];
            for(int i=0; i<2; i++){
                subtreeLikelihoods[i] = prune(map, treeModel.getChild(node,i));
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
*/

    /* Return the double time at which the given node occurred */

    public double getNodeTime(NodeRef node){
        double nodeHeight = getHeight(node);
        return estimatedLastSampleTime-nodeHeight;
    }

    public double heightToTime(double height){
        return estimatedLastSampleTime-height;
    }

    /**
     * Given a node, calculates the log likelihood of its parent branch and then goes down the tree and calculates the
     * log likelihoods of lower branches.
     */

    public boolean[] getRecalculationArray(){
        return nodeRecalculationNeeded;
    }

    private double getHeight(NodeRef node){
        return treeModel.getNodeHeight(node);
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
        NodeRef tip = treeModel.getNode(tipMap.get(thisCase));
        return getInfector(tip, branchMap);
    }

    public AbstractCase getInfector(NodeRef node, AbstractCase[] branchMap){
        if(treeModel.isRoot(node) || node.getNumber()== treeModel.getRoot().getNumber()){
            return null;
        } else {
            AbstractCase nodeCase = branchMap[node.getNumber()];
            if(branchMap[treeModel.getParent(node).getNumber()]!=nodeCase){
                return branchMap[treeModel.getParent(node).getNumber()];
            } else {
                return getInfector(treeModel.getParent(node), branchMap);
            }
        }
    }

    public boolean checkPartitions(){
        return checkPartitions(branchMap, true);
    }

    private boolean checkPartitions(AbstractCase[] map, boolean verbose){
        boolean out=true;
        for(int i=0; i< treeModel.getInternalNodeCount(); i++){
            NodeRef node = treeModel.getInternalNode(i);
            NodeRef firstChild = treeModel.getChild(node,0);
            NodeRef secondChild = treeModel.getChild(node,1);
            NodeRef parent = treeModel.getParent(node);
            if(map[node.getNumber()]!=map[firstChild.getNumber()] &&
                    map[node.getNumber()]!=map[secondChild.getNumber()] &&
                    (!extended || map[node.getNumber()]!=map[parent.getNumber()])){
                out = false;
                if(!verbose){
                    break;
                } else {
                    System.out.println("Node "+node.getNumber()+" failed partition connectedness check:");
                    System.out.println("Node partition: "+map[node.getNumber()].getName());
                    System.out.println("Parent partition: "+map[parent.getNumber()].getName());
                    System.out.println("Child 1 partition: "+map[firstChild.getNumber()].getName());
                    System.out.println("Child 2 partition: "+map[secondChild.getNumber()].getName());
                    System.out.println();
                }
            }
        }
        return out;
    }


    /* Return the partition of the parent of this node */

    public AbstractCase getParentCase(NodeRef node){
        return branchMap[treeModel.getParent(node).getNumber()];
    }

    /* Populates the branch map for external nodes */

    private AbstractCase[] prepareExternalNodeMap(AbstractCase[] map){
        for(int i=0; i< treeModel.getExternalNodeCount(); i++){
            TreeModel.Node currentExternalNode = (TreeModel.Node) treeModel.getExternalNode(i);
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

    private void partitionAccordingToSpecificTT(String networkFileName){
        System.out.println("Using specified starting transmission tree.");
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
            partitionAccordingToSpecificTT(specificParentMap);
        } catch(IOException e){
            throw new RuntimeException("Cannot read file: " + networkFileName );
        }
    }


    /* Takes a HashMap referring each case to its parent, and tries to paint the tree with it.
    @todo This only works on the non-extended version right now, watch it. */

    private void partitionAccordingToSpecificTT(HashMap<AbstractCase, AbstractCase> map){
        Arrays.fill(nodeRecalculationNeeded,true);
        branchMap = prepareExternalNodeMap(new AbstractCase[treeModel.getNodeCount()]);
        TreeModel.Node root = (TreeModel.Node) treeModel.getRoot();
        specificallyAssignNode(root, branchMap, map);
    }


    /* Assigns a phylogenetic tree node and its children to a partition according to a specified map of child to parent
    cases
    @todo This only works on the non-extended version right now, watch it. */

    private AbstractCase specificallyAssignNode(TreeModel.Node node, AbstractCase[] map,
                                                HashMap<AbstractCase, AbstractCase> parents){
        if(node.isExternal()){
            return map[node.getNumber()];
        } else {
            AbstractCase[] childPaintings = new AbstractCase[2];
            for(int i=0; i<node.getChildCount(); i++){
                childPaintings[i] = specificallyAssignNode(node.getChild(i), map, parents);
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
            for(int i=0; i< treeModel.getExternalNodeCount(); i++){
                TreeModel.Node extNode = (TreeModel.Node) treeModel.getExternalNode(i);
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
    * generates a random partition of the rest of the tree to start off with. If checkNonZero is true in
    * randomlyAssignNode then the network will be checked to prohibit links with zero (or rounded to zero)
    * likelihood first. @todo This always uses a non-extended partition. Check if this is OK.*/

    private void partitionAccordingToRandomTT(){
        boolean gotOne = false;
        int tries = 1;
        System.out.println("Generating a random starting painting of the tree (checking nonzero likelihood for all " +
                "branches and repeating up to 100 times until a start with nonzero likelihood is found)");
        System.out.print("Attempt: ");
        while(!gotOne){
            System.out.print(tries + "...");
            branchMap = prepareExternalNodeMap(new AbstractCase[treeModel.getNodeCount()]);
            //Warning - if the DisconnectedPartitionException in randomlyAssignNode might be caused by a bug rather than both
            //likelihoods rounding to zero, you want to stop catching this to investigate.
            try{
                partitionAccordingToRandomTT(branchMap, false);
            } catch(DisconnectedPartitionException ignored){}
            if(calculateLogLikelihood()!=Double.NEGATIVE_INFINITY){
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
    }



    /* Paints a phylogenetic tree with a random compatible painting; if checkNonZero is true, make sure all branch
    likelihoods are nonzero in the process (this sometimes still results in a zero likelihood for the whole tree, but
    is much less likely to).
    */

    private AbstractCase[] partitionAccordingToRandomTT(AbstractCase[] map, boolean checkNonZero){
        Arrays.fill(nodeRecalculationNeeded,true);
        TreeModel.Node root = (TreeModel.Node) treeModel.getRoot();
        randomlyAssignNode(root, map, checkNonZero);
        return map;
    }

    private AbstractCase randomlyAssignNode(TreeModel.Node node, AbstractCase[] map, boolean checkNonZero){
        //this makes a non-extended partition. This is OK, but if it keeps giving zero likelihoods then you could do
        //something else

        if(node.isExternal()){
            return map[node.getNumber()];
        } else {
            AbstractCase[] choices = new AbstractCase[2];
            for(int i=0; i<node.getChildCount(); i++){
                if((map[node.getChild(i).getNumber()]==null)){
                    choices[i] = randomlyAssignNode(node.getChild(i), map, checkNonZero);
                } else {
                    choices[i] = map[node.getChild(i).getNumber()];
                }
            }
            int randomSelection = MathUtils.nextInt(2);
            AbstractCase decision;
            if(checkNonZero){
                Double[] branchLogLs = new Double[2];
                for(int i=0; i<2; i++){
                    double infectionTime = getNodeTime(node)
                            - ((getNodeTime(treeModel.getParent(node)) - getNodeTime(node))
                            *infectionTimes.getParameterValue(cases.getCaseIndex(choices[1-i])));
                    branchLogLs[i]= cases.logProbXInfectedByYAtTimeT(choices[i], choices[1-i], infectionTime);
                }
                if(branchLogLs[0]==Double.NEGATIVE_INFINITY && branchLogLs[1]==Double.NEGATIVE_INFINITY){
                    throw new DisconnectedPartitionException("Both branch possibilities have zero likelihood: "
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

/*    private void debugOutputTree(){
        try{
            FlexibleTree treeCopy = new FlexibleTree(treeModel);
            for(int j=0; j<treeCopy.getNodeCount(); j++){
                FlexibleNode node = (FlexibleNode)treeCopy.getNode(j);
                node.setAttribute("Number", node.getNumber());
            }
            NexusExporter testTreesOut = new NexusExporter(new PrintStream("testTrees.nex"));
            testTreesOut.exportTree(treeCopy);
        } catch (IOException ignored) {System.out.println("IOException");}
        totalTreeLogLikelihood();
    }*/


    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        public static final String TREE_MODEL = "treeModel";
        public static final String STARTING_NETWORK = "startingNetwork";
        public static final String INFECTION_TIMES = "infectionTimes";
        public static final String EARLIEST_POSSIBLE_FIRST_INFECTION = "earliestPossibleFirstInfection";

        public String getParserName() {
            return CASE_TO_CASE_TREE_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel virusTree = (TreeModel) xo.getElementFirstChild(TREE_MODEL);

            String startingNetworkFileName=null;

            if(xo.hasChildNamed(STARTING_NETWORK)){
                startingNetworkFileName = (String) xo.getElementFirstChild(STARTING_NETWORK);
            }


            AbstractOutbreak caseSet = (AbstractOutbreak) xo.getChild(AbstractOutbreak.class);

            CaseToCaseTreeLikelihood likelihood;


            final boolean extended = xo.getBooleanAttribute("extended");

            Parameter infectionTimes = (Parameter) xo.getElementFirstChild(INFECTION_TIMES);
            Parameter earliestFirstInfection = (Parameter) xo.getElementFirstChild(EARLIEST_POSSIBLE_FIRST_INFECTION);

            try {
                likelihood = new CaseToCaseTreeLikelihood(virusTree, caseSet, startingNetworkFileName,
                        infectionTimes, earliestFirstInfection, extended);
            } catch (TaxonList.MissingTaxonException e) {
                throw new XMLParseException(e.toString());
            }

            return likelihood;
        }

        public String getParserDescription() {
            return "This element represents a likelihood function for case to case transmission.";
        }

        public Class getReturnType() {
            return CaseToCaseTreeLikelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newBooleanRule("extended"),
                new ElementRule("treeModel", Tree.class, "The tree"),
                new ElementRule(AbstractOutbreak.class, "The set of cases"),
                new ElementRule("startingNetwork", String.class, "A CSV file containing a specified starting network",
                        true)
        };
    };


    //************************************************************************
    // Loggable implementation
    //************************************************************************

    //@todo log the actual infection times, not the values of the infectionTimes parameter (probably in OMD, not here)

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

    // **************************************************************
    // TreeTraitProvider IMPLEMENTATION
    // **************************************************************

    public TreeTrait[] getTreeTraits() {
        return treeTraits.getTreeTraits();
    }

    public TreeTrait getTreeTrait(String key) {
        return treeTraits.getTreeTrait(key);
    }

    public String getNodePartition(Tree tree, NodeRef node) {
        if (tree != treeModel) {
            throw new RuntimeException("Can only reconstruct states on treeModel given to constructor");
        }
        if (!likelihoodKnown) {
            calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return branchMap[node.getNumber()].toString();
    }

    public class DisconnectedPartitionException extends RuntimeException{
        public DisconnectedPartitionException(String s){
            super(s);
        }
    }

    public class CaseToCaseLikelihoodCore extends GeneralLikelihoodCore{

        public CaseToCaseLikelihoodCore(){
            super(cases.size());
        }



        public void calculateLogLikelihoods(double[] partials, double outLogLikelihood){

            // this isn't a CTMC and doesn't have frequencies, and has only one pattern
            double[] dummyFreqArray = new double[stateCount];
            Arrays.fill(dummyFreqArray,1);
            double[] outLLArray = new double[1];
            outLLArray[0] = outLogLikelihood;

            calculateLogLikelihoods(partials, dummyFreqArray, outLLArray);
        }


    }
}








