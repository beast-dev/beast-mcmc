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
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.MathUtils;
import dr.math.distributions.ExponentialDistribution;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.*;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

import java.io.*;
import java.util.*;

/**
 * Handles manipulation of the tree partition, and likelihood of the infection times.
 *
 * @author Matthew Hall
 * @author Andrew Rambaut
 * @version $Id: $
 */

public class CaseToCaseTreeLikelihood extends AbstractTreeLikelihood implements Loggable, Citable, TreeTraitProvider {

    private final static boolean DEBUG = false;

    /* The phylogenetic tree. */

    private int noTips;

    /* Mapping of cases to branches on the tree; old version is stored before operators are applied */

    private AbstractCase[] branchMap;
    private AbstractCase[] storedBranchMap;


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
    private CaseToCaseLikelihoodCore core;

    // The tree needs to be traversed twice, first for the probability and second for the normalisation value. We need
    // a second array of flags for node updates, which corresponds to the former.

    private boolean[] updateNodeForSingleTraverse;

    // for individual partition calculations

    private double[] branchLogProbs;
    private double totalLogProb;
    private double storedTotalLogProb;

    // for normalisation

    private double[] normRootPartials;
    private double[] normProbs;
    private double normTotalProb;
    private double storedNormTotalProb;
    private boolean renormalisationNeeded;

    // where along the relevant branches the infections happen. IMPORTANT: if extended=false then this should be
    // all 0s.

    private Parameter infectionTimes;

    //because of the way the former works, we need a maximum value of the time from first infection to root node.

    private Parameter maxFirstInfToRoot;

    //@todo (experimental): mean of exponential prior on variance of infectious periods

    private ParametricDistributionModel variancePriorDistribution;

    // for extended version

    private boolean extended;

    // don't model the infectious periods

    private boolean noInfPeriodModels;

    // no need to normalise if the tree is fixed

    private boolean normalise;

    // PUBLIC STUFF

    // Name

    public static final String CASE_TO_CASE_TREE_LIKELIHOOD = "caseToCaseTreeLikelihood";
    public static final String PARTITIONS_KEY = "partition";


    // Basic constructor.

    public CaseToCaseTreeLikelihood(TreeModel virusTree, AbstractOutbreak caseData, String startingNetworkFileName,
                                    Parameter infectionTimes, ParametricDistributionModel variancePriorDistribution,
                                    Parameter maxFirstInfToRoot, boolean extended, boolean normalise,
                                    boolean noInfPeriodModels)
            throws TaxonList.MissingTaxonException {
        this(CASE_TO_CASE_TREE_LIKELIHOOD, virusTree, caseData, startingNetworkFileName, infectionTimes,
                variancePriorDistribution, maxFirstInfToRoot, extended, normalise, noInfPeriodModels);
    }

    // Constructor for an instance with a non-default name

    public CaseToCaseTreeLikelihood(String name, TreeModel virusTree, AbstractOutbreak caseData, String
            startingNetworkFileName, Parameter infectionTimes, ParametricDistributionModel variancePriorDistribution,
                                    Parameter maxFirstInfToRoot, boolean extended, boolean normalise,
                                    boolean noInfPeriodModels) {
        super(name, caseData, virusTree);

        updateNodeForSingleTraverse = new boolean[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            updateNodeForSingleTraverse[i] = true;
        }

        // @todo You should check the tree for zero branch lengths (or at least zero terminal branch lengths) and exit if you find any - it won't work

        if(stateCount!=treeModel.getExternalNodeCount()){
            throw new RuntimeException("There are duplicate tip cases.");
        }

        noTips = virusTree.getExternalNodeCount();

        cases = caseData;
        this.extended = extended;
        this.noInfPeriodModels = noInfPeriodModels;

        //@todo take this out if it doesn't work
        this.variancePriorDistribution = variancePriorDistribution;
        addModel(cases);
        verbose = false;

        this.normalise = normalise;

        Date lastSampleDate = getLatestTaxonDate();
        estimatedLastSampleTime = lastSampleDate.getTimeValue();

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
        addVariable(infectionTimes);

        if(DEBUG){
            for(int i=0; i<cases.size(); i++){
                if(!((CompoundParameter)infectionTimes).getParameter(i).getId().startsWith(cases.getCase(i).getName())){
                    throw new RuntimeException("Elements of outbreak and infectionTimes do not match up");
                }
            }
        }

        this.maxFirstInfToRoot = maxFirstInfToRoot;

        core = new CaseToCaseLikelihoodCore();
        core.initialize(nodeCount, 1, 1, false);

        for (int i = 0; i < noTips; i++) {
            // Find the id of tip i in the patternList
            String id = treeModel.getTaxonId(i);
            int index = patternList.getTaxonIndex(id);
            setStates(core, patternList, index, i);
        }
        for (int i = 0; i < treeModel.getInternalNodeCount(); i++) {
            core.createNodePartials(treeModel.getExternalNodeCount() + i);
        }


        branchLogProbs = new double[nodeCount];
        totalLogProb = 0;
        storedTotalLogProb = 0;

        normTotalProb = 0;
        storedTotalLogProb = 0;
        renormalisationNeeded = true;

        normRootPartials = new double[stateCount];
        normProbs = new double[stateCount*stateCount];


        descendantTipPartitions = new HashMap<Integer, HashSet<AbstractCase>>();
        descendantTipPartitions(virusTree.getRoot(), descendantTipPartitions);

        //paint the starting network onto the tree

        if(startingNetworkFileName==null){
            partitionAccordingToRandomTT();
        } else {
            partitionAccordingToSpecificTT(startingNetworkFileName);
        }

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

        if(DEBUG){
            treeTraits.addTrait("NodeNumber", new TreeTrait.S() {
                public String getTraitName() {
                    return "NodeNumber";
                }

                public Intent getIntent() {
                    return Intent.NODE;
                }

                public String getTrait(Tree tree, NodeRef node) {
                    return Integer.toString(node.getNumber());
                }
            });
        }

    }

    public AbstractOutbreak getOutbreak(){
        return cases;
    }

    /* Get the date of the last tip */

    private Date getLatestTaxonDate(){
        Date latestDate = new Date(new java.util.Date(-1000000), Units.Type.DAYS);
        for(AbstractCase thisCase : cases.getCases()){
            if(thisCase.getExamDate().after(latestDate)){
                latestDate = thisCase.getExamDate();
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

    /**
     * Set update flag for a node and its children
     */
    protected void updateNode(NodeRef node){
        updateNode(node, true);
    }


    protected void updateNode(NodeRef node, boolean forNormCalcToo) {

        if(forNormCalcToo){
            updateNode[node.getNumber()] = true;
        }
        updateNodeForSingleTraverse[node.getNumber()] = true;
        likelihoodKnown = false;
    }

    /**
     * Set update flag for a node and its direct children
     */

    protected void updateNodeAndChildren(NodeRef node){
        updateNodeAndChildren(node, true);
    }

    protected void updateNodeAndChildren(NodeRef node, boolean forNormCalcToo) {
        if(forNormCalcToo){
            updateNode[node.getNumber()] = true;
        }
        updateNodeForSingleTraverse[node.getNumber()] = true;

        for (int i = 0; i < treeModel.getChildCount(node); i++) {
            NodeRef child = treeModel.getChild(node, i);
            if(forNormCalcToo){
                updateNode[child.getNumber()] = true;
            }
            updateNodeForSingleTraverse[child.getNumber()] = true;
        }
        likelihoodKnown = false;
    }

    /**
     * Set update flag for a node and all its descendents
     */

    protected void updateNodeAndDescendents(NodeRef node){
        updateNodeAndDescendents(node, true);
    }

    protected void updateNodeAndDescendents(NodeRef node, boolean forNormCalcToo) {
        if(forNormCalcToo){
            updateNode[node.getNumber()] = true;
        }
        updateNodeForSingleTraverse[node.getNumber()] = true;

        for (int i = 0; i < treeModel.getChildCount(node); i++) {
            NodeRef child = treeModel.getChild(node, i);
            updateNodeAndDescendents(child, forNormCalcToo);
        }

        likelihoodKnown = false;
    }

    /**
     * Set update flag for all nodes
     */
    protected void updateAllNodes() {
        updateAllNodes(true);
    }

    protected void updateAllNodes(boolean forNormCalcToo) {
        for (int i = 0; i < nodeCount; i++) {
            if(forNormCalcToo){
                updateNode[i] = true;
            }
            updateNodeForSingleTraverse[i] = true;
        }
        likelihoodKnown = false;
    }



    // ************************************************************************************
    // EXTENDED VERSION METHODS
    // ************************************************************************************

    /* check if the given node is tip-linked under the current painting (the tip corresponding to its painting is
    a descendant of it
     */

    public boolean tipLinked(NodeRef node){
        return tipLinked(node, branchMap);
    }

    private boolean tipLinked(NodeRef node, AbstractCase[] map){
        NodeRef tip = treeModel.getNode(tipMap.get(map[node.getNumber()]));
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



    public int countChildrenInSamePartition(NodeRef node, AbstractCase[] map){
        if(treeModel.isExternal(node)){
            return -1;
        } else {
            int count = 0;
            AbstractCase parentCase = map[node.getNumber()];
            for(int i=0; i< treeModel.getChildCount(node); i++){
                if(map[treeModel.getChild(node,i).getNumber()]==parentCase){
                    count++;
                }
            }
            return count;
        }
    }

    public int countChildrenInSamePartition(NodeRef node){
        return countChildrenInSamePartition(node, branchMap);
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
        return samePartitionDownTree(node, branchMap, flagForRecalc);
    }

    private HashSet<Integer> samePartitionDownTree(NodeRef node, AbstractCase[] map, boolean flagForRecalc){
        if(!updateNodeForSingleTraverse[node.getNumber()] && flagForRecalc){
            flagForDescendantRecalculation(treeModel, node, true);
        }
        HashSet<Integer> out = new HashSet<Integer>();
        AbstractCase painting = map[node.getNumber()];
        NodeRef currentNode = node;
        NodeRef parentNode = treeModel.getParent(node);
        while(parentNode!=null && map[parentNode.getNumber()]==painting){
            out.add(parentNode.getNumber());
            if(extended){
                if(countChildrenInSamePartition(parentNode)==2){
                    NodeRef otherChild = sibling(treeModel, currentNode);
                    out.add(otherChild.getNumber());
                    out.addAll(samePartitionUpTree(otherChild, map, flagForRecalc));
                }
            }
            currentNode = parentNode;
            parentNode = treeModel.getParent(currentNode);
        }
        return out;
    }

    //Return a set of nodes that are descendants (and not equal to) the current node and are in the same parition as it.

    public HashSet<Integer> samePartitionUpTree(NodeRef node, boolean flagForRecalc){
        return samePartitionUpTree(node, branchMap, flagForRecalc);
    }

    private HashSet<Integer> samePartitionUpTree(NodeRef node, AbstractCase[] map, boolean flagForRecalc){
        HashSet<Integer> out = new HashSet<Integer>();
        AbstractCase painting = map[node.getNumber()];
        boolean creepsFurther = false;
        for(int i=0; i< treeModel.getChildCount(node); i++){
            if(map[treeModel.getChild(node,i).getNumber()]==painting){
                creepsFurther = true;
                out.add(treeModel.getChild(node,i).getNumber());
                out.addAll(samePartitionUpTree(treeModel.getChild(node, i), map, flagForRecalc));
            }
        }
        if(flagForRecalc && !creepsFurther){
            flagForDescendantRecalculation(treeModel, node, true);
        }
        return out;
    }

    // Return the node numbers of the entire subtree in the same partition as this node (including itself)

    public Integer[] samePartition(NodeRef node, boolean flagForRecalc){
        return samePartition(node, branchMap, flagForRecalc);
    }

    private Integer[] samePartition(NodeRef node, AbstractCase[] map, boolean flagForRecalc){
        HashSet<Integer> out = new HashSet<Integer>();
        out.add(node.getNumber());
        out.addAll(samePartitionDownTree(node, map, flagForRecalc));
        out.addAll(samePartitionUpTree(node, map, flagForRecalc));
        return out.toArray(new Integer[out.size()]);
    }

    public void changeMap(int node, AbstractCase partition){
        branchMap[node]=partition;
        likelihoodKnown = false;
        // @todo you could get efficiency savings here
        Arrays.fill(updateNodeForSingleTraverse, true);
        fireModelChanged();
    }

    public Integer[] getParentsArray(){
        Integer[] out = new Integer[cases.size()];
        for(AbstractCase thisCase : cases.getCases()){
            out[cases.getCaseIndex(thisCase)]=cases.getCaseIndex(getInfector(thisCase));
        }
        return out;
    }

    public HashMap<Integer, Double> getInfTimesMap(){
        HashMap<Integer, Double> out = new HashMap<Integer, Double>();
        for(AbstractCase thisCase : cases.getCases()){
            out.put(cases.getCaseIndex(thisCase), getInfectionTime(thisCase));
        }
        return out;
    }

    // change flags to indicate that something needs recalculation further down the tree

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

    public void flagForDescendantRecalculation(TreeModel tree, NodeRef node, boolean normCalcsToo){
        flagForDescendantRecalculation(tree, node, updateNodeForSingleTraverse);
        if(normCalcsToo){
            flagForDescendantRecalculation(tree, node, updateNode);
        }
    }

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************


    protected final void handleModelChangedEvent(Model model, Object object, int index) {

        if(model == cases){
            Arrays.fill(updateNode, true);
            Arrays.fill(updateNodeForSingleTraverse, true);
            renormalisationNeeded = true;
        }

        if (model == treeModel) {
            if(!extended){
                recalculateDescTipPartitions();
            }
            if (object instanceof TreeModel.TreeChangedEvent) {

                if (((TreeModel.TreeChangedEvent) object).isNodeChanged()) {
                    // If a node event occurs the node and its two child nodes
                    // are flagged for updating (this will result in everything
                    // above being updated as well. Node events occur when a node
                    // is added to a branch, removed from a branch or its height or
                    // rate changes.
                    updateNodeAndChildren(((TreeModel.TreeChangedEvent) object).getNode());
                } else if (((TreeModel.TreeChangedEvent) object).isTreeChanged()) {
                    // Full tree events result in a complete updating of the tree likelihood
                    updateAllNodes();
                } else {
                    // Other event types are ignored (probably trait changes).
                    //System.err.println("Another tree event has occured (possibly a trait change).");
                }
            }
            renormalisationNeeded = true;
        }

        fireModelChanged();

        likelihoodKnown = false;
    }


    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************


    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

        updateAllNodes(variable!=infectionTimes);

        renormalisationNeeded = variable!=infectionTimes;

        fireModelChanged();

        likelihoodKnown = false;
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    /**
     * Stores the precalculated state (in this case the node labels and subtree likelihoods)
     */

    protected final void storeState() {
        super.storeState();
        core.storeState();
        storedNormTotalProb = normTotalProb;
        storedTotalLogProb = totalLogProb;
        storedBranchMap = Arrays.copyOf(branchMap, branchMap.length);
        storedDescendantTipPartitions = new HashMap<Integer, HashSet<AbstractCase>>(descendantTipPartitions);
    }

    /**
     * Restores the precalculated state.
     */

    protected final void restoreState() {
        super.restoreState();
        core.restoreState();
        normTotalProb = storedNormTotalProb;
        totalLogProb = storedTotalLogProb;
        branchMap = storedBranchMap;
        descendantTipPartitions = storedDescendantTipPartitions;
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
        likelihoodKnown = false;
        // @todo you could get efficiency savings here
        Arrays.fill(updateNodeForSingleTraverse, true);
        fireModelChanged();
    }

    public final TreeModel getTree(){
        return treeModel;
    }

    public final void makeDirty() {
        likelihoodKnown = false;
        Arrays.fill(updateNode, true);
        Arrays.fill(updateNodeForSingleTraverse, true);
        renormalisationNeeded = true;
        if(!extended){
            recalculateDescTipPartitions();
        }
    }

    /**
     * Calculates the log likelihood of this set of node labels given the tree.
     */
    protected double calculateLogLikelihood() {

        double variance = getSummaryStatistics()[2];
        double logVariancePriorProb = variancePriorDistribution.logPdf(variance);

        if(noInfPeriodModels){
            return isAllowed() ? (0+logVariancePriorProb) : Double.NEGATIVE_INFINITY;
        }

        if(DEBUG && !checkPartitions(branchMap, true)){
            throw new RuntimeException("Partition rules are violated");
        }

        final NodeRef root = treeModel.getRoot();

        // unnormalised probability

        traverse(treeModel, root);

        // normalisation value
        if(normalise){
            if(renormalisationNeeded){
                traverseForNormalisation(treeModel, root);
            }
        } else {
            normTotalProb = 0;
        }

        double logL = totalLogProb - normTotalProb;

        // If logL is zero because the partition is actually impossible, you most certainly don't want to turn scaling
        // on, hence the okToRescale

        if (renormalisationNeeded && normTotalProb == Double.NEGATIVE_INFINITY) {
//            Logger.getLogger("dr.evomodel").info("TreeLikelihood, " + this.getId() + ", turning on partial " +
//                    "likelihood scaling to avoid precision loss");

            // We probably had an underflow... turn on scaling
            core.setUseScaling(true);

            // and try again...
            updateAllNodes();
            updateAllPatterns();
            traverseForNormalisation(treeModel, root);

            logL = totalLogProb - (normTotalProb) * patternWeights[0];

        }


        for (int i = 0; i < nodeCount; i++) {
            updateNode[i] = false;
            updateNodeForSingleTraverse[i] = false;
        }

        renormalisationNeeded = false;

        // If the normalisation value rounds to zero it is an almighty pain, but perhaps it can't be helped

        if(logL==Double.POSITIVE_INFINITY){
            return Double.NEGATIVE_INFINITY;
        }

        return logL;

    }

    // if no infectious models, just need to check whether any infections occur after the infector was no
    // longer infectious

    private boolean isAllowed(){
        return isAllowed(treeModel.getRoot());
    }

    private boolean isAllowed(NodeRef node){
        if(!treeModel.isRoot(node)){
            AbstractCase childCase = branchMap[node.getNumber()];
            AbstractCase parentCase = branchMap[treeModel.getParent(node).getNumber()];
            if(childCase!=parentCase){
                double infectionTime = getInfectionTime(childCase);
                if(infectionTime>=parentCase.getCullDate().getTimeValue()){
                    return false;
                }
            }
        }
        if(!treeModel.isExternal(node)){
            return isAllowed(treeModel.getChild(node,0)) && isAllowed(treeModel.getChild(node,1));
        } else {
            return true;
        }
    }

    // unnormalised probability of this branchMap and corresponding infTimes

    private boolean traverse(Tree tree, NodeRef node){

        boolean update = false;

        int nodeNum = node.getNumber();

        NodeRef parent = tree.getParent(node);

        if (updateNodeForSingleTraverse[nodeNum]) {
            if(parent!=null){
                AbstractCase nodeCase = branchMap[node.getNumber()];
                AbstractCase parentCase = branchMap[parent.getNumber()];
                if(nodeCase!=parentCase){
                    double infectionTime = getInfectionTime(nodeCase);
                    branchLogProbs[node.getNumber()]
                            = cases.logProbXInfectedByYAtTimeT(nodeCase, parentCase, infectionTime);
                } else {
                    branchLogProbs[node.getNumber()] = 0;
                }
            } else {
                AbstractCase nodeCase = branchMap[node.getNumber()];
                final double infectionTime = getRootInfectionTime();
                branchLogProbs[node.getNumber()]
                        = cases.logProbXInfectedAtTimeT(nodeCase, infectionTime);
            }
            update = true;
        }
        if (!tree.isExternal(node)) {
            NodeRef child1 = tree.getChild(node, 0);
            final boolean update1 = traverse(tree, child1);
            NodeRef child2 = tree.getChild(node, 1);
            final boolean update2 = traverse(tree, child2);
            if (update1 || update2) {
                if (parent == null) {
                    totalLogProb = 0;
                    for(int i=0; i<nodeCount; i++){
                        totalLogProb += branchLogProbs[i];
                    }
                }
                update = true;
            }
        }
        return(update);
    }



    // normalisation value

    private boolean traverseForNormalisation(Tree tree, NodeRef node){

        if(DEBUG){
            debugOutputTree("traverseProblem.nex");
        }

        boolean update = false;

        int nodeNum = node.getNumber();

        NodeRef parent = tree.getParent(node);

        if (updateNode[nodeNum] && parent!=null) {

            if(treeModel.isExternal(node)){
                // vast majority of entries are zero in this case
                Arrays.fill(normProbs, 0);
                // kind of hacky to use branchMap here to be honest
                AbstractCase destination = branchMap[node.getNumber()];
                int j = cases.getCaseIndex(destination);
                for(int i=0; i<stateCount; i++){
                    AbstractCase origin = cases.getCase(i);
                    // always compatible if non-extended
                    boolean treeCompatibilityCheck = true;
                    if(!extended){
                        HashSet<AbstractCase> parentDescTips = descendantTipPartitions.get(parent.getNumber());
                        treeCompatibilityCheck = parentDescTips.contains(origin);
                    }
                    if(!treeCompatibilityCheck){
                        normProbs[stateCount*i + j]=0;
                    } else if (origin==destination){
                        normProbs[stateCount*i + j]=1;
                    } else {
                        if(!extended){
                            normProbs[stateCount*i + j]=cases.probXInfectedByYAtTimeT(destination, origin,
                                    getNodeTime(parent));
                        } else {
                            normProbs[stateCount*i + j]=cases.probXInfectedByYBetweenTandU(destination, origin,
                                    getNodeTime(parent), getNodeTime(node));
                        }
                    }
                }
            } else {
                HashSet<AbstractCase> nodeDescTips = descendantTipPartitions.get(node.getNumber());
                HashSet<AbstractCase> parentDescTips = descendantTipPartitions.get(parent.getNumber());
                for(int i=0; i<stateCount; i++){
                    AbstractCase origin = cases.getCase(i);
                    for(int j=0; j<stateCount; j++){
                        AbstractCase destination = cases.getCase(j);

                        // is the tip in parent's partition a descendant of this node? If so, the node must be in
                        // this partition also
                        boolean paintingForcedByParent = nodeDescTips.contains(origin);
                        boolean treeCompatibilityCheck;
                        if(!extended){
                            // valid combinations:
                            // 1) paintingForcedByParent = false
                            // 2) paintingForcedByParent = true and both nodes in same partition
                            treeCompatibilityCheck = nodeDescTips.contains(destination)
                                    && parentDescTips.contains(origin)
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
                            treeCompatibilityCheck = ((!nodeCreep && !parentCreep
                                    && (!paintingForcedByParent || origin == destination))
                                    || (parentCreep && !nodeCreep) || (nodeCreep && (origin == destination)));
                        }
                        if(!treeCompatibilityCheck){
                            normProbs[stateCount*i + j]=0;
                        } else if(origin==destination) {
                            normProbs[stateCount*i + j]=1;
                        } else {
                            if(!extended){
                                normProbs[stateCount*i + j]=cases.probXInfectedByYAtTimeT(destination, origin,
                                        getNodeTime(parent));
                            } else {
                                normProbs[stateCount*i + j]=cases.probXInfectedByYBetweenTandU(destination, origin,
                                        getNodeTime(parent), getNodeTime(node));
                            }
                        }
                    }

                }
            }
            core.setNodeMatrix(nodeNum, 0, normProbs);
            update = true;
        }
        if (!tree.isExternal(node)) {

            // Traverse down the two child nodes
            NodeRef child1 = tree.getChild(node, 0);
            final boolean update1 = traverseForNormalisation(tree, child1);

            NodeRef child2 = tree.getChild(node, 1);
            final boolean update2 = traverseForNormalisation(tree, child2);

            if (update1 || update2) {

                final int childNum1 = child1.getNumber();
                final int childNum2 = child2.getNumber();

                core.setNodePartialsForUpdate(nodeNum);

                core.calculatePartials(childNum1, childNum2, nodeNum);

                if (parent == null) {
                    // No parent so this is the root of the tree -
                    // calculate the pattern likelihoods

                    double[] partials = getRootPartials();

                    core.calculateLogLikelihoods(partials);
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

            normRootPartials[i] *= cases.probXInfectedBetweenTandU
                    (cases.getCase(i), getNodeTime(treeModel.getRoot()) - maxFirstInfToRoot.getParameterValue(0),
                            getNodeTime(treeModel.getRoot()));
        }


        return normRootPartials;
    }

    /* Return the double time at which the given node occurred */

    public double getNodeTime(NodeRef node){
        double nodeHeight = getHeight(node);
        return estimatedLastSampleTime-nodeHeight;
    }

    public double heightToTime(double height){
        return estimatedLastSampleTime-height;
    }

    public double timeToHeight(double time){
        return estimatedLastSampleTime-time;
    }

    private double getHeight(NodeRef node){
        return treeModel.getNodeHeight(node);
    }

    /* Return the case which infected this case */

    public AbstractCase getInfector(AbstractCase thisCase){
        return getInfector(thisCase, branchMap);
    }

    public AbstractCase getInfector(int i){
        return getInfector(getOutbreak().getCase(i), branchMap);
    }

    /* Return the case which was the infector in the infection event represented by this node */

    public AbstractCase getInfector(NodeRef node){
        return getInfector(node, branchMap);
    }

    public AbstractCase getInfector(AbstractCase thisCase, AbstractCase[] branchMap){
        NodeRef tip = treeModel.getNode(tipMap.get(thisCase));
        return getInfector(tip, branchMap);
    }

    public double getInfectionTime(AbstractCase thisCase, AbstractCase[] branchMap){
        NodeRef child = treeModel.getNode(tipMap.get(thisCase));
        NodeRef parent = treeModel.getParent(child);
        boolean transmissionFound = false;
        while(!transmissionFound){
            if(branchMap[child.getNumber()]!=branchMap[parent.getNumber()]){
                transmissionFound = true;
            } else {
                child = parent;
                parent = treeModel.getParent(child);
                if(parent == null){
                    transmissionFound = true;
                }
            }
        }
        if(parent!=null){
            AbstractCase parentCase = branchMap[parent.getNumber()];
            double min = heightToTime(treeModel.getNodeHeight(parent));

            // the infection must have taken place on this branch, and before the cull of the parent case

            double max = Math.min(heightToTime(treeModel.getNodeHeight(child)),
                    parentCase.getCullDate().getTimeValue());


            return getInfectionTime(min, max, thisCase);
        } else {
            return getRootInfectionTime(branchMap);
        }
    }

    public double getInfectiousPeriod(AbstractCase thisCase){
        return getInfectiousPeriod(thisCase, branchMap);
    }

    public double getInfectiousPeriod(AbstractCase thisCase, AbstractCase[] branchMap){
        double infectionTime = getInfectionTime(thisCase, branchMap);
        double cullTime = thisCase.getCullTime();
        return cullTime - infectionTime;
    }

    private double getInfectionTime(double min, double max, AbstractCase infected){
        final double branchLength = max-min;

        return max - branchLength*infectionTimes.getParameterValue(cases.getCaseIndex(infected));
    }

    // return an array of the mean, median, variance and standard deviation of infectious periods

    public double[] getSummaryStatistics(){
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for(AbstractCase thisCase: cases.getCases()){
            stats.addValue(getInfectiousPeriod(thisCase));
        }
        double[] out = new double[4];
        out[0] = stats.getMean();
        out[1] = stats.getPercentile(50);
        out[2] = stats.getVariance();
        out[3] = stats.getStandardDeviation();
        return out;
    }


    private double getRootInfectionTime(AbstractCase[] branchMap){
        NodeRef root = treeModel.getRoot();
        AbstractCase rootCase = branchMap[root.getNumber()];
        final double branchLength = maxFirstInfToRoot.getParameterValue(0);

        return heightToTime(treeModel.getNodeHeight(root)
                + branchLength*infectionTimes.getParameterValue(cases.getCaseIndex(rootCase)));

    }

    private double getRootInfectionTime(){
        return getRootInfectionTime(branchMap);
    }

    public double getInfectionTime(AbstractCase thisCase){
        return getInfectionTime(thisCase, branchMap);
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
        boolean foundProblem = false;
        for(int i=0; i<treeModel.getInternalNodeCount(); i++){
            boolean foundTip = false;
            for(Integer nodeNumber : samePartition(treeModel.getInternalNode(i), map, false)){
                if(treeModel.isExternal(treeModel.getNode(nodeNumber))){
                    foundTip = true;
                }
            }
            if(!foundProblem && !foundTip){
                foundProblem = true;
                if(verbose){
                    System.out.println("Node "+(i+treeModel.getExternalNodeCount()) + " is not connected to a tip");
                }
            }

        }
        if(foundProblem){
            debugOutputTree(map, "checkPartitionProblem");
        }
        return !foundProblem;
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

/*  The CSV file should have a header, and then lines matching each case to its infector*/

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

    private void partitionAccordingToSpecificTT(HashMap<AbstractCase, AbstractCase> map){
        branchMap = prepareExternalNodeMap(new AbstractCase[treeModel.getNodeCount()]);

        AbstractCase firstCase=null;
        for(AbstractCase aCase : cases.getCases()){
            if(map.get(aCase)==null){
                firstCase = aCase;
            }
        }
        if(firstCase==null){
            throw new RuntimeException("Given starting network is not compatible with the starting tree");
        }
        NodeRef root = treeModel.getRoot();
        specificallyPartitionUpwards(root, firstCase, map);
        if(!checkPartitions()){
            throw new RuntimeException("Given starting network is not compatible with the starting tree");
        }

    }

    private void specificallyPartitionUpwards(NodeRef node, AbstractCase thisCase,
                                              HashMap<AbstractCase, AbstractCase> map){
        if(treeModel.isExternal(node)){
            return;
        }
        branchMap[node.getNumber()]=thisCase;
        if(tipLinked(node)){
            for(int i=0; i<treeModel.getChildCount(node); i++){
                specificallyPartitionUpwards(treeModel.getChild(node, i), thisCase, map);
            }
        } else {
            branchMap[node.getNumber()]=null;
            HashSet<AbstractCase> children = new HashSet<AbstractCase>();
            for(AbstractCase aCase : cases.getCases()){
                if(map.get(aCase)==thisCase){
                    children.add(aCase);
                }
            }
            HashSet<AbstractCase> relevantChildren = new HashSet<AbstractCase>(children);
            for(AbstractCase child: children){
                int tipNo = tipMap.get(child);
                NodeRef currentNode = treeModel.getExternalNode(tipNo);
                while(currentNode!=node && currentNode!=null){
                    currentNode = treeModel.getParent(currentNode);
                }
                if(currentNode==null){
                    relevantChildren.remove(child);
                }
            }
            if(relevantChildren.size()==1){
                //no creep
                AbstractCase child = relevantChildren.iterator().next();
                branchMap[node.getNumber()]=child;
            } else {
                branchMap[node.getNumber()]=thisCase;
            }
            for(int i=0; i<treeModel.getChildCount(node); i++){
                specificallyPartitionUpwards(treeModel.getChild(node, i), branchMap[node.getNumber()], map);
            }
        }

    }


    /* Assigns a phylogenetic tree node and its children to a partition according to a specified map of child to parent
    cases. This only works on the non-extended version right now, watch it. */

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
    * likelihood first. This always uses a non-extended partition. */

    private void partitionAccordingToRandomTT(){
        boolean gotOne = false;
        int tries = 1;
        System.out.println("Generating a random starting painting of the tree (checking nonzero likelihood for all " +
                "branches and repeating up to 100 times until a start with nonzero likelihood is found)");
        System.out.print("Attempt: ");
        while(!gotOne){
            boolean failed = false;
            System.out.print(tries + "...");
            branchMap = prepareExternalNodeMap(new AbstractCase[treeModel.getNodeCount()]);
            //Warning - if the BadPartitionException in randomlyAssignNode might be caused by a bug rather than both
            //likelihoods rounding to zero, you want to stop catching this to investigate.

            try{
                partitionAccordingToRandomTT(branchMap);
            } catch(BadPartitionException e){
                failed = true;
            }

            if(!failed && calculateLogLikelihood()!=Double.NEGATIVE_INFINITY){
                gotOne = true;
                System.out.println("found.");
            }
            tries++;
            if(tries==101){
                System.out.println("giving " +
                        "up.");
                throw new RuntimeException("Failed to find a starting transmission tree with nonzero likelihood");
            }
        }
    }



    /* Paints a phylogenetic tree with a random compatible painting; if checkNonZero is true, make sure all branch
    likelihoods are nonzero in the process (this sometimes still results in a zero likelihood for the whole tree, but
    is much less likely to).
    */

    private AbstractCase[] partitionAccordingToRandomTT(AbstractCase[] map){
        Arrays.fill(updateNode, true);
        Arrays.fill(updateNodeForSingleTraverse, true);
        TreeModel.Node root = (TreeModel.Node) treeModel.getRoot();
        randomlyAssignNode(root, map);
        return map;
    }

    private AbstractCase randomlyAssignNode(TreeModel.Node node, AbstractCase[] map){
        //this makes a non-extended partition. This is OK, but if it keeps giving zero likelihoods then you could do
        //something else

        if(node.isExternal()){
            return map[node.getNumber()];
        } else {
            AbstractCase[] choices = new AbstractCase[2];
            for(int i=0; i<node.getChildCount(); i++){
                if((map[node.getChild(i).getNumber()]==null)){
                    choices[i] = randomlyAssignNode(node.getChild(i), map);
                } else {
                    choices[i] = map[node.getChild(i).getNumber()];
                }
            }
            int randomSelection = MathUtils.nextInt(2);
            int decision;
            Double[] branchLogLs = new Double[2];
            for(int i=0; i<2; i++){
                branchLogLs[i]= cases.logProbXInfectedByYBetweenTandU(choices[1-i], choices[i],
                        getNodeTime(node), getNodeTime(treeModel.getChild(node, 1-i)));
            }
            if(branchLogLs[0]==Double.NEGATIVE_INFINITY && branchLogLs[1]==Double.NEGATIVE_INFINITY){
                throw new BadPartitionException("Both branch possibilities have zero likelihood: "
                        +node.toString()+", cases " + choices[0].getName() + " and " + choices[1].getName() + ".");
            } else if(branchLogLs[0]==Double.NEGATIVE_INFINITY || branchLogLs[1]==Double.NEGATIVE_INFINITY){
                if(branchLogLs[0]==Double.NEGATIVE_INFINITY){
                    decision = 1;
                } else {
                    decision = 0;
                }
            } else {
                decision = randomSelection;
            }
            AbstractCase winner = choices[decision];
            map[node.getNumber()]=winner;
            return winner;
        }
    }

    public void debugOutputTree(String fileName){
        debugOutputTree(branchMap, fileName);
    }


    public void debugOutputTree(AbstractCase[] map, String fileName){
        try{
            FlexibleTree treeCopy = new FlexibleTree(treeModel);
            for(int j=0; j<treeCopy.getNodeCount(); j++){
                FlexibleNode node = (FlexibleNode)treeCopy.getNode(j);
                node.setAttribute("Number", node.getNumber());
                node.setAttribute("Time", heightToTime(node.getHeight()));
                node.setAttribute("Partition", map[node.getNumber()]);
            }
            NexusExporter testTreesOut = new NexusExporter(new PrintStream(fileName));
            testTreesOut.exportTree(treeCopy);
        } catch (IOException ignored) {System.out.println("IOException");}
    }




    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        public static final String STARTING_NETWORK = "startingNetwork";
        public static final String INFECTION_TIMES = "infectionTimes";
        public static final String MAX_FIRST_INF_TO_ROOT= "maxFirstInfToRoot";
        public static final String EXTENDED = "extended";
        public static final String NORMALISE = "normalise";
        public static final String NO_INF_PERIOD_MODELS = "noInfPeriodModels";
        public static final String VARIANCE_PRIOR_DISTRIBUTION = "variancePriorDistribution";

        public String getParserName() {
            return CASE_TO_CASE_TREE_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel virusTree = (TreeModel) xo.getChild(TreeModel.class);

            String startingNetworkFileName=null;

            if(xo.hasChildNamed(STARTING_NETWORK)){
                startingNetworkFileName = (String) xo.getElementFirstChild(STARTING_NETWORK);
            }


            AbstractOutbreak caseSet = (AbstractOutbreak) xo.getChild(AbstractOutbreak.class);

            CaseToCaseTreeLikelihood likelihood;


            final boolean extended = xo.getBooleanAttribute(EXTENDED);

            final boolean normalise = xo.getBooleanAttribute(NORMALISE);

            final boolean noInfPeriodModels = xo.getBooleanAttribute(NO_INF_PERIOD_MODELS);


            ParametricDistributionModel variancePriorDistribution = (ParametricDistributionModel)
                    xo.getElementFirstChild(VARIANCE_PRIOR_DISTRIBUTION);
            Parameter infectionTimes = (Parameter) xo.getElementFirstChild(INFECTION_TIMES);
            Parameter earliestFirstInfection = (Parameter) xo.getElementFirstChild(MAX_FIRST_INF_TO_ROOT);

            try {
                likelihood = new CaseToCaseTreeLikelihood(virusTree, caseSet, startingNetworkFileName, infectionTimes,
                        variancePriorDistribution, earliestFirstInfection, extended, normalise, noInfPeriodModels);
            } catch (TaxonList.MissingTaxonException e) {
                throw new XMLParseException(e.toString());
            }

            return likelihood;
        }

        public String getParserDescription() {
            return "This element represents a probability distribution for the infection dates for cases of an outbreak"
                    +"given a phylogenetic tree";
        }

        public Class getReturnType() {
            return CaseToCaseTreeLikelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newBooleanRule(EXTENDED),
                AttributeRule.newBooleanRule(NORMALISE),
                AttributeRule.newBooleanRule(NO_INF_PERIOD_MODELS),
                new ElementRule(TreeModel.class, "The tree"),
                new ElementRule(AbstractOutbreak.class, "The set of cases"),
                new ElementRule("startingNetwork", String.class, "A CSV file containing a specified starting network",
                        true),
                new ElementRule(MAX_FIRST_INF_TO_ROOT, Parameter.class, "The maximum time from the first infection to" +
                        "the root node"),
                new ElementRule(INFECTION_TIMES, Parameter.class),
                new ElementRule(VARIANCE_PRIOR_DISTRIBUTION, ParametricDistributionModel.class)
        };
    };


    //************************************************************************
    // Loggable implementation
    //************************************************************************


    public LogColumn[] getColumns(){
        LogColumn[] columns = new LogColumn[cases.size()];
        for(int i=0; i<cases.size(); i++){
            final AbstractCase infected = cases.getCase(i);
            columns[i] = new LogColumn.Abstract(infected.toString()+"_infector"){
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

    public class CaseToCaseLikelihoodCore extends GeneralLikelihoodCore{

        public CaseToCaseLikelihoodCore(){
            super(cases.size());
        }

        public void calculateLogLikelihoods(double[] partials){

            // this isn't a CTMC and doesn't have frequencies, and has only one pattern
            double[] dummyFreqArray = new double[stateCount];
            Arrays.fill(dummyFreqArray,1);
            double[] outLLArray = new double[1];
            calculateLogLikelihoods(partials, dummyFreqArray, outLLArray);
            normTotalProb = outLLArray[0];
        }

    }
}








