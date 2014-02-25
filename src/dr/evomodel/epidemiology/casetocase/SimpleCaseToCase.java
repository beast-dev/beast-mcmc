package dr.evomodel.epidemiology.casetocase;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.GeneralLikelihoodCore;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.*;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by mhall on 17/12/2013.
 */
public class SimpleCaseToCase extends CaseToCaseTreeLikelihood {

    public static final String SIMPLE_CASE_TO_CASE = "simpleCaseToCase";

    /* the partitions that all descendant tips of a current node belong to */

    private HashMap<Integer, HashSet<AbstractCase>> descendantTipPartitions;
    private HashMap<Integer, HashSet<AbstractCase>> storedDescendantTipPartitions;

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


    private boolean traversalProbKnown = false;
    private boolean storedTraversalProbKnown = false;


    public SimpleCaseToCase(PartitionedTreeModel virusTree, AbstractOutbreak caseData,
                            String startingNetworkFileName, Parameter infectionTimeBranchPositions,
                            Parameter maxFirstInfToRoot) {
        super(SIMPLE_CASE_TO_CASE, virusTree, caseData, infectionTimeBranchPositions, null, maxFirstInfToRoot);
        if(caseData.hasLatentPeriods()){
            throw new RuntimeException("Latent periods not currently implemeted for SimpleCaseToCase");
        }
        updateNodeForSingleTraverse = new boolean[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            updateNodeForSingleTraverse[i] = true;
        }


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

        prepareTree(startingNetworkFileName);
    }

    protected double calculateLogLikelihood() {
        double treeLogL;

        if(!traversalProbKnown){

            if(DEBUG && !checkPartitions(branchMap, true)){
                throw new RuntimeException("Partition rules are violated");
            }

            if(DEBUG){
                debugOutputTree("treeTest.nex", false);
            }

            final NodeRef root = treeModel.getRoot();

            // unnormalised probability

            traverse(treeModel, root);

            // normalisation value

            if(renormalisationNeeded){
                traverseForNormalisation(treeModel, root);
            }


            treeLogL =  totalLogProb - normTotalProb;

            if (renormalisationNeeded && normTotalProb == Double.NEGATIVE_INFINITY) {
                //            Logger.getLogger("dr.evomodel").info("TreeLikelihood, " + this.getId() + ", turning on partial " +
                //                    "likelihood scaling to avoid precision loss");

                // We probably had an underflow... turn on scaling
                core.setUseScaling(true);

                // and try again...
                updateAllNodes();
                updateAllPatterns();
                traverseForNormalisation(treeModel, root);

                treeLogL = totalLogProb - (normTotalProb) * patternWeights[0];

            }


            for (int i = 0; i < nodeCount; i++) {
                updateNode[i] = false;
                updateNodeForSingleTraverse[i] = false;
            }

            renormalisationNeeded = false;

            // If the normalisation value rounds to zero it is an almighty pain, but perhaps it can't be helped

            if(treeLogL==Double.POSITIVE_INFINITY){
                return Double.NEGATIVE_INFINITY;
            }
        } else {
            treeLogL =  totalLogProb - normTotalProb;
        }

        return treeLogL;

    }

    //   unnormalised probability of this branchMap and corresponding infTimes

    private boolean traverse(Tree tree, NodeRef node){

        boolean update = false;

        int nodeNum = node.getNumber();

        NodeRef parent = tree.getParent(node);

        if (updateNodeForSingleTraverse[nodeNum]) {
            if(parent!=null){
                AbstractCase nodeCase = branchMap.get(node.getNumber());
                AbstractCase parentCase = branchMap.get(parent.getNumber());
                if(nodeCase!=parentCase){
                    double infectionTime = infectionTimes[outbreak.getCaseIndex(nodeCase)];
                    branchLogProbs[node.getNumber()]
                            = ((SimpleOutbreak) outbreak).logProbXInfectedByYAtTimeT(nodeCase, parentCase, infectionTime);
                } else {
                    branchLogProbs[node.getNumber()] = 0;
                }
            } else {
                AbstractCase nodeCase = branchMap.get(node.getNumber());
                final double infectionTime = getRootInfectionTime();
                branchLogProbs[node.getNumber()]
                        = ((SimpleOutbreak) outbreak).logProbXInfectedAtTimeT(nodeCase, infectionTime);
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

        boolean update = false;

        int nodeNum = node.getNumber();

        NodeRef parent = tree.getParent(node);

        if (updateNode[nodeNum] && parent!=null) {

            if(treeModel.isExternal(node)){
                // vast majority of entries are zero in this case
                Arrays.fill(normProbs, 0);
                // kind of hacky to use branchMap here to be honest
                AbstractCase destination = branchMap.get(node.getNumber());
                int j = outbreak.getCaseIndex(destination);
                for(int i=0; i<stateCount; i++){
                    AbstractCase origin = outbreak.getCase(i);
                    // always compatible if non-extended
                    boolean treeCompatibilityCheck = true;
                    if(!treeCompatibilityCheck){
                        normProbs[stateCount*i + j]=0;
                    } else if (origin==destination){
                        normProbs[stateCount*i + j]=1;
                    } else {
                        normProbs[stateCount*i + j]=
                                ((SimpleOutbreak) outbreak).probXInfectedByYBetweenTandU(destination, origin,
                                        getNodeTime(parent), getNodeTime(node));
                    }
                }
            } else {
                HashSet<AbstractCase> nodeDescTips = descendantTipPartitions.get(node.getNumber());
                HashSet<AbstractCase> parentDescTips = descendantTipPartitions.get(parent.getNumber());
                for(int i=0; i<stateCount; i++){
                    AbstractCase origin = outbreak.getCase(i);
                    for(int j=0; j<stateCount; j++){
                        AbstractCase destination = outbreak.getCase(j);

                        // is the tip in parent's partition a descendant of this node? If so, the node must be in
                        // this partition also
                        boolean paintingForcedByParent = nodeDescTips.contains(origin);
                        boolean treeCompatibilityCheck;

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

                        if(!treeCompatibilityCheck){
                            normProbs[stateCount*i + j]=0;
                        } else if(origin==destination) {
                            normProbs[stateCount*i + j]=1;
                        } else {
                            normProbs[stateCount*i + j]=
                                    ((SimpleOutbreak) outbreak).probXInfectedByYBetweenTandU(destination, origin,
                                            getNodeTime(parent), getNodeTime(node));
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

            normRootPartials[i] *= ((SimpleOutbreak) outbreak).probXInfectedBetweenTandU
                    (outbreak.getCase(i), getNodeTime(treeModel.getRoot()) - maxFirstInfToRoot.getParameterValue(0),
                            getNodeTime(treeModel.getRoot()));
        }


        return normRootPartials;
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



    public void recalculateDescTipPartitions(){
        descendantTipPartitions.clear();
        descendantTipPartitions(treeModel.getRoot(), descendantTipPartitions);
    }

    public void flagForDescendantRecalculation(TreeModel tree, NodeRef node){
        // todo Single traversals only? Watch this.
        flagForDescendantRecalculation(tree, node, updateNodeForSingleTraverse);
    }

    protected final void handleModelChangedEvent(Model model, Object object, int index) {
        super.handleModelChangedEvent(model, object, index);

        if(model==treeModel){
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

            traversalProbKnown = false;
            renormalisationNeeded = true;
        }

        if(model==branchMap){
            Arrays.fill(updateNodeForSingleTraverse, true);
            traversalProbKnown = false;
        }
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        super.handleVariableChangedEvent(variable, index, type);

        updateAllNodes(variable != infectionTimeBranchPositions && variable != infectiousTimePositions);

        if(variable == infectionTimeBranchPositions || variable == infectiousTimePositions){

            traversalProbKnown = false;
            renormalisationNeeded = false;
        }

    }

    protected void storeState() {
        super.storeState();
        core.storeState();
        storedNormTotalProb = normTotalProb;
        storedTotalLogProb = totalLogProb;
        storedDescendantTipPartitions = new HashMap<Integer, HashSet<AbstractCase>>(descendantTipPartitions);
        storedTraversalProbKnown = traversalProbKnown;
    }

    protected void restoreState() {
        super.restoreState();
        core.restoreState();
        normTotalProb = storedNormTotalProb;
        totalLogProb = storedTotalLogProb;
        descendantTipPartitions = storedDescendantTipPartitions;
        traversalProbKnown = storedTraversalProbKnown;
    }


    public final void makeDirty() {
        super.makeDirty();
        likelihoodKnown = false;
        traversalProbKnown = false;
        Arrays.fill(updateNode, true);
        Arrays.fill(updateNodeForSingleTraverse, true);
        renormalisationNeeded = true;
    }


    public class CaseToCaseLikelihoodCore extends GeneralLikelihoodCore {

        public CaseToCaseLikelihoodCore(){
            super(outbreak.size());
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

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        public static final String STARTING_NETWORK = "startingNetwork";
        public static final String INFECTION_TIMES = "infectionTimeBranchPositions";
        public static final String MAX_FIRST_INF_TO_ROOT = "maxFirstInfToRoot";

        public String getParserName() {
            return SIMPLE_CASE_TO_CASE;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            PartitionedTreeModel virusTree = (PartitionedTreeModel) xo.getChild(TreeModel.class);

            String startingNetworkFileName=null;

            if(xo.hasChildNamed(STARTING_NETWORK)){
                startingNetworkFileName = (String) xo.getElementFirstChild(STARTING_NETWORK);
            }

            AbstractOutbreak caseSet = (AbstractOutbreak) xo.getChild(AbstractOutbreak.class);

            CaseToCaseTreeLikelihood likelihood;

            Parameter infectionTimes = (Parameter) xo.getElementFirstChild(INFECTION_TIMES);

            Parameter earliestFirstInfection = (Parameter) xo.getElementFirstChild(MAX_FIRST_INF_TO_ROOT);

            likelihood = new SimpleCaseToCase(virusTree, caseSet, startingNetworkFileName, infectionTimes,
                    earliestFirstInfection);

            return likelihood;
        }

        public String getParserDescription() {
            return "This element provides a tree prior for a partitioned tree, with each partitioned tree generated" +
                    "by a coalescent process";
        }

        public Class getReturnType() {
            return WithinCaseCoalescent.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(PartitionedTreeModel.class, "The tree"),
                new ElementRule(WithinCaseCategoryOutbreak.class, "The set of outbreak"),
                new ElementRule("startingNetwork", String.class, "A CSV file containing a specified starting network",
                        true),
                new ElementRule(MAX_FIRST_INF_TO_ROOT, Parameter.class, "The maximum time from the first infection to" +
                        "the root node"),
                new ElementRule(INFECTION_TIMES, Parameter.class)
        };
    };

}
