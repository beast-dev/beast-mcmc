package dr.evomodel.epidemiology.casetocase;
import dr.app.tools.NexusExporter;
import dr.evolution.tree.*;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.epidemiology.casetocase.periodpriors.AbstractPeriodPriorDistribution;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.AbstractTreeLikelihood;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.MathUtils;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
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

public abstract class CaseToCaseTreeLikelihood extends AbstractTreeLikelihood implements Loggable, Citable,
        TreeTraitProvider {

    protected static final boolean DEBUG = false;

    /* The phylogenetic tree. */

    protected int noTips;
    protected int noCases;

    /* Mapping of outbreak to branches on the tree; old version is stored before operators are applied */

    protected BranchMapModel branchMap;

    /* Matches outbreak to external nodes */

    protected HashMap<AbstractCase, Integer> tipMap;
    private double estimatedLastSampleTime;
    protected TreeTraitProvider.Helper treeTraits = new Helper();

    /**
     * The set of cases
     */
    protected AbstractOutbreak outbreak;

    protected double[] infectionTimes;
    private double[] storedInfectionTimes;
    protected double[] infectiousPeriods;
    private double[] storedInfectiousPeriods;
    protected double[] infectiousTimes;
    private double[] storedInfectiousTimes;
    protected double[] latentPeriods;
    private double[] storedLatentPeriods;
    protected boolean[] recalculateCaseFlags;

    //because of the way the former works, we need a maximum value of the time from first infection to root node.

    protected Parameter maxFirstInfToRoot;

    // latent periods

    protected boolean hasLatentPeriods;


    // PUBLIC STUFF

    // Name

    public static final String CASE_TO_CASE_TREE_LIKELIHOOD = "caseToCaseTreeLikelihood";
    public static final String PARTITIONS_KEY = "partition";


    // Basic constructor.

    public CaseToCaseTreeLikelihood(PartitionedTreeModel virusTree, AbstractOutbreak caseData,
                                    Parameter maxFirstInfToRoot)
            throws TaxonList.MissingTaxonException {
        this(CASE_TO_CASE_TREE_LIKELIHOOD, virusTree, caseData, maxFirstInfToRoot);
    }

    // Constructor for an instance with a non-default name

    public CaseToCaseTreeLikelihood(String name, PartitionedTreeModel virusTree, AbstractOutbreak caseData,
                                    Parameter maxFirstInfToRoot) {
        super(name, caseData, virusTree);


        if(stateCount!=treeModel.getExternalNodeCount()){
            throw new RuntimeException("There are duplicate tip outbreak.");
        }

        noTips = virusTree.getExternalNodeCount();


        //subclasses should add outbreak as a model if it contains any information that ever changes

        outbreak = caseData;

        noCases = outbreak.getCases().size();

        addModel(outbreak);

        estimatedLastSampleTime = getLatestTaxonTime();

        //map outbreak to tips

        branchMap = virusTree.getBranchMap();

        addModel(branchMap);

        tipMap = new HashMap<AbstractCase, Integer>();

        //map the outbreak to the external nodes
        for(int i=0; i<virusTree.getExternalNodeCount(); i++){
            TreeModel.Node currentExternalNode = (TreeModel.Node)virusTree.getExternalNode(i);
            Taxon currentTaxon = currentExternalNode.taxon;
            for(AbstractCase thisCase : outbreak.getCases()){
                if(thisCase.wasEverInfected()) {
                    for (Taxon caseTaxon : thisCase.getAssociatedTaxa()) {
                        if (caseTaxon.equals(currentTaxon)) {
                            tipMap.put(thisCase, currentExternalNode.getNumber());
                        }
                    }
                }
            }
        }

        hasLatentPeriods = outbreak.hasLatentPeriods();

        infectionTimes = new double[outbreak.size()];
        infectiousPeriods = new double[outbreak.size()];

        if(hasLatentPeriods){
            infectiousTimes = new double[outbreak.size()];
            latentPeriods = new double[outbreak.size()];
        }


        recalculateCaseFlags = new boolean[outbreak.size()];
        Arrays.fill(recalculateCaseFlags, true);



        this.maxFirstInfToRoot = maxFirstInfToRoot;


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

        likelihoodKnown = false;
    }

    protected void prepareTree(String startingNetworkFileName){
        if(startingNetworkFileName==null){
            partitionAccordingToRandomTT(true);
        } else {
            partitionAccordingToSpecificTT(startingNetworkFileName);
        }

        prepareTimings();
        likelihoodKnown = false;
    }

    public AbstractOutbreak getOutbreak(){
        return outbreak;
    }

    public boolean hasLatentPeriods(){
        return hasLatentPeriods;
    }

    /* Get the date of the last tip */

    private double getLatestTaxonTime(){
        double latestTime = Double.NEGATIVE_INFINITY;
        for(AbstractCase thisCase : outbreak.getCases()){
            if (thisCase.wasEverInfected() && thisCase.getExamTime() > latestTime) {
                latestTime = thisCase.getExamTime();
            }
        }
        return latestTime;
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
        return tipLinked(node, branchMap);
    }

    private boolean tipLinked(NodeRef node, BranchMapModel map){
        NodeRef tip = treeModel.getNode(tipMap.get(map.get(node.getNumber())));
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


    //Counts the children of the current node which have the same painting as itself under the current map.
    //This will always be 1 if extended==false.



    public int countChildrenInSamePartition(NodeRef node, BranchMapModel map){
        if(treeModel.isExternal(node)){
            return -1;
        } else {
            int count = 0;
            AbstractCase parentCase = map.get(node.getNumber());
            for(int i=0; i< treeModel.getChildCount(node); i++){
                if(map.get(treeModel.getChild(node,i).getNumber())==parentCase){
                    count++;
                }
            }
            return count;
        }
    }

    public int countChildrenInSamePartition(NodeRef node){
        return countChildrenInSamePartition(node, branchMap);
    }


    public static NodeRef sibling(TreeModel tree, NodeRef node){
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

    // find all partitions of the descendant tips of the current node. If map is specified then it makes a map of node
    // number to possible partitions; map can be null.

    public HashSet<AbstractCase> descendantTipPartitions(NodeRef node, HashMap<Integer, HashSet<AbstractCase>> map){
        HashSet<AbstractCase> out = new HashSet<AbstractCase>();
        if(treeModel.isExternal(node)){
            out.add(branchMap.get(node.getNumber()));
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

    // change flags to indicate that something needs recalculation further down the tree

    protected static void flagForDescendantRecalculation(TreeModel tree, NodeRef node, boolean[] flags){
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
        flagForDescendantRecalculation(tree, node, updateNode);
    }

    //Return a set of nodes that are not descendants of (or equal to) the current node and are in the same partition as
    // it. If flagForRecalc is true, then this also sets the flags for likelihood recalculation for all these nodes
    // to true

    public HashSet<Integer> samePartitionDownTree(NodeRef node, boolean flagForRecalc){
        return samePartitionDownTree(node, branchMap, flagForRecalc);
    }

    public HashSet<Integer> samePartitionDownTree(NodeRef node, BranchMapModel map, boolean flagForRecalc){
        if(flagForRecalc){
            flagForDescendantRecalculation(treeModel, node);
        }
        HashSet<Integer> out = new HashSet<Integer>();
        AbstractCase painting = map.get(node.getNumber());
        NodeRef currentNode = node;
        NodeRef parentNode = treeModel.getParent(node);
        while(parentNode!=null && map.get(parentNode.getNumber())==painting){
            out.add(parentNode.getNumber());
            if(countChildrenInSamePartition(parentNode)==2){
                NodeRef otherChild = sibling(treeModel, currentNode);
                out.add(otherChild.getNumber());
                out.addAll(samePartitionUpTree(otherChild, map, flagForRecalc));
            }
            currentNode = parentNode;
            parentNode = treeModel.getParent(currentNode);
        }
        return out;
    }

    //Return a set of nodes that are descendants (and not equal to) the current node and are in the same partition as
    // it.

    public HashSet<Integer> samePartitionUpTree(NodeRef node, boolean flagForRecalc){
        return samePartitionUpTree(node, branchMap, flagForRecalc);
    }

    public HashSet<Integer> samePartitionUpTree(NodeRef node, BranchMapModel map, boolean flagForRecalc){
        HashSet<Integer> out = new HashSet<Integer>();
        AbstractCase painting = map.get(node.getNumber());
        boolean creepsFurther = false;
        for(int i=0; i< treeModel.getChildCount(node); i++){
            if(map.get(treeModel.getChild(node,i).getNumber())==painting){
                creepsFurther = true;
                out.add(treeModel.getChild(node,i).getNumber());
                out.addAll(samePartitionUpTree(treeModel.getChild(node, i), map, flagForRecalc));
            }
        }
        if(flagForRecalc && !creepsFurther){
            flagForDescendantRecalculation(treeModel, node);
        }
        return out;
    }

    public Integer[] samePartition(NodeRef node, boolean flagForRecalc){
        return samePartition(node, branchMap, flagForRecalc);
    }

    private Integer[] samePartition(NodeRef node, BranchMapModel map, boolean flagForRecalc){
        HashSet<Integer> out = new HashSet<Integer>();
        out.add(node.getNumber());
        out.addAll(samePartitionDownTree(node, map, flagForRecalc));
        out.addAll(samePartitionUpTree(node, map, flagForRecalc));
        return out.toArray(new Integer[out.size()]);
    }

    // returns all nodes that are the earliest nodes in the partitions corresponding to this outbreak' children in
    // the _transmission_ tree.

    private Integer[] getAllChildInfectionNodes(AbstractCase thisCase){
        HashSet<Integer> out = new HashSet<Integer>();
        NodeRef tip = treeModel.getNode(tipMap.get(thisCase));
        Integer[] partition = samePartition(tip, false);
        for (Integer i : partition) {
            NodeRef node = treeModel.getNode(i);
            if (!treeModel.isExternal(node)) {
                for (int j = 0; j < treeModel.getChildCount(node); j++) {
                    NodeRef child = treeModel.getChild(node, j);
                    if(branchMap.get(child.getNumber())!=thisCase){
                        out.add(child.getNumber());
                    }
                }
            }
        }
        return out.toArray(new Integer[out.size()]);
    }

    public NodeRef getTip(AbstractCase thisCase){
        return treeModel.getNode(tipMap.get(thisCase));
    }

    public NodeRef getEarliestNodeInPartition(AbstractCase thisCase, BranchMapModel branchMap){
        if(thisCase.wasEverInfected()) {
            NodeRef child = treeModel.getNode(tipMap.get(thisCase));
            NodeRef parent = treeModel.getParent(child);
            boolean transmissionFound = false;
            while (!transmissionFound) {
                if (branchMap.get(child.getNumber()) != branchMap.get(parent.getNumber())) {
                    transmissionFound = true;
                } else {
                    child = parent;
                    parent = treeModel.getParent(child);
                    if (parent == null) {
                        transmissionFound = true;
                    }
                }
            }
            return child;
        }
        return null;
    }

    public NodeRef getEarliestNodeInPartition(AbstractCase thisCase){
        return getEarliestNodeInPartition(thisCase, branchMap);
    }

    public HashSet<AbstractCase> getDescendants(AbstractCase thisCase){
        HashSet<AbstractCase> out = new HashSet<AbstractCase>(getInfectees(thisCase));

        if(thisCase.wasEverInfected()) {
            for (AbstractCase child : out) {
                out.addAll(getDescendants(child));
            }
        }
        return out;
    }


    public Integer[] getParentsArray(){
        Integer[] out = new Integer[outbreak.size()];
        for(AbstractCase thisCase : outbreak.getCases()){
            if(thisCase.wasEverInfected()) {
                out[outbreak.getCaseIndex(thisCase)] = outbreak.getCaseIndex(getInfector(thisCase));
            } else {
                out[outbreak.getCaseIndex(thisCase)] = null;
            }
        }
        return out;
    }


    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************


    protected void handleModelChangedEvent(Model model, Object object, int index) {

        if(!(model instanceof AbstractPeriodPriorDistribution)){

            if (model == treeModel) {

                if(object instanceof PartitionedTreeModel.PartitionsChangedEvent){
                    HashSet<AbstractCase> changedPartitions =
                            ((PartitionedTreeModel.PartitionsChangedEvent)object).getCasesToRecalculate();
                    for(AbstractCase aCase : changedPartitions){
                        recalculateCase(aCase);

                    }
                }
            } else if (model == branchMap){
                if(object instanceof ArrayList){

                    for(int i=0; i<((ArrayList) object).size(); i++){
                        BranchMapModel.BranchMapChangedEvent event
                                =  (BranchMapModel.BranchMapChangedEvent)((ArrayList) object).get(i);

                        recalculateCase(event.getOldCase());
                        recalculateCase(event.getNewCase());

                        NodeRef node = treeModel.getNode(event.getNodeToRecalculate());
                        NodeRef parent = treeModel.getParent(node);

                        if(parent!=null){
                            recalculateCase(branchMap.get(parent.getNumber()));
                        }
                    }
                } else {
                    throw new RuntimeException("Unanticipated model changed event from BranchMapModel");
                }
            } else if (model == outbreak){

                if(object instanceof AbstractCase){
                    recalculateCase((AbstractCase)object);
                } else {
                    for (AbstractCase aCase : outbreak.getCases()) {
                        recalculateCase(aCase);
                    }
                }
            }

            fireModelChanged(model);

            likelihoodKnown = false;
        }
    }

    protected void recalculateCase(int index){
        recalculateCaseFlags[index] = true;
    }

    protected void recalculateCase(AbstractCase aCase){
        if(aCase.wasEverInfected()) {
            recalculateCase(outbreak.getCaseIndex(aCase));
        }
    }


    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************


    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

        fireModelChanged();

        likelihoodKnown = false;
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    /**
     * Stores the precalculated state (in this case the node labels and subtree likelihoods)
     */

    protected void storeState() {
        super.storeState();
        storedInfectionTimes = Arrays.copyOf(infectionTimes, infectionTimes.length);
        storedInfectiousPeriods = Arrays.copyOf(infectiousPeriods, infectiousPeriods.length);
        if(hasLatentPeriods){
            storedInfectiousTimes = Arrays.copyOf(infectiousTimes, infectionTimes.length);
            storedLatentPeriods = Arrays.copyOf(latentPeriods, latentPeriods.length);
        }
    }

    /**
     * Restores the precalculated state.
     */

    protected void restoreState() {
        super.restoreState();
        infectionTimes = storedInfectionTimes;
        infectiousPeriods = storedInfectiousPeriods;
        if(hasLatentPeriods){
            infectiousTimes = storedInfectiousTimes;
            latentPeriods = storedLatentPeriods;
        }
    }

    protected final void acceptState() {
    } // nothing to do

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    public final BranchMapModel getBranchMap(){
        return branchMap;
    }

    public final PartitionedTreeModel getTreeModel(){
        return (PartitionedTreeModel)treeModel;
    }


    public void makeDirty() {
        likelihoodKnown = false;
        Arrays.fill(recalculateCaseFlags, true);
    }


    protected void prepareTimings(){

        infectionTimes = getInfectionTimes(true);

        if(hasLatentPeriods){
            infectiousTimes = getInfectiousTimes(true);
        }

        infectiousPeriods = getInfectiousPeriods(true);

        if(hasLatentPeriods){
            latentPeriods = getLatentPeriods(true);
        }

        Arrays.fill(recalculateCaseFlags, false);
    }

    /**
     * Calculates the log likelihood of this set of node labels given the tree.
     */

    protected abstract double calculateLogLikelihood();



    // if no infectious models, just need to check whether any infections occur after the infector was no
    // longer infectious

    protected boolean isAllowed(){
        return isAllowed(treeModel.getRoot());
    }

    private boolean isAllowed(NodeRef node){
        if(!treeModel.isRoot(node)){
            AbstractCase childCase = branchMap.get(node.getNumber());
            AbstractCase parentCase = branchMap.get(treeModel.getParent(node).getNumber());
            if(childCase!=parentCase){
                double infectionTime = infectionTimes[outbreak.getCaseIndex(childCase)];
                if(infectionTime>parentCase.getCullTime()
                        || (hasLatentPeriods && infectionTime<infectiousTimes[outbreak.getCaseIndex(parentCase)])){
                    return false;
                }
            }
        }
        return treeModel.isExternal(node) ||
                (isAllowed(treeModel.getChild(node, 0)) && isAllowed(treeModel.getChild(node, 1)));
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

    public AbstractCase getInfector(AbstractCase thisCase, BranchMapModel branchMap){
        if(thisCase.wasEverInfected()) {
            NodeRef tip = treeModel.getNode(tipMap.get(thisCase));
            return getInfector(tip, branchMap);
        }
        return null;
    }

    public AbstractCase getRootCase(){
        return branchMap.get(treeModel.getRoot().getNumber());
    }

    public HashSet<AbstractCase> getInfectees(AbstractCase thisCase){
        return getInfectees(thisCase, branchMap);
    }

    public HashSet<AbstractCase> getInfectees(AbstractCase thisCase, BranchMapModel branchMap){
        if(thisCase.wasEverInfected()) {
            return getInfecteesInClade(getEarliestNodeInPartition(thisCase), branchMap);
        }
        return new HashSet<AbstractCase>();
    }

    public HashSet<AbstractCase> getInfecteesInClade(NodeRef node, BranchMapModel branchMap){
        HashSet<AbstractCase> out = new HashSet<AbstractCase>();
        if(treeModel.isExternal(node)){
            return out;
        } else {
            AbstractCase thisCase = branchMap.get(node.getNumber());
            for(int i=0; i<treeModel.getChildCount(node); i++){
                NodeRef child = treeModel.getChild(node, i);
                AbstractCase childCase = branchMap.get(child.getNumber());
                if(childCase!=thisCase){
                    out.add(childCase);
                } else {
                    out.addAll(getInfecteesInClade(child, branchMap));
                }
            }
            return out;
        }
    }

    public double getInfectionTime(AbstractCase thisCase){

        if(!recalculateCaseFlags[outbreak.getCaseIndex(thisCase)]){
            return infectionTimes[outbreak.getCaseIndex(thisCase)];
        } else {
            if(thisCase.wasEverInfected()) {
                NodeRef child = getEarliestNodeInPartition(thisCase);
                NodeRef parent = treeModel.getParent(child);

                if (parent != null) {

                    double min = heightToTime(treeModel.getNodeHeight(parent));

                    // Let the likelihood evaluate to zero due to culling dates if it must...

                    double max = heightToTime(treeModel.getNodeHeight(child));


                    return getInfectionTime(min, max, thisCase);
                } else {
                    return getRootInfectionTime(branchMap);
                }
            } else {
                return Double.POSITIVE_INFINITY;
            }
        }
    }

    private double getInfectionTime(double min, double max, AbstractCase infected){
        final double branchLength = max-min;
        return max - branchLength*infected.getInfectionBranchPosition().getParameterValue(0);
    }


    public double[] getInfectionTimes(boolean recalculate){
        if(recalculate) {
            for(int i=0; i<noCases; i++){
                if(recalculateCaseFlags[i]){
                    infectionTimes[i] = getInfectionTime(outbreak.getCase(i));
                }
            }
        }
        return infectionTimes;
    }

    public void setInfectionTime(AbstractCase thisCase, double time){

        setInfectionHeight(thisCase, timeToHeight(time));

    }

    public void setInfectionHeight(AbstractCase thisCase, double height){
        if(thisCase.wasEverInfected()) {
            NodeRef child = getEarliestNodeInPartition(thisCase);
            NodeRef parent = treeModel.getParent(child);

            double minHeight = treeModel.getNodeHeight(child);
            double maxHeight = parent != null ? treeModel.getNodeHeight(parent)
                    : minHeight + maxFirstInfToRoot.getParameterValue(0);

            if (height < minHeight || height > maxHeight) {
                throw new RuntimeException("Trying to set an infection time outside the branch on which it must occur");
            }

            double branchPosition = (height - minHeight) / (maxHeight - minHeight);

            thisCase.setInfectionBranchPosition(branchPosition);

        }

    }

    public double getInfectiousTime(AbstractCase thisCase){
        if(!hasLatentPeriods){
            return getInfectionTime(thisCase);
        } else {
            if (recalculateCaseFlags[outbreak.getCaseIndex(thisCase)]) {
                if(thisCase.wasEverInfected()) {

                    String latentCategory = ((CategoryOutbreak) outbreak).getLatentCategory(thisCase);
                    Parameter latentPeriod = ((CategoryOutbreak) outbreak).getLatentPeriod(latentCategory);
                    infectiousTimes[outbreak.getCaseIndex(thisCase)] = getInfectionTime(thisCase)
                            + latentPeriod.getParameterValue(0);
                } else {
                    infectiousTimes[outbreak.getCaseIndex(thisCase)] = Double.POSITIVE_INFINITY;
                }
            }
            return infectiousTimes[outbreak.getCaseIndex(thisCase)];
        }

    }

    public double[] getInfectiousTimes(boolean recalculate){
        if(recalculate){
            for(int i=0; i<noCases; i++){
                if(recalculateCaseFlags[i]){
                    infectiousTimes[i] = getInfectiousTime(outbreak.getCase(i));
                }
            }
        }
        return infectiousTimes;
    }


    public double getInfectiousPeriod(AbstractCase thisCase){
        if(recalculateCaseFlags[outbreak.getCaseIndex(thisCase)]){
            if(thisCase.wasEverInfected()) {

                if (!hasLatentPeriods) {
                    double infectionTime = getInfectionTime(thisCase);
                    double cullTime = thisCase.getCullTime();
                    infectiousPeriods[outbreak.getCaseIndex(thisCase)] = cullTime - infectionTime;
                } else {
                    double infectiousTime = getInfectiousTime(thisCase);
                    double cullTime = thisCase.getCullTime();
                    infectiousPeriods[outbreak.getCaseIndex(thisCase)] = cullTime - infectiousTime;
                }
            } else {
                infectiousPeriods[outbreak.getCaseIndex(thisCase)] = 0;
            }
        }
        return infectiousPeriods[outbreak.getCaseIndex(thisCase)];
    }

    public double[] getInfectiousPeriods(boolean recalculate){
        if(recalculate){
            for(int i=0; i<noCases; i++){
                if(recalculateCaseFlags[i]){
                    infectiousPeriods[i] = getInfectiousPeriod(outbreak.getCase(i));
                }
            }
        }
        return infectiousPeriods;
    }

    public Double[] getNonzeroInfectiousPeriods(){
        ArrayList<Double> out = new ArrayList<Double>();


        for(int i=0; i<noCases; i++){
            AbstractCase thisCase = outbreak.getCase(i);

            if(thisCase.wasEverInfected()){
                out.add(getInfectiousPeriod(thisCase));
            }
        }

        return out.toArray(new Double[out.size()]);
    }


    public double getLatentPeriod(AbstractCase thisCase){
        if(!hasLatentPeriods || !thisCase.wasEverInfected()){
            return 0.0;
        }
        if(recalculateCaseFlags[outbreak.getCaseIndex(thisCase)]){
            latentPeriods[outbreak.getCaseIndex(thisCase)] = getInfectiousTime(thisCase) - getInfectionTime(thisCase);
        }
        return latentPeriods[outbreak.getCaseIndex(thisCase)];
    }

    public double[] getLatentPeriods(boolean recalculate){
        if(recalculate){
            for(int i=0; i<noCases; i++){
                if(recalculateCaseFlags[i]){
                    latentPeriods[i] = getLatentPeriod(outbreak.getCase(i));
                }
            }
        }
        return latentPeriods;
    }

    public Double[] getNonzeroLatentPeriods(){
        ArrayList<Double> out = new ArrayList<Double>();


        for(int i=0; i<noCases; i++){
            AbstractCase thisCase = outbreak.getCase(i);

            if(thisCase.wasEverInfected()){
                out.add(getLatentPeriod(thisCase));
            }
        }

        return out.toArray(new Double[out.size()]);
    }



    public double[] getInfectedPeriods(boolean recalculate){
        if(!hasLatentPeriods){
            return getInfectiousPeriods(recalculate);
        } else {
            double[] out = new double[noCases];
            for(int i=0; i<noCases; i++){
                out[i] = getInfectedPeriod(outbreak.getCase(i));
            }
            return out;
        }
    }


    public Double[] getNonzeroInfectedPeriods(){
        ArrayList<Double> out = new ArrayList<Double>();


        for(int i=0; i<noCases; i++){
            AbstractCase thisCase = outbreak.getCase(i);

            if(thisCase.wasEverInfected()){
                out.add(getInfectedPeriod(thisCase));
            }
        }

        return out.toArray(new Double[out.size()]);
    }




    public double getInfectedPeriod(AbstractCase thisCase){
        if(thisCase.wasEverInfected) {
            return thisCase.getCullTime() - getInfectionTime(thisCase);
        }
        return 0;
    }


    // return an array of the mean, median, variance and standard deviation of the given array
    // @todo this is pretty wasteful since it gets called so many times per log entry

    public static Double[] getSummaryStatistics(Double[] variable){

        double[] primitiveVariable = new double[variable.length];
        for(int i=0; i<variable.length; i++){
            primitiveVariable[i] = variable[i];
        }

        DescriptiveStatistics stats = new DescriptiveStatistics(primitiveVariable);
        Double[] out = new Double[4];
        out[0] = stats.getMean();
        out[1] = stats.getPercentile(50);
        out[2] = stats.getVariance();
        out[3] = stats.getStandardDeviation();
        return out;
    }



    private double getRootInfectionTime(BranchMapModel branchMap){
        NodeRef root = treeModel.getRoot();
        AbstractCase rootCase = branchMap.get(root.getNumber());
        final double branchLength = maxFirstInfToRoot.getParameterValue(0);

        return heightToTime(treeModel.getNodeHeight(root)
                + branchLength * rootCase.getInfectionBranchPosition().getParameterValue(0));

    }

    protected double getRootInfectionTime(){
        AbstractCase rootCase = branchMap.get(treeModel.getRoot().getNumber());
        return getInfectionTime(rootCase);
    }

    public AbstractCase getInfector(NodeRef node, BranchMapModel branchMap){
        if(treeModel.isRoot(node) || node.getNumber() == treeModel.getRoot().getNumber()){
            return null;
        } else {
            AbstractCase nodeCase = branchMap.get(node.getNumber());
            if(branchMap.get(treeModel.getParent(node).getNumber())!=nodeCase){
                return branchMap.get(treeModel.getParent(node).getNumber());
            } else {
                return getInfector(treeModel.getParent(node), branchMap);
            }
        }
    }

    public boolean checkPartitions(){
        return checkPartitions(branchMap, true);
    }

    protected boolean checkPartitions(BranchMapModel map, boolean verbose){
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
            debugOutputTree(map, "checkPartitionProblem", false);
            throw new BadPartitionException("Tree is not partitioned properly");
        }
        return !foundProblem;
    }


    /* Return the partition of the parent of this node */

    public AbstractCase getParentCase(NodeRef node){
        return branchMap.get(treeModel.getParent(node).getNumber());
    }

    /* Populates the branch map for external nodes */

    private AbstractCase[] prepareExternalNodeMap(AbstractCase[] map){
        for(int i=0; i< treeModel.getExternalNodeCount(); i++){
            TreeModel.Node currentExternalNode = (TreeModel.Node) treeModel.getExternalNode(i);
            Taxon currentTaxon = currentExternalNode.taxon;
            for(AbstractCase thisCase : outbreak.getCases()){
                if(thisCase.wasEverInfected()) {
                    for (Taxon caseTaxon : thisCase.getAssociatedTaxa()) {
                        if (caseTaxon.equals(currentTaxon)) {
                            map[currentExternalNode.getNumber()] = thisCase;
                        }
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
                    specificParentMap.put(outbreak.getCase(splitLine[0]), outbreak.getCase(splitLine[1]));
                } else {
                    specificParentMap.put(outbreak.getCase(splitLine[0]), null);
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
        branchMap.setAll(prepareExternalNodeMap(new AbstractCase[treeModel.getNodeCount()]), true);

        AbstractCase firstCase=null;
        for(AbstractCase aCase : outbreak.getCases()){
            if(aCase.wasEverInfected() && map.get(aCase)==null){
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
        branchMap.set(node.getNumber(), thisCase, true);
        if(tipLinked(node)){
            for(int i=0; i<treeModel.getChildCount(node); i++){
                specificallyPartitionUpwards(treeModel.getChild(node, i), thisCase, map);
            }
        } else {
            branchMap.set(node.getNumber(), null, true);
            HashSet<AbstractCase> children = new HashSet<AbstractCase>();
            for(AbstractCase aCase : outbreak.getCases()){
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
                branchMap.set(node.getNumber(), child, true);
            } else {
                branchMap.set(node.getNumber(), thisCase, true);
            }
            for(int i=0; i<treeModel.getChildCount(node); i++){
                specificallyPartitionUpwards(treeModel.getChild(node, i), branchMap.get(node.getNumber()), map);
            }
        }

    }


    /* Assigns a phylogenetic tree node and its children to a partition according to a specified map of child to parent
    outbreak. This only works on the non-extended version right now, watch it. */

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

    private void partitionAccordingToRandomTT(boolean checkNonZero){
        boolean gotOne = false;
        int tries = 1;
        System.out.println("Generating a random starting partition of the tree (checking nonzero likelihood for all " +
                "branches and repeating up to 100 times until a start with nonzero likelihood is found)");
        System.out.print("Attempt: ");
        while(!gotOne){

            likelihoodKnown = false;
            boolean failed = false;
            System.out.print(tries + "...");
            branchMap.setAll(prepareExternalNodeMap(new AbstractCase[treeModel.getNodeCount()]), true);
            //Warning - if the BadPartitionException in randomlyAssignNode might be caused by a bug rather than both
            //likelihoods rounding to zero, you want to stop catching this to investigate.

            try{
                partitionAccordingToRandomTT(branchMap, checkNonZero);
            } catch(BadPartitionException e){
                failed = true;
            }

            makeDirty();

            infectionTimes = getInfectionTimes(true);
            if(hasLatentPeriods){
                infectiousTimes = getInfectiousTimes(true);
            }

            infectiousPeriods = getInfectiousPeriods(true);
            if(hasLatentPeriods){
                latentPeriods = getLatentPeriods(true);
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



    /* Partitions a phylogenetic tree randomly; if checkNonZero is true, make sure all branch likelihoods are nonzero
    in the process (this sometimes still results in a zero likelihood for the whole tree, but is much less likely to).
    */

    private BranchMapModel partitionAccordingToRandomTT(BranchMapModel map, boolean checkNonZero){
        makeDirty();
        TreeModel.Node root = (TreeModel.Node) treeModel.getRoot();
        randomlyAssignNode(root, map, checkNonZero);

        return map;
    }

    private AbstractCase randomlyAssignNode(TreeModel.Node node, BranchMapModel map, boolean checkNonZero){
        //this makes a non-extended partition. This is OK, but if it keeps giving zero likelihoods then you could do
        //something else

        if(node.isExternal()){
            return map.get(node.getNumber());
        } else {

            AbstractCase[] choices = new AbstractCase[2];
            for(int i=0; i<node.getChildCount(); i++){
                if((map.get(node.getChild(i).getNumber())==null)){
                    choices[i] = randomlyAssignNode(node.getChild(i), map, checkNonZero);
                } else {
                    choices[i] = map.get(node.getChild(i).getNumber());
                }
            }
            int randomSelection = MathUtils.nextInt(2);
            int decision;
            if(checkNonZero){
                Boolean[] branchLogLs = new Boolean[2];
                for(int i=0; i<2; i++){
                    double nodeTime = getNodeTime(node);
                    double branchLength = getNodeTime(treeModel.getChild(node, 1-i)) - getNodeTime(node);
                    AbstractCase infector = choices[i];
                    AbstractCase infectee = choices[1-i];

                    branchLogLs[i] = !infector.culledYet(nodeTime
                            + infectee.getInfectionBranchPosition().getParameterValue(0)*branchLength);
                }
                if(!branchLogLs[0] && !branchLogLs[1]){
                    throw new BadPartitionException("Both branch possibilities have zero likelihood: "
                            +node.toString()+", outbreak " + choices[0].getName() + " and " + choices[1].getName() + ".");
                } else if(!branchLogLs[0] || !branchLogLs[1]){
                    if(!branchLogLs[0]){
                        decision = 1;
                    } else {
                        decision = 0;
                    }
                } else {
                    decision = randomSelection;
                }
            } else {
                decision = randomSelection;
            }
            AbstractCase winner = choices[decision];
            map.getArray()[node.getNumber()]=winner;
            return winner;
        }
    }

    public void debugOutputTree(String fileName, boolean rewire){
        debugOutputTree(branchMap, fileName, rewire);
    }


    public void debugOutputTree(BranchMapModel map, String fileName, boolean rewire){
        try{
            FlexibleTree treeCopy;
            if(!rewire){
                treeCopy = new FlexibleTree(treeModel);
                for(int j=0; j<treeCopy.getNodeCount(); j++){
                    FlexibleNode node = (FlexibleNode)treeCopy.getNode(j);
                    node.setAttribute("Number", node.getNumber());
                    node.setAttribute("Time", heightToTime(node.getHeight()));
                    node.setAttribute(PARTITIONS_KEY, map.get(node.getNumber()));
                }
            } else {
                treeCopy = rewireTree(treeModel);
            }
            NexusExporter testTreesOut = new NexusExporter(new PrintStream(fileName));
            testTreesOut.exportTree(treeCopy);
        } catch (IOException ignored) {System.out.println("IOException");}
    }

    public FlexibleTree rewireTree(Tree tree){
        prepareTimings();

        FlexibleTree outTree = new FlexibleTree(tree, true);

        for(int j=0; j<outTree.getNodeCount(); j++){
            FlexibleNode node = (FlexibleNode)outTree.getNode(j);
            node.setAttribute("Number", node.getNumber());
            node.setAttribute("Time", heightToTime(node.getHeight()));
            node.setAttribute(PARTITIONS_KEY, branchMap.get(node.getNumber()));
        }

        for(AbstractCase aCase : outbreak.getCases()){
            NodeRef originalNode = getEarliestNodeInPartition(aCase);
            int infectionNodeNo = originalNode.getNumber();
            if(!treeModel.isRoot(originalNode)){
                NodeRef originalParent = treeModel.getParent(originalNode);
                double nodeTime = getNodeTime(originalNode);
                double infectionTime = getInfectionTime(aCase);
                double heightToBreakBranch = getHeight(originalNode) +  (nodeTime - infectionTime);
                FlexibleNode newNode = (FlexibleNode)outTree.getNode(infectionNodeNo);
                FlexibleNode oldParent = (FlexibleNode)outTree.getParent(newNode);

                outTree.beginTreeEdit();
                outTree.removeChild(oldParent, newNode);
                FlexibleNode infectionNode = new FlexibleNode();
                infectionNode.setHeight(heightToBreakBranch);
                infectionNode.setLength(oldParent.getHeight() - heightToBreakBranch);
                infectionNode.setAttribute(PARTITIONS_KEY, getNodePartition(treeModel, originalParent));
                newNode.setLength(nodeTime - infectionTime);

                outTree.addChild(oldParent, infectionNode);
                outTree.addChild(infectionNode, newNode);
                outTree.endTreeEdit();
            } else {
                double nodeTime = getNodeTime(originalNode);
                double infectionTime = getInfectionTime(aCase);
                double heightToInstallRoot = getHeight(originalNode) +  (nodeTime - infectionTime);
                FlexibleNode newNode = (FlexibleNode)outTree.getNode(infectionNodeNo);
                outTree.beginTreeEdit();
                FlexibleNode infectionNode = new FlexibleNode();
                infectionNode.setHeight(heightToInstallRoot);
                infectionNode.setAttribute(PARTITIONS_KEY, "The_Ether");
                outTree.addChild(infectionNode, newNode);
                newNode.setLength(heightToInstallRoot - getHeight(originalNode));
                outTree.setRoot(infectionNode);
                outTree.endTreeEdit();
            }
        }

        outTree = new FlexibleTree((FlexibleNode)outTree.getRoot());


        for(int i=0; i<outTree.getNodeCount(); i++){
            NodeRef node = outTree.getNode(i);
            NodeRef parent = outTree.getParent(node);
            if(parent!=null && outTree.getNodeHeight(node)>outTree.getNodeHeight(parent)){
                try{
                    NexusExporter exporter = new NexusExporter(new PrintStream("fancyProblem.nex"));
                    exporter.exportTree(outTree);
                } catch(IOException e){
                    e.printStackTrace();
                }
                try{
                    checkPartitions();
                } catch(BadPartitionException e){
                    System.out.print("Rewiring messed up because of partition problem.");
                }


                throw new RuntimeException("Rewiring messed up; investigate");
            }

        }


        return outTree;
    }

    //************************************************************************
    // Loggable implementation
    //************************************************************************


    public LogColumn[] getColumns(){
        LogColumn[] columns = new LogColumn[outbreak.infectedSize()];
        for(int i=0; i< outbreak.size(); i++){
            final AbstractCase infected = outbreak.getCase(i);
            if(infected.wasEverInfected()) {
                columns[i] = new LogColumn.Abstract(infected.toString() + "_infector") {
                    protected String getFormattedValue() {
                        if (getInfector(infected) == null) {
                            return "Start";
                        } else {
                            return getInfector(infected).toString();
                        }
                    }
                };
            }
        }
        return columns;
    }

    public LogColumn[] passColumns(){
        ArrayList<LogColumn> columns = new ArrayList<LogColumn>();
        for(int i=0; i< outbreak.size(); i++){
            final AbstractCase infected = outbreak.getCase(i);
            if(infected.wasEverInfected()) {
                columns.add(new LogColumn.Abstract(infected.toString() + "_infection_date") {
                    protected String getFormattedValue() {
                        return String.valueOf(getInfectionTime(infected));
                    }
                });
            }
        }
        if(hasLatentPeriods){
            for(int i=0; i< outbreak.size(); i++){
                final AbstractCase infected = outbreak.getCase(i);
                if(infected.wasEverInfected()) {
                    columns.add(new LogColumn.Abstract(infected.toString() + "_infectious_date") {
                        protected String getFormattedValue() {
                            return String.valueOf(getInfectiousTime(infected));
                        }
                    });
                }
            }
            for(int i=0; i< outbreak.size(); i++){
                final AbstractCase infected = outbreak.getCase(i);
                if(infected.wasEverInfected()) {
                    columns.add(new LogColumn.Abstract(infected.toString() + "_latent_period") {
                        protected String getFormattedValue() {
                            return String.valueOf(getLatentPeriod(infected));
                        }
                    });
                }
            }
        }
        for(int i=0; i< outbreak.size(); i++){
            final AbstractCase infected = outbreak.getCase(i);
            if(infected.wasEverInfected()) {
                columns.add(new LogColumn.Abstract(infected.toString() + "_infectious_period") {
                    protected String getFormattedValue() {
                        return String.valueOf(getInfectiousPeriod(infected));
                    }
                });
            }
        }
        if(hasLatentPeriods){
            for(int i=0; i< outbreak.size(); i++){
                final AbstractCase infected = outbreak.getCase(i);
                if(infected.wasEverInfected()) {
                    columns.add(new LogColumn.Abstract(infected.toString() + "_infected_period") {
                        protected String getFormattedValue() {
                            return String.valueOf(
                                    getInfectiousPeriod(infected) + getLatentPeriod(infected));
                        }
                    });
                }
            }
        }

        columns.add(new LogColumn.Abstract("infectious_period.mean"){
            protected String getFormattedValue() {
                return String.valueOf(CaseToCaseTreeLikelihood
                        .getSummaryStatistics(getNonzeroInfectiousPeriods())[0]);
            }
        });
        columns.add(new LogColumn.Abstract("infectious_period.median"){
            protected String getFormattedValue() {
                return String.valueOf(CaseToCaseTreeLikelihood
                        .getSummaryStatistics(getNonzeroInfectiousPeriods())[1]);
            }
        });
        columns.add(new LogColumn.Abstract("infectious_period.var") {
            protected String getFormattedValue() {
                return String.valueOf(CaseToCaseTreeLikelihood
                        .getSummaryStatistics(getNonzeroInfectiousPeriods())[2]);
            }
        });
        columns.add(new LogColumn.Abstract("infectious_period.stdev"){
            protected String getFormattedValue() {
                return String.valueOf(CaseToCaseTreeLikelihood
                        .getSummaryStatistics(getNonzeroInfectiousPeriods())[3]);
            }
        });
        if(hasLatentPeriods){
            columns.add(new LogColumn.Abstract("latent_period.mean"){
                protected String getFormattedValue() {
                    return String.valueOf(CaseToCaseTreeLikelihood
                            .getSummaryStatistics(getNonzeroLatentPeriods())[0]);
                }
            });
            columns.add(new LogColumn.Abstract("latent_period.median"){
                protected String getFormattedValue() {
                    return String.valueOf(CaseToCaseTreeLikelihood
                            .getSummaryStatistics(getNonzeroLatentPeriods())[1]);
                }
            });
            columns.add(new LogColumn.Abstract("latent_period.var") {
                protected String getFormattedValue() {
                    return String.valueOf(CaseToCaseTreeLikelihood
                            .getSummaryStatistics(getNonzeroLatentPeriods())[2]);
                }
            });
            columns.add(new LogColumn.Abstract("latent_period.stdev"){
                protected String getFormattedValue() {
                    return String.valueOf(CaseToCaseTreeLikelihood
                            .getSummaryStatistics(getNonzeroLatentPeriods())[3]);
                }
            });
            columns.add(new LogColumn.Abstract("infected_period.mean"){
                protected String getFormattedValue() {
                    return String.valueOf(CaseToCaseTreeLikelihood
                            .getSummaryStatistics(getNonzeroInfectedPeriods())[0]);
                }
            });
            columns.add(new LogColumn.Abstract("infected_period.median"){
                protected String getFormattedValue() {
                    return String.valueOf(CaseToCaseTreeLikelihood
                            .getSummaryStatistics(getNonzeroInfectedPeriods())[1]);
                }
            });
            columns.add(new LogColumn.Abstract("infected_period.var") {
                protected String getFormattedValue() {
                    return String.valueOf(CaseToCaseTreeLikelihood
                            .getSummaryStatistics(getNonzeroInfectedPeriods())[2]);
                }
            });
            columns.add(new LogColumn.Abstract("infected_period.stdev"){
                protected String getFormattedValue() {
                    return String.valueOf(CaseToCaseTreeLikelihood
                            .getSummaryStatistics(getNonzeroInfectedPeriods())[3]);
                }
            });
            for(int i=0; i< outbreak.size(); i++){
                final AbstractCase infected = outbreak.getCase(i);
                if(infected.wasEverInfected()) {
                    columns.add(new LogColumn.Abstract(infected.toString() + "_ibp") {
                        protected String getFormattedValue() {
                            return String.valueOf(infected.getInfectionBranchPosition().getParameterValue(0));
                        }
                    });
                }
            }
        }

        return columns.toArray(new LogColumn[columns.size()]);


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
            // we're trying to annotate a partitioned tree, we hope
            try{
                NodeRef oldNode = treeModel.getNode((Integer)tree.getNodeAttribute(node,"Number"));
                if(treeModel.getNodeHeight(oldNode)!=tree.getNodeHeight(node)){
                    throw new RuntimeException("Can only reconstruct states on treeModel given to constructor or a " +
                            "partitioned tree derived from it");
                } else {
                    return branchMap.get(oldNode.getNumber()).toString();
                }
            } catch(NullPointerException e){
                if(tree.isRoot(node)){
                    return "Start";
                } else {
                    NodeRef parent = tree.getParent(node);
                    int originalParentNumber = (Integer)tree.getNodeAttribute(parent,"Number");
                    return branchMap.get(originalParentNumber).toString();
                }
            }
        } else {
            if (!likelihoodKnown) {
                calculateLogLikelihood();
                likelihoodKnown = true;
            }
            return branchMap.get(node.getNumber()).toString();
        }
    }

}








