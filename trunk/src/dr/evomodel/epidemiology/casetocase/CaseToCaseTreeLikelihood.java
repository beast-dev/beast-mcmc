package dr.evomodel.epidemiology.casetocase;
import dr.app.tools.NexusExporter;
import dr.evolution.tree.*;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.AbstractTreeLikelihood;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.MathUtils;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;

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

    protected final static boolean DEBUG = false;

    /* The phylogenetic tree. */

    protected int noTips;

    /* Mapping of cases to branches on the tree; old version is stored before operators are applied */

    protected BranchMapModel branchMap;

    /* Matches cases to external nodes */

    protected HashMap<AbstractCase, Integer> tipMap;
    private double estimatedLastSampleTime;
    protected TreeTraitProvider.Helper treeTraits = new Helper();

    /**
     * The set of cases
     */
    protected AbstractOutbreak cases;

    // where along the relevant branches the infections happen. IMPORTANT: if extended=false then this should be
    // all 0s.

    protected Parameter infectionTimeBranchPositions;

    // where between infection and first child infection infectiousness happens

    protected Parameter infectiousTimePositions;

    protected double[] infectionTimes;
    private double[] storedInfectionTimes;
    protected double[] infectiousPeriods;
    private double[] storedInfectiousPeriods;
    protected double[] infectiousTimes;
    private double[] storedInfectiousTimes;
    protected double[] latentPeriods;
    private double[] storedLatentPeriods;

    //because of the way the former works, we need a maximum value of the time from first infection to root node.

    protected Parameter maxFirstInfToRoot;

    // latent periods

    protected boolean hasLatentPeriods;


    // PUBLIC STUFF

    // Name

    public static final String CASE_TO_CASE_TREE_LIKELIHOOD = "caseToCaseTreeLikelihood";
    public static final String PARTITIONS_KEY = "partition";


    // Basic constructor.

    public CaseToCaseTreeLikelihood(TreeModel virusTree, AbstractOutbreak caseData,
                                    Parameter infectionTimeBranchPositions, Parameter infectiousTimePositions,
                                    Parameter maxFirstInfToRoot)
            throws TaxonList.MissingTaxonException {
        this(CASE_TO_CASE_TREE_LIKELIHOOD, virusTree, caseData, infectionTimeBranchPositions, infectiousTimePositions,
                maxFirstInfToRoot);
    }

    // Constructor for an instance with a non-default name

    public CaseToCaseTreeLikelihood(String name, TreeModel virusTree, AbstractOutbreak caseData,
                                    Parameter infectionTimeBranchPositions, Parameter infectiousTimePositions,
                                    Parameter maxFirstInfToRoot) {
        super(name, caseData, virusTree);


        if(stateCount!=treeModel.getExternalNodeCount()){
            throw new RuntimeException("There are duplicate tip cases.");
        }

        noTips = virusTree.getExternalNodeCount();

        cases = caseData;

        addModel(cases);


        Date lastSampleDate = getLatestTaxonDate();
        estimatedLastSampleTime = lastSampleDate.getTimeValue();

        //map cases to tips

        branchMap = new BranchMapModel(virusTree);
        branchMap.setAll(new AbstractCase[treeModel.getNodeCount()]);
        addModel(branchMap);

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

        this.infectiousTimePositions = infectiousTimePositions;
        hasLatentPeriods = cases.hasLatentPeriods();
        if(hasLatentPeriods){
            addVariable(infectiousTimePositions);
        }

        this.infectionTimeBranchPositions = infectionTimeBranchPositions;
        addVariable(infectionTimeBranchPositions);

        if(DEBUG){
            for(int i=0; i<cases.size(); i++){
                if(!((CompoundParameter) infectionTimeBranchPositions).getParameter(i).getId()
                        .startsWith(cases.getCase(i).getName())){
                    throw new RuntimeException("Elements of outbreak and infectionTimeBranchPositions do not match up");
                }
            }
        }

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
    }

    protected void prepareTree(String startingNetworkFileName){
        if(startingNetworkFileName==null){
            partitionAccordingToRandomTT(true);
        } else {
            partitionAccordingToSpecificTT(startingNetworkFileName);
        }

        infectionTimes = getInfectionTimes(branchMap);
        infectiousPeriods = getInfectiousPeriods(branchMap);

        if(hasLatentPeriods){
            infectiousTimes = getInfectiousTimes(branchMap);
            latentPeriods = getLatentPeriods(branchMap);
        }
    }

    public AbstractOutbreak getOutbreak(){
        return cases;
    }

    public boolean hasLatentPeriods(){
        return hasLatentPeriods;
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

    private HashSet<Integer> samePartitionDownTree(NodeRef node, BranchMapModel map, boolean flagForRecalc){
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

    private HashSet<Integer> samePartitionUpTree(NodeRef node, BranchMapModel map, boolean flagForRecalc){
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

    // returns all nodes that are the earliest nodes in the partitions corresponding to this cases' children in
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
        NodeRef child = treeModel.getNode(tipMap.get(thisCase));
        NodeRef parent = treeModel.getParent(child);
        boolean transmissionFound = false;
        while(!transmissionFound){
            if(branchMap.get(child.getNumber())!=branchMap.get(parent.getNumber())){
                transmissionFound = true;
            } else {
                child = parent;
                parent = treeModel.getParent(child);
                if(parent == null){
                    transmissionFound = true;
                }
            }
        }
        return child;
    }

    protected NodeRef getEarliestNodeInPartition(AbstractCase thisCase){
        return getEarliestNodeInPartition(thisCase, branchMap);
    }

    public HashSet<AbstractCase> getDescendants(AbstractCase thisCase){
        HashSet<AbstractCase> out = new HashSet<AbstractCase>(getInfectees(thisCase));
        for(AbstractCase child : out){
            out.addAll(getDescendants(child));
        }
        return out;
    }


    public Integer[] getParentsArray(){
        Integer[] out = new Integer[cases.size()];
        for(AbstractCase thisCase : cases.getCases()){
            out[cases.getCaseIndex(thisCase)]=cases.getCaseIndex(getInfector(thisCase));
        }
        return out;
    }

    // @todo these currently return maps because C2CTransL may need them for the sorting, but I don't like it

    public HashMap<Integer, Double> getInfTimesMap(){
        if(infectionTimes==null){
            infectionTimes = getInfectionTimes(branchMap);
        }
        HashMap<Integer, Double> out = new HashMap<Integer, Double>();
        for(int i=0; i<cases.size(); i++){
            out.put(i, infectionTimes[i]);
        }
        return out;
    }

    public HashMap<Integer, Double> getInfnsTimesMap(){
        if(infectiousTimes==null){
            infectiousTimes = getInfectiousTimes(branchMap);
        }
        HashMap<Integer, Double> out = new HashMap<Integer, Double>();
        for(int i=0; i<cases.size(); i++){
            out.put(i, infectiousTimes[i]);
        }
        return out;
    }

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************


    protected void handleModelChangedEvent(Model model, Object object, int index) {

        if (model == treeModel) {

            // todo actually, most of these don't change

            infectionTimes = null;
            infectiousPeriods = null;
            if(hasLatentPeriods){
                infectiousTimes = null;
                latentPeriods = null;
            }
        }

        // todo this can be better

        if (model == branchMap){
            infectionTimes = null;
            infectiousPeriods = null;
            if(hasLatentPeriods){
                infectiousTimes = null;
                latentPeriods = null;
            }
        }

        fireModelChanged(object);

        likelihoodKnown = false;
    }


    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************


    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

        if(variable == infectionTimeBranchPositions){
            infectionTimes = null;
            infectiousPeriods = null;
            if(hasLatentPeriods){
                infectiousTimes = null;
                latentPeriods = null;
            }

        } else if(variable == infectiousTimePositions){
            infectiousPeriods = null;
            infectiousTimes = null;
            latentPeriods = null;
        }

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
        storedInfectionTimes = infectionTimes;
        storedInfectiousPeriods = infectiousPeriods;
        if(hasLatentPeriods){
            storedInfectiousTimes = infectiousTimes;
            storedLatentPeriods = latentPeriods;
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

    public final TreeModel getTreeModel(){
        return treeModel;
    }

    // todo this should now be part of handleModelChanged for the branchmap

    public void makeDirty() {
        likelihoodKnown = false;
        infectionTimes = null;
        infectiousPeriods = null;
        if(hasLatentPeriods){
            infectiousTimes = null;
            latentPeriods = null;
        }
    }


    protected void prepareTimings(){
        if(infectionTimes==null){
            infectionTimes = getInfectionTimes(branchMap);
        }

        if(hasLatentPeriods){
            if(infectiousTimes == null){
                infectiousTimes = getInfectiousTimes(branchMap);
            }
        }

        if(infectiousPeriods==null){
            infectiousPeriods = getInfectiousPeriods(branchMap);
        }

        if(hasLatentPeriods){
            if(latentPeriods == null){
                latentPeriods = getLatentPeriods(branchMap);
            }
        }
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
                double infectionTime = infectionTimes[cases.getCaseIndex(childCase)];
                if(infectionTime>parentCase.getCullDate().getTimeValue()
                        || hasLatentPeriods && infectionTime<infectiousTimes[cases.getCaseIndex(parentCase)]){
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
        NodeRef tip = treeModel.getNode(tipMap.get(thisCase));
        return getInfector(tip, branchMap);
    }

    public HashSet<AbstractCase> getInfectees(AbstractCase thisCase){
        return getInfectees(thisCase, branchMap);
    }

    public HashSet<AbstractCase> getInfectees(AbstractCase thisCase, BranchMapModel branchMap){
        return getInfecteesInClade(getEarliestNodeInPartition(thisCase), branchMap);
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

    // todo all these could be better


    public double getInfectionTime(AbstractCase thisCase){
        return infectionTimes[cases.getCaseIndex(thisCase)];
    }

    public double getInfectionTime(AbstractCase thisCase, BranchMapModel branchMap){
        if(infectionTimes!=null){
            return infectionTimes[cases.getCaseIndex(thisCase)];
        } else {
            NodeRef child = getEarliestNodeInPartition(thisCase);
            NodeRef parent = treeModel.getParent(child);

            if(parent!=null){
                AbstractCase parentCase = branchMap.get(parent.getNumber());
                double min = heightToTime(treeModel.getNodeHeight(parent));

                // Let the likelihood evaluate to zero due to culling dates if it must...

                double max = heightToTime(treeModel.getNodeHeight(child));


                return getInfectionTime(min, max, thisCase);
            } else {
                return getRootInfectionTime(branchMap);
            }
        }
    }

    private double getInfectionTime(double min, double max, AbstractCase infected){
        final double branchLength = max-min;
        return max - branchLength*infectionTimeBranchPositions.getParameterValue(cases.getCaseIndex(infected));
    }

    public double[] getInfectionTimes(){
        return infectionTimes;
    }

    public double[] getInfectionTimes(BranchMapModel branchMap){
        double[] out = new double[noTips];
        for(int i=0; i<noTips; i++){
            out[i] = getInfectionTime(cases.getCase(i), branchMap);
        }
        return out;
    }

    public double getInfectiousTime(AbstractCase thisCase){
        if(infectionTimes == null){
            throw new RuntimeException("Trying to get infectious times with null infection times array");
        }
        if(!hasLatentPeriods){
            return getInfectionTime(thisCase);
        } else if(infectiousTimes!=null){
            return infectiousTimes[cases.getCaseIndex(thisCase)];
        } else {
            HashSet<AbstractCase> infectees = getInfectees(thisCase);
            // needn't assume infectious at time of sampling, I don't think, but upper limit on infectious time is
            // cull time
            double latestInfectiousTime = thisCase.getCullTime();
            for(AbstractCase infectee: infectees){
                if(getInfectionTime(infectee)<latestInfectiousTime){
                    latestInfectiousTime = getInfectionTime(infectee);
                }
            }
            double infectionTime = getInfectionTime(thisCase);
            return infectionTime
                    + infectiousTimePositions.getParameterValue(cases.getCaseIndex(thisCase))
                    *(latestInfectiousTime-infectionTime);
        }
    }

    public double[] getInfectiousTimes(BranchMapModel branchMap){
        if(infectionTimes == null){
            infectionTimes = getInfectionTimes(branchMap);
        }
        double[] out = new double[noTips];
        for(int i=0; i<noTips; i++){
            out[i] = getInfectiousTime(cases.getCase(i));
        }
        return out;
    }

    public double[] getInfectiousTimes(){
        return infectiousTimes;
    }


    public double getInfectiousPeriod(AbstractCase thisCase){
        if(infectionTimes == null){
            throw new RuntimeException("Trying to get infectious periods before infection times");
        }
        if(!hasLatentPeriods){
            double infectionTime = getInfectionTime(thisCase);
            double cullTime = thisCase.getCullTime();
            return cullTime - infectionTime;
        } else {
            double infectiousTime = getInfectiousTime(thisCase);
            double cullTime = thisCase.getCullTime();
            return cullTime - infectiousTime;
        }
    }

    public double[] getInfectiousPeriods(BranchMapModel branchMap){
        if(infectionTimes == null){
            infectionTimes = getInfectionTimes(branchMap);
        }
        if(hasLatentPeriods && infectiousTimes == null){
            infectiousTimes = getInfectiousTimes(branchMap);
        }
        double[] out = new double[noTips];
        for(int i=0; i<noTips; i++){
            out[i] = getInfectiousPeriod(cases.getCase(i));
        }
        return out;
    }

    public double[] getInfectiousPeriods(){
        if(infectiousPeriods == null){
            infectiousPeriods = getInfectiousPeriods(branchMap);
        }
        return infectiousPeriods;
    }

    public double getLatentPeriod(AbstractCase thisCase){
        if(!hasLatentPeriods){
            return 0;
        } else if(latentPeriods!=null){
            return latentPeriods[cases.getCaseIndex(thisCase)];
        } else {
            int number = cases.getCaseIndex(thisCase);
            if(infectionTimes == null){
                infectionTimes = getInfectionTimes(branchMap);
            }
            if(infectiousTimes == null){
                infectiousTimes = getInfectiousTimes(branchMap);
            }
            return infectiousTimes[number] - infectionTimes[number];
        }
    }

    public double[] getLatentPeriods(BranchMapModel branchMap){
        if(infectionTimes == null){
            infectionTimes = getInfectionTimes(branchMap);
        }
        if(hasLatentPeriods && infectiousTimes == null){
            infectiousTimes = getInfectiousTimes(branchMap);
        }
        double[] out = new double[noTips];
        for(int i=0; i<noTips; i++){
            out[i] = getLatentPeriod(cases.getCase(i));
        }
        return out;
    }

    public double[] getLatentPeriods(){
        if(latentPeriods == null){
            latentPeriods = getLatentPeriods(branchMap);
        }
        return latentPeriods;
    }

    public double[] getInfectedPeriods(){
        if(infectionTimes == null){
            throw new RuntimeException("Trying to get infected periods with null infection times array");
        }
        if(!hasLatentPeriods){
            return infectiousPeriods;
        } else {
            double[] out = new double[noTips];
            for(int i=0; i<noTips; i++){
                out[i] = getInfectedPeriod(cases.getCase(i));
            }
            return out;
        }
    }

    public double getInfectedPeriod(int caseIndex){
        return cases.getCase(caseIndex).getCullTime() - infectionTimes[caseIndex];
    }

    public double getInfectedPeriod(AbstractCase thisCase){
        return getInfectedPeriod(cases.getCaseIndex(thisCase));
    }


    // return an array of the mean, median, variance and standard deviation of the given array
    // @todo this is pretty wasteful since it gets called so many times per log entry

    public static double[] getSummaryStatistics(double[] variable){
        DescriptiveStatistics stats = new DescriptiveStatistics(variable);
        double[] out = new double[4];
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
                + branchLength * infectionTimeBranchPositions.getParameterValue(cases.getCaseIndex(rootCase)));

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
            throw new RuntimeException("Tree is not partitioned properly");
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
        branchMap.setAll(prepareExternalNodeMap(new AbstractCase[treeModel.getNodeCount()]));

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
        branchMap.set(node.getNumber(), thisCase);
        if(tipLinked(node)){
            for(int i=0; i<treeModel.getChildCount(node); i++){
                specificallyPartitionUpwards(treeModel.getChild(node, i), thisCase, map);
            }
        } else {
            branchMap.set(node.getNumber(),null);
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
                branchMap.set(node.getNumber(), child);
            } else {
                branchMap.set(node.getNumber(), thisCase);
            }
            for(int i=0; i<treeModel.getChildCount(node); i++){
                specificallyPartitionUpwards(treeModel.getChild(node, i), branchMap.get(node.getNumber()), map);
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

    private void partitionAccordingToRandomTT(boolean checkNonZero){
        boolean gotOne = false;
        int tries = 1;
        System.out.println("Generating a random starting painting of the tree (checking nonzero likelihood for all " +
                "branches and repeating up to 100 times until a start with nonzero likelihood is found)");
        System.out.print("Attempt: ");
        while(!gotOne){

            likelihoodKnown = false;
            boolean failed = false;
            System.out.print(tries + "...");
            branchMap.setAll(prepareExternalNodeMap(new AbstractCase[treeModel.getNodeCount()]));
            //Warning - if the BadPartitionException in randomlyAssignNode might be caused by a bug rather than both
            //likelihoods rounding to zero, you want to stop catching this to investigate.

            try{
                partitionAccordingToRandomTT(branchMap, checkNonZero);
            } catch(BadPartitionException e){
                failed = true;
            }

            infectionTimes = getInfectionTimes(branchMap);
            if(hasLatentPeriods){
                infectiousTimes = getInfectiousTimes(branchMap);
            }

            infectiousPeriods = getInfectiousPeriods(branchMap);
            if(hasLatentPeriods){
                latentPeriods = getLatentPeriods(branchMap);
            }

            makeDirty();

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

                    branchLogLs[i]= !infector.culledYet(nodeTime + infectionTimeBranchPositions.getParameterValue(
                            cases.getCaseIndex(infectee))*branchLength);
                }
                if(!branchLogLs[0] && !branchLogLs[1]){
                    throw new BadPartitionException("Both branch possibilities have zero likelihood: "
                            +node.toString()+", cases " + choices[0].getName() + " and " + choices[1].getName() + ".");
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
            map.set(node.getNumber(), winner);
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

        for(AbstractCase aCase : cases.getCases()){
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
                throw new RuntimeException("Rewiring messed up; investigate");
            }

        }


        return outTree;
    }

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








