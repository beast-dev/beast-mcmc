/*
 * CaseToCaseTreeLikelihood.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.epidemiology.casetocase;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

import dr.app.tools.NexusExporter;
import dr.evolution.tree.FlexibleNode;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.epidemiology.casetocase.periodpriors.AbstractPeriodPriorDistribution;
import dr.evomodel.tree.TreeModel;
import dr.oldevomodel.treelikelihood.AbstractTreeLikelihood;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import org.apache.commons.math.stat.descriptive.moment.Mean;
import org.apache.commons.math.stat.descriptive.moment.Variance;
import org.apache.commons.math.stat.descriptive.rank.Median;

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

    protected static double tolerance = 1E-10;

    /* The phylogenetic tree. */

    protected int noTips;
    protected int noCases;

    /* Mapping of outbreak to branches on the tree; old version is stored before operators are applied */


    /* Matches outbreak to external nodes */

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

    protected HashMap<AbstractCase,Treelet> elementsAsTrees;
    protected HashMap<AbstractCase,Treelet> storedElementsAsTrees;

    //because of the way the former works, we need a maximum value of the time from first infection to root node.

    protected Parameter maxFirstInfToRoot;

    // latent periods

    protected boolean hasLatentPeriods;


    // PUBLIC STUFF

    // Name

    public static final String CASE_TO_CASE_TREE_LIKELIHOOD = "caseToCaseTreeLikelihood";
    public static final String PARTITIONS_KEY = "partition";


    // Basic constructor.

    public CaseToCaseTreeLikelihood(PartitionedTreeModel tree, AbstractOutbreak caseData,
                                    Parameter maxFirstInfToRoot)
            throws TaxonList.MissingTaxonException {
        this(CASE_TO_CASE_TREE_LIKELIHOOD, tree, caseData, maxFirstInfToRoot);
    }

    // Constructor for an instance with a non-default name

    public CaseToCaseTreeLikelihood(String name, PartitionedTreeModel tree, AbstractOutbreak caseData,
                                    Parameter maxFirstInfToRoot) {


        super(name, caseData, tree);


        if(stateCount!=treeModel.getExternalNodeCount()){
            throw new RuntimeException("There are duplicate tip outbreak.");
        }

        noTips = tree.getExternalNodeCount();


        //subclasses should add outbreak as a model if it contains any information that ever changes

        outbreak = caseData;

        noCases = outbreak.getCases().size();

        addModel(outbreak);

        estimatedLastSampleTime = getLatestTaxonTime();

        //map outbreak to tips

        addModel(tree.getBranchMap());

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



    public AbstractOutbreak getOutbreak(){
        return outbreak;
    }

    public boolean hasLatentPeriods(){
        return hasLatentPeriods;
    }

    /* Get the date of the last tip */

    private double getLatestTaxonTime(){
        double latestTime = Double.NEGATIVE_INFINITY;
        for(int i=0; i<treeModel.getExternalNodeCount(); i++){
            Taxon taxon = treeModel.getNodeTaxon(treeModel.getExternalNode(i));
            if(taxon.getDate().getTimeValue() > latestTime){
                latestTime = taxon.getDate().getTimeValue();
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



    protected void explodeTree(){

        for(int i=0; i<outbreak.size(); i++){
            AbstractCase aCase = outbreak.getCase(i);
            if(aCase.wasEverInfected() && elementsAsTrees.get(aCase)==null){

                NodeRef partitionRoot = ((PartitionedTreeModel)treeModel).getEarliestNodeInElement(aCase);

                double extraHeight;

                if(treeModel.isRoot(partitionRoot)){
                    extraHeight = maxFirstInfToRoot.getParameterValue(0)
                            * aCase.getInfectionBranchPosition().getParameterValue(0);
                } else {
                    extraHeight = treeModel.getBranchLength(partitionRoot)
                            * aCase.getInfectionBranchPosition().getParameterValue(0);
                }

                FlexibleNode newRoot = new FlexibleNode();

                FlexibleTree littleTree = new FlexibleTree(newRoot);
                littleTree.beginTreeEdit();

                if (!treeModel.isExternal(partitionRoot)) {
                    for (int j = 0; j < treeModel.getChildCount(partitionRoot); j++) {
                        copyElementToTreelet(littleTree, treeModel.getChild(partitionRoot, j), newRoot, aCase);
                    }
                }

                littleTree.endTreeEdit();

                littleTree.resolveTree();

                Treelet treelet = new Treelet(littleTree,
                        littleTree.getRootHeight() + extraHeight);

                elementsAsTrees.put(aCase, treelet);
            }
        }
    }

    private void copyElementToTreelet(FlexibleTree littleTree, NodeRef oldNode, NodeRef newParent,
                                      AbstractCase element){
        if(element.wasEverInfected()) {
            if (getBranchMap().get(oldNode.getNumber()) == element) {
                if (treeModel.isExternal(oldNode)) {
                    NodeRef newTip = new FlexibleNode(new Taxon(treeModel.getNodeTaxon(oldNode).getId()));
                    littleTree.addChild(newParent, newTip);
                    littleTree.setBranchLength(newTip, treeModel.getBranchLength(oldNode));
                } else {
                    NodeRef newChild = new FlexibleNode();
                    littleTree.addChild(newParent, newChild);
                    littleTree.setBranchLength(newChild, treeModel.getBranchLength(oldNode));
                    for (int i = 0; i < treeModel.getChildCount(oldNode); i++) {
                        copyElementToTreelet(littleTree, treeModel.getChild(oldNode, i), newChild, element);
                    }
                }
            } else {
                // we need a new tip
                NodeRef transmissionTip = new FlexibleNode(
                        new Taxon("Transmission_" + getBranchMap().get(oldNode.getNumber()).getName()));
                double parentTime = getNodeTime(treeModel.getParent(oldNode));
                double childTime = getInfectionTime(getBranchMap().get(oldNode.getNumber()));
                littleTree.addChild(newParent, transmissionTip);
                littleTree.setBranchLength(transmissionTip, childTime - parentTime);
            }
        }
    }

    protected class Treelet extends FlexibleTree {

        private double zeroHeight;

        protected Treelet(FlexibleTree tree, double zeroHeight){
            super(tree);
            this.zeroHeight = zeroHeight;

        }

        protected double getZeroHeight(){
            return zeroHeight;
        }

        protected void setZeroHeight(double rootBranchLength){
            this.zeroHeight = zeroHeight;
        }
    }



    // find all partitions of the descendant tips of the current node. If map is specified then it makes a map of node
    // number to possible partitions; map can be null.

    public HashSet<AbstractCase> descendantTipPartitions(NodeRef node, HashMap<Integer, HashSet<AbstractCase>> map){
        HashSet<AbstractCase> out = new HashSet<AbstractCase>();
        if(treeModel.isExternal(node)){
            out.add(getBranchMap().get(node.getNumber()));
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
            } else if (model == getBranchMap()){
                if(object instanceof ArrayList){

                    for(int i=0; i<((ArrayList) object).size(); i++){
                        BranchMapModel.BranchMapChangedEvent event
                                =  (BranchMapModel.BranchMapChangedEvent)((ArrayList) object).get(i);

                        recalculateCase(event.getOldCase());
                        recalculateCase(event.getNewCase());

                        NodeRef node = treeModel.getNode(event.getNodeToRecalculate());
                        NodeRef parent = treeModel.getParent(node);

                        if(parent!=null){
                            recalculateCase(getBranchMap().get(parent.getNumber()));
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
        return ((PartitionedTreeModel)treeModel).getBranchMap();
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
            AbstractCase childCase = getBranchMap().get(node.getNumber());
            AbstractCase parentCase = getBranchMap().get(treeModel.getParent(node).getNumber());
            if(childCase!=parentCase){
                double infectionTime = infectionTimes[outbreak.getCaseIndex(childCase)];
                if(infectionTime>parentCase.getEndTime()
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



    public double getInfectionTime(AbstractCase thisCase){

        if(!recalculateCaseFlags[outbreak.getCaseIndex(thisCase)]){
            return infectionTimes[outbreak.getCaseIndex(thisCase)];
        } else {
            if(thisCase.wasEverInfected()) {
                NodeRef child = ((PartitionedTreeModel)treeModel).getEarliestNodeInElement(thisCase);
                NodeRef parent = treeModel.getParent(child);

                if (parent != null) {

                    double min = heightToTime(treeModel.getNodeHeight(parent));

                    // Let the likelihood evaluate to zero due to culling dates if it must...

                    double max = heightToTime(treeModel.getNodeHeight(child));


                    return getInfectionTime(min, max, thisCase);
                } else {
                    return getRootInfectionTime(getBranchMap());
                }
            } else {
                return Double.POSITIVE_INFINITY;
            }
        }
    }

    private double getInfectionTime(double min, double max, AbstractCase infected){
        final double branchLength = max-min;
        return min + branchLength*(1-infected.getInfectionBranchPosition().getParameterValue(0));
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
            NodeRef child = ((PartitionedTreeModel)treeModel).getEarliestNodeInElement(thisCase);
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
                    double cullTime = thisCase.getEndTime();
                    infectiousPeriods[outbreak.getCaseIndex(thisCase)] = cullTime - infectionTime;
                } else {
                    double infectiousTime = getInfectiousTime(thisCase);
                    double cullTime = thisCase.getEndTime();
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
            return thisCase.getEndTime() - getInfectionTime(thisCase);
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

        Double[] out = new Double[4];
        out[0] = (new Mean()).evaluate(primitiveVariable);
        out[1] = (new Median()).evaluate(primitiveVariable);
        out[2] = (new Variance()).evaluate(primitiveVariable);
        out[3] = Math.sqrt(out[2]);
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
        AbstractCase rootCase = getBranchMap().get(treeModel.getRoot().getNumber());
        return getInfectionTime(rootCase);
    }

    public void outputTreeToFile(String fileName, boolean includeTransmissionNodes){
        outputTreeToFile(getBranchMap(), fileName, includeTransmissionNodes);
    }


    public void outputTreeToFile(BranchMapModel map, String fileName, boolean includeTransmissionNodes){
        try{
            FlexibleTree treeCopy;
            if(!includeTransmissionNodes){
                treeCopy = new FlexibleTree(treeModel);
                for(int j=0; j<treeCopy.getNodeCount(); j++){
                    FlexibleNode node = (FlexibleNode)treeCopy.getNode(j);
                    node.setAttribute("Number", node.getNumber());
                    node.setAttribute("Time", heightToTime(node.getHeight()));
                    node.setAttribute(PARTITIONS_KEY, map.get(node.getNumber()));
                }
            } else {
                treeCopy = addTransmissionNodes(treeModel);
            }
            NexusExporter testTreesOut = new NexusExporter(new PrintStream(fileName));
            testTreesOut.exportTree(treeCopy);
        } catch (IOException ignored) {System.out.println("IOException");}
    }

    public FlexibleTree addTransmissionNodes(Tree tree){
        prepareTimings();

        FlexibleTree outTree = new FlexibleTree(tree, true);

        for(int j=0; j<outTree.getNodeCount(); j++){
            FlexibleNode node = (FlexibleNode)outTree.getNode(j);
            node.setAttribute("Number", node.getNumber());
            node.setAttribute("Time", heightToTime(node.getHeight()));
            node.setAttribute(PARTITIONS_KEY, getBranchMap().get(node.getNumber()));
        }

        for(AbstractCase aCase : outbreak.getCases()){
            if(aCase.wasEverInfected()) {
                NodeRef originalNode = ((PartitionedTreeModel)treeModel).getEarliestNodeInElement(aCase);

                int infectionNodeNo = originalNode.getNumber();
                if (!treeModel.isRoot(originalNode)) {
                    NodeRef originalParent = treeModel.getParent(originalNode);
                    double nodeTime = getNodeTime(originalNode);
                    double infectionTime = getInfectionTime(aCase);
                    double heightToBreakBranch = getHeight(originalNode) + (nodeTime - infectionTime);
                    FlexibleNode newNode = (FlexibleNode) outTree.getNode(infectionNodeNo);
                    FlexibleNode oldParent = (FlexibleNode) outTree.getParent(newNode);

                    outTree.beginTreeEdit();
                    outTree.removeChild(oldParent, newNode);
                    FlexibleNode infectionNode = new FlexibleNode();
                    infectionNode.setHeight(heightToBreakBranch);
                    infectionNode.setLength(oldParent.getHeight() - heightToBreakBranch);
                    infectionNode.setAttribute(PARTITIONS_KEY, getNodePartition(treeModel, originalParent));
                    infectionNode.setAttribute("Time", heightToTime(heightToBreakBranch));
                    newNode.setLength(nodeTime - infectionTime);

                    outTree.addChild(oldParent, infectionNode);
                    outTree.addChild(infectionNode, newNode);
                    outTree.endTreeEdit();
                } else {
                    double nodeTime = getNodeTime(originalNode);
                    double infectionTime = getInfectionTime(aCase);
                    double heightToInstallRoot = getHeight(originalNode) + (nodeTime - infectionTime);
                    FlexibleNode newNode = (FlexibleNode) outTree.getNode(infectionNodeNo);
                    outTree.beginTreeEdit();
                    FlexibleNode infectionNode = new FlexibleNode();
                    infectionNode.setHeight(heightToInstallRoot);
                    infectionNode.setAttribute("Time", heightToTime(heightToInstallRoot));
                    infectionNode.setAttribute(PARTITIONS_KEY, "Origin");
                    outTree.addChild(infectionNode, newNode);
                    newNode.setLength(heightToInstallRoot - getHeight(originalNode));
                    outTree.setRoot(infectionNode);
                    outTree.endTreeEdit();
                }
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
                    ((PartitionedTreeModel)treeModel).checkPartitions();
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
        int count = 0;
        for(int i=0; i<outbreak.size(); i++){
            final AbstractCase infected = outbreak.getCase(i);
            if(infected.wasEverInfected()) {
                columns[count] = new LogColumn.Abstract(infected.toString() + "_infector") {
                    protected String getFormattedValue() {
                        if (((PartitionedTreeModel)treeModel).getInfector(infected) == null) {
                            return "Start";
                        } else {
                            return ((PartitionedTreeModel)treeModel).getInfector(infected).toString();
                        }
                    }
                };
                count++;
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

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "Case to Case Transmission Tree model";
    }

    public List<Citation> getCitations() {
        return Arrays.asList(new Citation(
                        new Author[]{new Author("M", "Hall"), new Author("M", "Woolhouse"), new Author("A", "Rambaut")},
            "Epidemic Reconstruction in a Phylogenetics Framework: Transmission Trees as Partitions of the Node Set",
            2016, "PLOS Comput Biol",
                11,
                0, 0, "10.1371/journal.pcbi.1004613",
                Citation.Status.PUBLISHED));
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
                    return getBranchMap().get(oldNode.getNumber()).toString();
                }
            } catch(NullPointerException e){
                if(tree.isRoot(node)){
                    return "Start";
                } else {
                    NodeRef parent = tree.getParent(node);
                    int originalParentNumber = (Integer)tree.getNodeAttribute(parent,"Number");
                    return getBranchMap().get(originalParentNumber).toString();
                }
            }
        } else {
            return getBranchMap().get(node.getNumber()).toString();
        }
    }

    public Integer[] getParentsArray(){
        Integer[] out = new Integer[outbreak.size()];
        for(AbstractCase thisCase : outbreak.getCases()){
            if(thisCase.wasEverInfected()) {
                out[outbreak.getCaseIndex(thisCase)] = outbreak.getCaseIndex(((PartitionedTreeModel)treeModel).getInfector(thisCase));
            } else {
                out[outbreak.getCaseIndex(thisCase)] = null;
            }
        }
        return out;
    }

    public AbstractCase getInfector(int i){
        return ((PartitionedTreeModel)treeModel).getInfector(getOutbreak().getCase(i));
    }

}








