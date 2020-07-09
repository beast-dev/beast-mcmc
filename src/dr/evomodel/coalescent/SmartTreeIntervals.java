/*
 * AbstractCoalescentLikelihood.java
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

package dr.evomodel.coalescent;

import dr.evolution.coalescent.IntervalList;
import dr.evolution.coalescent.IntervalType;
import dr.evolution.coalescent.Intervals;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;

import java.util.*;


/**
 * Forms a base class for a number of coalescent likelihood calculators.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author JT McCrone
 * @version $Id: CoalescentLikelihood.java,v 1.43 2006/07/28 11:27:32 rambaut Exp $
 */
public class SmartTreeIntervals extends AbstractModel implements Units, IntervalList {

    // PUBLIC STUFF

    public SmartTreeIntervals(
            Tree tree,
            TaxonList includeSubtree,
            List<TaxonList> excludeSubtrees) throws TreeUtils.MissingTaxonException {

        super("SmartTreeIntervals");

        this.tree = tree;

        if (includeSubtree != null) {
            includedLeafSet = TreeUtils.getLeavesForTaxa(tree, includeSubtree);
        }

        if (excludeSubtrees != null) {
            excludedLeafSets = new Set[excludeSubtrees.size()];
            for (int i = 0; i < excludeSubtrees.size(); i++) {
                excludedLeafSets[i] = TreeUtils.getLeavesForTaxa(tree, excludeSubtrees.get(i));
            }
        }
//        } else {
//            excludedLeafSets = new Set[0];
//        }

        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }

        intervals = new Intervals(tree.getNodeCount());
        storedIntervals = new Intervals(tree.getNodeCount());
        eventsKnown = false;
        fullUpdate = true;

        sortedNodeList = new ArrayList<>();
        storedSortedNodeList = new ArrayList<>();

        updateNode = new boolean[tree.getNodeCount()];
        storedUpdateNode = new boolean[tree.getNodeCount()];
        for (int i = 0; i < tree.getNodeCount(); i++) {
            updateNode[i] = true;
            sortedNodeList.add(i);
        }

        updateSortedNodeList();
        addStatistic(new DeltaStatistic());
    }

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {


        if (model == tree) {
            if (object instanceof TreeChangedEvent) {

                if (((TreeChangedEvent) object).isNodeChanged()) {
                    // If a node event occurs the node and its two child nodes
                    // are flagged for updating (this will result in everything
                    // above being updated as well. Node events occur when a node
                    // is added to a branch, removed from a branch or its height or
                    // rate changes.
                    NodeRef node = ((TreeChangedEvent) object).getNode();
                    updateNodeEvent(node);

                } else if (((TreeChangedEvent) object).isTreeChanged()) {
                    // Full tree events result in a complete updating of the tree likelihood
                    // This event type is now used for EmpiricalTreeDistributions.
//                    System.err.println("Full tree update event - these events currently aren't used\n" +
//                            "so either this is in error or a new feature is using them so remove this message.");
                    updateAllNodeEvents();
                }  // Other event types are ignored (probably trait changes).
                //System.err.println("Another tree event has occured (possibly a trait change).");

            }

            fireModelChanged();
        }
    }

    private void updateNodeEvent(NodeRef node) {
        updateNode[node.getNumber()] = true;
        eventsKnown = false;
    }

    private void updateAllNodeEvents() {
        for (int i = 0; i < tree.getNodeCount(); i++) {
            updateNode[i] = false;
        }
        fullUpdate = true;
        eventsKnown = true;
    }

    public void makeDirty(){
        updateAllNodeEvents();
    }

    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
    } // No parameters to respond to

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    /**
     * Stores the precalculated state: in this case the intervals
     */
    protected void storeState() {
        // copy the intervals into the storedIntervals
        storedIntervals.copyIntervals(intervals);

        storedSortedNodeList = new ArrayList<>(sortedNodeList);
        System.arraycopy(updateNode, 0, storedUpdateNode, 0, updateNode.length);

        storedEventsKnown = eventsKnown;
        storedFullUpdate = fullUpdate;
    }

    /**
     * Restores the precalculated state: that is the intervals of the tree.
     */
    protected void restoreState() {
        // swap the intervals back
        Intervals tmp = storedIntervals;
        storedIntervals = intervals;
        intervals = tmp;

        boolean[] tmp2 = storedUpdateNode;
        storedUpdateNode = updateNode;
        updateNode = tmp2;

        List<Integer> tmp3 = storedSortedNodeList;
        storedSortedNodeList = sortedNodeList;
        sortedNodeList = tmp3;


        eventsKnown = storedEventsKnown;
        fullUpdate = storedFullUpdate;
    }

    protected final void acceptState() {
    } // nothing to do

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    public final Model getModel() {
        return this;
    }

    /**
     * @return the node ref of the MRCA of this coalescent prior in the given tree.
     */
    private NodeRef getIncludedMRCA(Tree tree) {
        if (includedLeafSet != null) {
            return TreeUtils.getCommonAncestorNode(tree, includedLeafSet);
        } else {
            return tree.getRoot();
        }
    }

    /**
     * @return an array of noderefs that represent the MRCAs of subtrees to exclude from coalescent prior.
     * May return null if no subtrees should be excluded.
     */
    private Set<NodeRef> getExcludedMRCAs(Tree tree) {

        if (excludedLeafSets.length == 0) return null;

        Set<NodeRef> excludeNodesBelow = new HashSet<NodeRef>();
        for (Set<String> excludedLeafSet : excludedLeafSets) {
            excludeNodesBelow.add(TreeUtils.getCommonAncestorNode(tree, excludedLeafSet));
        }
        return excludeNodesBelow;
    }

    public Tree getTree() {
        return tree;
    }

    public IntervalList getIntervals() {
        return intervals;
    }

    /**
     * Recalculates all the intervals from the tree model.
     */
    public final void calculateIntervals() {

        intervals.resetEvents();
        if (includedLeafSet != null || excludedLeafSets != null) {
            collectTimes(tree, getIncludedMRCA(tree), getExcludedMRCAs(tree), intervals);
        } else {
            collectTimes(tree, intervals);
        }

        // force a calculation of the intervals...
        intervals.calculateIntervals(true);

        eventsKnown = true;
    }


    /**
     * extract coalescent times and tip information into ArrayList times from tree.
     *
     * @param tree      the tree
     * @param node      the node to start from
     * @param intervals the intervals object to store the events
     */
    private void collectTimes(Tree tree, NodeRef node, Set<NodeRef> excludeNodesBelow, Intervals intervals) {

        intervals.addCoalescentEvent(tree.getNodeHeight(node));

        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef child = tree.getChild(node, i);

            // check if this subtree is included in the coalescent density
            boolean include = true;

            if (excludeNodesBelow != null && excludeNodesBelow.contains(child)) {
                include = false;
            }

            if (!include || tree.isExternal(child)) {
                intervals.addSampleEvent(tree.getNodeHeight(child));
            } else {
                collectTimes(tree, child, excludeNodesBelow, intervals);
            }
        }

    }

    /**
     * An alternative non-recursive version (doesn't exclude nodes).
     *
     * @param tree      the tree
     * @param intervals the intervals object to store the events
     */
    private void collectTimes(Tree tree, Intervals intervals) {
        updateSortedNodeList();

        for (Integer nodeIndex : sortedNodeList) {
            NodeRef node = tree.getNode(nodeIndex);
            if (tree.isExternal(node)) {
                intervals.addSampleEvent(tree.getNodeHeight(node));
            } else {
                intervals.addCoalescentEvent(tree.getNodeHeight(node));
            }
        }
    }


    private void updateSortedNodeList() {
        if (fullUpdate) {
            //TODO add type for ties (.thenComparing)
            sortedNodeList.sort(Comparator.comparing(n -> tree.getNodeHeight(tree.getNode(n))));
            fullUpdate = false;
            for (int i = 0; i < updateNode.length; i++) {
                updateNode[i] = false;
            }
            return;
        }

        List<Integer> updatedNodes = new ArrayList<>();
        List<Integer> nodeIndexes = new ArrayList<>();

        for (int i = 0; i < updateNode.length; i++) {
            if (updateNode[i]) {
                int nodeIndex = sortedNodeList.indexOf(i);
                sortedNodeList.remove(nodeIndex);
                updatedNodes.add(i);
                nodeIndexes.add(nodeIndex);

                for (int j = 0; j < updatedNodes.size(); j++) {
                    Integer originalIndex = nodeIndexes.get(j);
                    if (nodeIndex < nodeIndexes.get(j)) {
                        originalIndex = originalIndex - 1;
                    }
                    if (originalIndex >= sortedNodeList.size()) {
                        originalIndex = originalIndex - 1;
                    }
                    nodeIndexes.set(j, originalIndex);
                }
                updateNode[i] = false;
            }
        }

        for (int i = 0; i < updatedNodes.size(); i++) {
            int currentIndex = nodeIndexes.get(i);
            Integer nodeIndex = updatedNodes.get(i);

            double height = tree.getNodeHeight(tree.getNode(nodeIndex));

            double currentHeight = tree.getNodeHeight(tree.getNode(sortedNodeList.get(currentIndex)));
            while (currentHeight < height && currentIndex < sortedNodeList.size()) {
                currentHeight = tree.getNodeHeight(tree.getNode(sortedNodeList.get(currentIndex)));
                currentIndex++;
            }

            if(currentHeight>height){
                while (currentHeight > height && currentIndex > -1) {
                    currentIndex--;
                    currentHeight = tree.getNodeHeight(tree.getNode(sortedNodeList.get(currentIndex)));
                }
                currentIndex++;
            }


            sortedNodeList.add(currentIndex, nodeIndex);
        }

    }

    @Override
    public double getStartTime() {
        return startTime;
    }

    @Override
    public int getIntervalCount() {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.getIntervalCount();
    }

    @Override
    public int getSampleCount() {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.getSampleCount();
    }

    @Override
    public double getInterval(int i) {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.getInterval(i);
    }

    @Override
    public double getIntervalTime(int i) {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.getIntervalTime(i);
    }

    @Override
    public int getLineageCount(int i) {
        if (!eventsKnown) {
            calculateIntervals();
        }
        if (i >= getIntervalCount()) throw new IllegalArgumentException();
        return intervals.getLineageCount(i);
    }

    @Override
    public int getCoalescentEvents(int i) {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.getCoalescentEvents(i);
    }

    @Override
    public IntervalType getIntervalType(int i) {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.getIntervalType(i);
    }

    @Override
    public double getTotalDuration() {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.getTotalDuration();
    }

    @Override
    public boolean isBinaryCoalescent() {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.isBinaryCoalescent();
    }

    @Override
    public boolean isCoalescentOnly() {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.isCoalescentOnly();
    }

    @Override
    public Type getUnits() {
        return intervals.getUnits();
    }

    @Override
    public void setUnits(Type units) {
        intervals.setUnits(units);
    }

    // ****************************************************************
    // Inner classes
    // ****************************************************************

    public class DeltaStatistic extends Statistic.Abstract {

        public DeltaStatistic() {
            super("delta");
        }

        public int getDimension() {
            return 1;
        }

        public double getStatisticValue(int i) {
            throw new RuntimeException("Not implemented");
//			return IntervalList.Utils.getDelta(intervals);
        }

    }

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    /**
     * The tree.
     */
    private Tree tree = null;
    private Set<String> includedLeafSet = null;
    private Set[] excludedLeafSets = null;

    private double startTime;
    private List<Integer> sortedNodeList;

    /**
     * The intervals.
     */
    private Intervals intervals = null;

    /**
     * The stored values for intervals.
     */
    private Intervals storedIntervals = null;
    private List<Integer> storedSortedNodeList;
    private boolean[] updateNode;
    private boolean[] storedUpdateNode;

    private boolean eventsKnown = false;
    private boolean storedEventsKnown = false;
    private boolean fullUpdate = true;
    private boolean storedFullUpdate = true;

}