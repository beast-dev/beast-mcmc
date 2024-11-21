/*
 * TreeIntervals.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.coalescent;

import dr.evolution.coalescent.*;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import java.util.*;


/**
 * Forms a base class for a number of coalescent likelihood calculators.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class TreeIntervals extends AbstractModel implements Units, TreeIntervalList {

    // PUBLIC STUFF


    public TreeIntervals(Tree tree) {
            this(tree, true);
    }
    public TreeIntervals(Tree tree, boolean buildIntervalNodeMapping) {
        super("TreeIntervals");

        this.buildIntervalNodeMapping = buildIntervalNodeMapping;

        setup(tree);
    }

    public TreeIntervals(
            Tree tree,
            TaxonList includeSubtree,
            List<TaxonList> excludeSubtrees) throws TreeUtils.MissingTaxonException {

        this(tree, includeSubtree, excludeSubtrees, false);
    }

    public TreeIntervals(
            Tree tree,
            TaxonList includeSubtree,
            List<TaxonList> excludeSubtrees,
            boolean buildIntervalNodeMapping) throws TreeUtils.MissingTaxonException {

        super("TreeIntervals");

        this.buildIntervalNodeMapping = buildIntervalNodeMapping;

        if (includeSubtree != null) {
            includedLeafSet = TreeUtils.getLeavesForTaxa(tree, includeSubtree);
        }

        if (excludeSubtrees != null && excludeSubtrees.size() > 0) {
            excludedLeafSets = new Set[excludeSubtrees.size()];
            for (int i = 0; i < excludeSubtrees.size(); i++) {
                excludedLeafSets[i] = TreeUtils.getLeavesForTaxa(tree, excludeSubtrees.get(i));
            }
        } else {
            excludedLeafSets = null;
        }

        setup(tree);
    }

    private void setup(Tree tree) {
        this.tree = tree;

        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }

        if (this.buildIntervalNodeMapping) {
            intervals = new Intervals(tree.getNodeCount());
            storedIntervals = new Intervals(tree.getNodeCount());
        } else {
            intervals = new FastIntervals(tree.getExternalNodeCount(), tree.getInternalNodeCount());
            storedIntervals = new FastIntervals(tree.getExternalNodeCount(), tree.getInternalNodeCount());
        }

        eventsKnown = false;
        this.intervalNodeMapping = buildIntervalNodeMapping ?new IntervalNodeMapProvider.IntervalNodeMapping.Default(tree.getNodeCount(), tree):new IntervalNodeMapProvider.IntervalNodeMapping.None();

        addStatistic(new DeltaStatistic());
    }

   public void setBuildIntervalNodeMapping(boolean buildIntervalNodeMapping){
       this.buildIntervalNodeMapping = buildIntervalNodeMapping;
       if (this.buildIntervalNodeMapping) {
        intervals = new Intervals(tree.getNodeCount());
        storedIntervals = new Intervals(tree.getNodeCount());
        } else {
            intervals = new FastIntervals(tree.getExternalNodeCount(), tree.getInternalNodeCount());
            storedIntervals = new FastIntervals(tree.getExternalNodeCount(), tree.getInternalNodeCount());
        }
       //Force a recalculation here
        eventsKnown = false;
        this.intervalNodeMapping = buildIntervalNodeMapping ?new IntervalNodeMapProvider.IntervalNodeMapping.Default(tree.getNodeCount(), tree):new IntervalNodeMapProvider.IntervalNodeMapping.None();
   }
    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == tree) {
            // treeModel has changed so recalculate the intervals
            eventsKnown = false;
        }

        fireModelChanged();
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
        storedEventsKnown = eventsKnown;
        this.intervalNodeMapping.storeMapping();
    }

    /**
     * Restores the precalculated state: that is the intervals of the tree.
     */
    protected void restoreState() {
        // swap the intervals back
        MutableIntervalList tmp = storedIntervals;
        storedIntervals = intervals;
        intervals = tmp;

        eventsKnown = storedEventsKnown;
        this.intervalNodeMapping.restoreMapping();

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

        if (excludedLeafSets == null || excludedLeafSets.length == 0) return null;

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
    //TODO pull node interval map stuff out of old Abtract coalescentLikelihood into here
    public final void calculateIntervals() {

        intervals.resetEvents();
        if (includedLeafSet != null || excludedLeafSets != null) {
            collectTimes(tree, getIncludedMRCA(tree), getExcludedMRCAs(tree), intervals);
        } else {
            collectTimes(tree, intervals);
        }

        // force a calculation of the intervals...
        intervals.getIntervalCount();

        if (buildIntervalNodeMapping){
            this.intervalNodeMapping.initializeMaps();
            for(int i=0; i<intervals.getIntervalCount()+1;i++){
                intervalNodeMapping.addNode(intervals.getNodeForEvent(i));
                if(i>0&& i<intervals.getIntervalCount()){ //If the event is not the first but not the last add it again for the start of the next interval
                    intervalNodeMapping.addNode(intervals.getNodeForEvent(i));
                }
            }
            intervalNodeMapping.setIntervalStartIndices(intervals.getIntervalCount());
        }
        eventsKnown = true;
    }

    /**
     * extract coalescent times and tip information into ArrayList times from tree.
     *
     * @param tree      the tree
     * @param node      the node to start from
     * @param intervals the intervals object to store the events
     */
    private void collectTimes(Tree tree, NodeRef node, Set<NodeRef> excludeNodesBelow, MutableIntervalList intervals) {

        intervals.addCoalescentEvent(tree.getNodeHeight(node));

        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef child = tree.getChild(node, i);

            // check if this subtree is included in the coalescent density
            boolean include = true;

            if (excludeNodesBelow != null && excludeNodesBelow.contains(child)) {
                include = false;
            }

            if (!include || tree.isExternal(child)) {
                // the mrca of the clade below that is being excluded becomes a sampling event
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
    private void collectTimes(Tree tree, MutableIntervalList intervals) {

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            intervals.addSampleEvent(tree.getNodeHeight(tree.getExternalNode(i)), tree.getExternalNode(i).getNumber());
        }
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            intervals.addCoalescentEvent(tree.getNodeHeight(tree.getInternalNode(i)), tree.getInternalNode(i).getNumber());
        }
    }

    @Override
    public double getStartTime() {
        return intervals.getStartTime();
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

    public IntervalNodeMapProvider getIntervalNodeMap() {
        return intervalNodeMapping;
    }

    // Interval Node mapping


    @Override
    public int[] getIntervalsForNode(int nodeNumber) {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return this.intervalNodeMapping.getIntervalsForNode(nodeNumber);
    }

    @Override
    public int[] getNodeNumbersForInterval(int interval) {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return this.intervalNodeMapping.getNodeNumbersForInterval(interval);
    }

    @Override
    public boolean isBuildIntervalNodeMapping() {
        return buildIntervalNodeMapping;
    }

    @Override
    public NodeRef getCoalescentNode(int interval) {
        if (!eventsKnown) {
            calculateIntervals();
        }
        if(getIntervalType(interval)!=IntervalType.COALESCENT) throw new IllegalArgumentException("Not a coalescent interval");
        return tree.getNode(this.intervalNodeMapping.getNodeNumbersForInterval(interval)[1]); //TODO verify index;
    }

    @Override
    public double[] sortByNodeNumbers(double[] unSortedNodeHeightGradient) {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return this.intervalNodeMapping.sortByNodeNumbers(unSortedNodeHeightGradient);
    }

    @Override
    public double[] getCoalescentIntervals() {
        double[] coalIntervals = new double[tree.getInternalNodeCount()];
        int currentIndex = 0;
        for (int i = 0; i < this.intervals.getIntervalCount(); i++) {
            if (this.getIntervalType(i)==IntervalType.COALESCENT) {
                coalIntervals[currentIndex] = this.getInterval(i);
                currentIndex+=1;
            }
        }
        return coalIntervals;
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
    private boolean buildIntervalNodeMapping;

    /**
     * The intervals.
     */
    private MutableIntervalList intervals = null;
    private IntervalNodeMapProvider.IntervalNodeMapping intervalNodeMapping;
    /**
     * The stored values for intervals.
     */
    private MutableIntervalList storedIntervals = null;

    private boolean eventsKnown = false;
    private boolean storedEventsKnown = false;

}
