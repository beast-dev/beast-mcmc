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
        super("TreeIntervals");
        setup(tree);
    }

    public TreeIntervals(
            Tree tree,
            TaxonList includeSubtree,
            List<TaxonList> excludeSubtrees) throws TreeUtils.MissingTaxonException {

        super("TreeIntervals");

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

        sampleCount = tree.getExternalNodeCount();
        coalescentCount = tree.getInternalNodeCount();

//            intervals = new Intervals(tree.getNodeCount());
//            storedIntervals = new Intervals(tree.getNodeCount());
        // unlike the Intervals class, FastIntervals needs the exact number of sample and coalescent
        // events at construction to create fixed size arrays.
        if (excludedLeafSets != null) {
            for (Set excludedLeafSet : excludedLeafSets) {
                // remove all the tip sample events but add one back for
                // where the excluded subtree joins the full tree
                sampleCount = sampleCount - excludedLeafSet.size() + 1;
                // remove all the coalescent events in the excluded subtree
                coalescentCount = coalescentCount - (excludedLeafSet.size() - 1);
            }
        }
        if (includedLeafSet != null) {
            sampleCount = includedLeafSet.size();
            coalescentCount = includedLeafSet.size() - 1;
        }
        intervals = new FastIntervals(sampleCount, coalescentCount);
        storedIntervals = new FastIntervals(sampleCount, coalescentCount);

        eventsKnown = false;

        addStatistic(new DeltaStatistic());
    }

    // This option is set in the constructor as final
//    public void setBuildIntervalNodeMapping(boolean buildIntervalNodeMapping){
//        this.buildIntervalNodeMapping = buildIntervalNodeMapping;
//        this.intervalNodeMapping = buildIntervalNodeMapping ? new IntervalNodeMapping.Default(tree.getNodeCount(),tree):new IntervalNodeMapping.None();
//        //Force a recalculation here
//        eventsKnown = false;
//    }
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

    public boolean isMonophyly() {
        boolean monophyly = true;
        if (includedLeafSet != null) {
            if (!TreeUtils.isMonophyletic(tree, includedLeafSet)) {
                monophyly = false;
            }
        }
        if (excludedLeafSets != null) {
            for (Set<String> leafSet : excludedLeafSets) {
                if (!TreeUtils.isMonophyletic(tree, leafSet)) {
                    monophyly = false;
                }
            }
        }
        return monophyly;
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
     * public method for unit tests
     */
    public void testStoreState() {
        storeState();
    }


    /**
     * Stores the precalculated state: in this case the intervals
     */
    protected void storeState() {
        // copy the intervals into the storedIntervals
        storedIntervals.copyIntervals(intervals);
        storedEventsKnown = eventsKnown;
    }

    /**
     * Restores the precalculated state: that is the intervals of the tree.
     */
    protected void restoreState() {
        // swap the intervals back
        FastIntervals tmp = storedIntervals;
        storedIntervals = intervals;
        intervals = tmp;

        eventsKnown = storedEventsKnown;

        assert isMonophyly();
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
            NodeRef mrca = TreeUtils.getCommonAncestorNode(tree, includedLeafSet);
            return mrca;
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

        assert monophyly;

        intervals.resetEvents();

        includedCoalescentCount = 0;

        if (includedLeafSet != null || excludedLeafSets != null) {
            collectTimes(tree, getIncludedMRCA(tree), getExcludedMRCAs(tree), intervals);
        } else {
            collectTimes(tree, intervals);
        }

        // force a calculation of the intervals...
        intervals.getIntervalCount();

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

        includedCoalescentCount += 1;
        intervals.addCoalescentEvent(tree.getNodeHeight(node));

        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef child = tree.getChild(node, i);

            // check if this subtree is included in the coalescent density
            boolean include = excludeNodesBelow == null || !excludeNodesBelow.contains(child);

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

    public boolean hasIntervalChanged(int i) {
        if (!eventsKnown) {
            calculateIntervals();
        }
        boolean intervalChanged = intervals.getIntervalEndTime(i) != storedIntervals.getIntervalEndTime(i) ||
                intervals.getIntervalType(i) != storedIntervals.getIntervalType(i) ||
                intervals.getLineageCount(i) != storedIntervals.getLineageCount(i) ||
                intervals.getIntervalNodeNumber(i) != storedIntervals.getIntervalNodeNumber(i);

        return intervalChanged;
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
        return getIntervalEndTime(i);
    }

    @Override
    public double getIntervalStartTime(int i) {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.getIntervalStartTime(i);
    }

    @Override
    public double getIntervalEndTime(int i) {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.getIntervalEndTime(i);
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
    public int getCoalescentEventCount(int i) {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.getCoalescentEventCount(i);
    }

    @Override
    public IntervalType getIntervalType(int i) {
        return getIntervalEndType(i);
    }

    @Override
    public IntervalType getIntervalStartType(int i) {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.getIntervalStartType(i);
    }

    @Override
    public IntervalType getIntervalEndType(int i) {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.getIntervalEndType(i);
    }

    @Override
    public int getIntervalNodeNumber(int interval) {
        return intervals.getIntervalNodeNumber(interval);
    }

    @Override
    public int getIntervalStartNodeNumber(int interval) {
        return intervals.getIntervalStartNodeNumber(interval);
    }

    @Override
    public int getIntervalEndNodeNumber(int interval) {
        return intervals.getIntervalEndNodeNumber(interval);
    }

    @Override
    public NodeRef getIntervalNode(int interval) {
        return tree.getNode(intervals.getIntervalStartNodeNumber(interval));
    }

    @Override
    public NodeRef getIntervalStartNode(int interval) {
        return tree.getNode(intervals.getIntervalStartNodeNumber(interval));
    }

    @Override
    public NodeRef getIntervalEndNode(int interval) {
        return tree.getNode(intervals.getIntervalStartNodeNumber(interval));
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

//    @Override
//    public double[] getCoalescentIntervals() {
//        double[] coalIntervals = new double[tree.getInternalNodeCount()];
//        int currentIndex = 0;
//        for (int i = 0; i < this.intervals.getIntervalCount(); i++) {
//            if (this.getIntervalEndType(i)==IntervalType.COALESCENT) {
//                coalIntervals[currentIndex] = this.getInterval(i);
//                currentIndex+=1;
//            }
//        }
//        return coalIntervals;
//    }

    @Override
    public Type getUnits() {
        return intervals.getUnits();
    }

    @Override
    public void setUnits(Type units) {
        intervals.setUnits(units);
    }

    public void setIntervalsUnknown() {
        eventsKnown = false;
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
    private Set<String>[] excludedLeafSets = null;

    private int sampleCount;
    private int coalescentCount;

    private int includedCoalescentCount;

    private boolean monophyly = true;

    /**
     * The intervals.
     */
    private FastIntervals intervals = null;
    /**
     * The stored values for intervals.
     */
    private FastIntervals storedIntervals = null;

    private boolean eventsKnown = false;
    private boolean storedEventsKnown = false;
}
