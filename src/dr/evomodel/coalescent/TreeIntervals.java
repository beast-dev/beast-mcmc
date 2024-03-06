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
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;

import java.util.HashSet;
import java.util.List;
import java.util.Set;


/**
 * Forms a base class for a number of coalescent likelihood calculators.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: CoalescentLikelihood.java,v 1.43 2006/07/28 11:27:32 rambaut Exp $
 */
public class TreeIntervals extends AbstractModel implements Units, IntervalList {

    // PUBLIC STUFF

    public TreeIntervals(
            Tree tree,
            TaxonList includeSubtree,
            List<TaxonList> excludeSubtrees) throws TreeUtils.MissingTaxonException {

        super("TreeIntervals");

        this.tree = tree;

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

        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }

        intervals = new Intervals(tree.getNodeCount());
        storedIntervals = new Intervals(tree.getNodeCount());
        eventsKnown = false;

        addStatistic(new DeltaStatistic());
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
    }

    /**
     * Restores the precalculated state: that is the intervals of the tree.
     */
    protected void restoreState() {
        // swap the intervals back
        Intervals tmp = storedIntervals;
        storedIntervals = intervals;
        intervals = tmp;

        eventsKnown = storedEventsKnown;
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
     *         May return null if no subtrees should be excluded.
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
    public final void calculateIntervals() {

        intervals.resetEvents();
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
    private void collectTimes(Tree tree, Intervals intervals) {

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            intervals.addSampleEvent(tree.getNodeHeight(tree.getExternalNode(i)));
        }
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            intervals.addCoalescentEvent(tree.getNodeHeight(tree.getInternalNode(i)));
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
        if(!eventsKnown){
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

    /**
     * The intervals.
     */
    private Intervals intervals = null;

    /**
     * The stored values for intervals.
     */
    private Intervals storedIntervals = null;

    private boolean eventsKnown = false;
    private boolean storedEventsKnown = false;
}
