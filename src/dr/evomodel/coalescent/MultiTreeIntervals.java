/*
 * TreeIntervals.java
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
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Extracts the intervals from a set of trees over a fixed time period
 *
 * @author Andrew Rambaut
 */
public class MultiTreeIntervals extends AbstractModel implements IntervalList {

    public MultiTreeIntervals(Tree tree, boolean includeStems, double cutoffTime) {
        super("MultiTreeIntervals");
        throw new UnsupportedOperationException("Not implemented yet");
    }


    public MultiTreeIntervals(Collection<Tree> trees, Taxa singletonTaxa, boolean includeStems, double cutoffTime) {
        super("MultiTreeIntervals");

        int maxEventCount = 0;
        for (Tree tree : trees) {
            if (tree instanceof Model) {
                addModel((Model) tree);
            }
            // one event for each tip, internal node and one extra for the top of the stem.
            maxEventCount += tree.getNodeCount() + 1;
        }
        if (includeStems && singletonTaxa != null) {
            // two events each for the singletons (one at the top and bottom).
            maxEventCount += singletonTaxa.getTaxonCount() * 2;
        }


        this.trees = new ArrayList<Tree>(trees);
        this.singletonTaxa = singletonTaxa;
        this.includeStems = includeStems;
        this.cutoffTime = cutoffTime;
        this.units = this.trees.get(0).getUnits();

        this.intervals = new Intervals(maxEventCount);
        this.storedIntervals = new Intervals(maxEventCount);
    }


    @Override
    public double getStartTime() {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.getStartTime();
    }

    @Override
    public int getSampleCount() {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.getSampleCount();
    }

    /**
     * get number of intervals
     */
    @Override
    public int getIntervalCount() {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.getIntervalCount();
    }

    /**
     * Gets an interval.
     */
    @Override
    public double getInterval(int i) {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.getInterval(i);
    }

    /**
     * Returns the time of the start of an interval
     */
    @Override
    public double getIntervalTime(int i) {
        if(!eventsKnown){
            calculateIntervals();
        }
        return intervals.getIntervalTime(i);
    }

    /**
     * Returns the number of uncoalesced lineages within this interval.
     * Required for s-coalescents, where new lineages are added as
     * earlier samples are come across.
     */
    @Override
    public int getLineageCount(int i) {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.getLineageCount(i);
    }

    /**
     * Returns the number coalescent events in an interval
     */
    @Override
    public int getCoalescentEvents(int i) {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.getCoalescentEvents(i);
    }

    /**
     * Returns the type of interval observed.
     */
    @Override
    public IntervalType getIntervalType(int i) {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.getIntervalType(i);
    }

    /**
     * get the total height of the genealogy represented by these
     * intervals.
     */
    @Override
    public double getTotalDuration() {

        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.getTotalDuration();
    }

    /**
     * Checks whether this set of coalescent intervals is fully resolved
     * (i.e. whether is has exactly one coalescent event in each
     * subsequent interval)
     */
    @Override
    public boolean isBinaryCoalescent() {
        return intervals.isBinaryCoalescent();
    }

    /**
     * Checks whether this set of coalescent intervals coalescent only
     * (i.e. whether is has exactly one or more coalescent event in each
     * subsequent interval)
     */
    @Override
    public boolean isCoalescentOnly() {
        if (!eventsKnown) {
            calculateIntervals();
        }
        return intervals.isCoalescentOnly();
    }

    /**
     * Recalculates all the intervals for the given tree.
     */
    public void calculateIntervals() {

        intervals.resetEvents();
        for (Tree tree : trees) {
            for (int i = 0; i < tree.getExternalNodeCount(); i++) {
                NodeRef node = tree.getExternalNode(i);
                intervals.addSampleEvent(tree.getNodeTaxon(node).getHeight());
            }
            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                NodeRef node = tree.getInternalNode(i);
                intervals.addCoalescentEvent(tree.getNodeHeight(node));
            }
            if (includeStems) {
                // add the nothing event at the top of the stem of the root of each subtree.
                intervals.addNothingEvent(cutoffTime);
            } else {
                // add nothing event at the root of the subtree
                NodeRef node = tree.getInternalNode(tree.getInternalNodeCount() - 1);
                intervals.addNothingEvent(tree.getNodeHeight(node));
            }
        }
        if(includeStems && singletonTaxa != null){
            for (Taxon taxon : singletonTaxa) {
                intervals.addSampleEvent(taxon.getHeight());
                intervals.addNothingEvent(cutoffTime);
            }
        }


        // call this to sort and calculate the intervals in the inner
        // object.
        intervals.getIntervalCount();

        fireModelChanged();

        eventsKnown = true;
    }

    /**
     * Return the units that this tree is expressed in.
     */
    public final Type getUnits() {
        return units;
    }

    /**
     * Sets the units that this tree is expressed in.
     */
    public final void setUnits(Type units) {
        throw new IllegalArgumentException("Can't set interval's units");
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        eventsKnown = false;
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // nothing to do.
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
        Intervals tmp = storedIntervals;
        storedIntervals = intervals;
        intervals = tmp;

        eventsKnown = storedEventsKnown;
    }

    @Override
    protected void acceptState() {
        // do nothing
    }

    /**
     * The tree.
     */
    private final List<Tree> trees;
    private final Taxa singletonTaxa;
    private final boolean includeStems;
    private final double cutoffTime;
    private final Type units;

    /**
     * The intervals.
     */
    private Intervals intervals = null;

    /**
     * The stored values for intervals.
     */
    private Intervals storedIntervals = null;

    /**
     * are the events known?
     */
    private boolean eventsKnown = false;
    private boolean storedEventsKnown = false;

}