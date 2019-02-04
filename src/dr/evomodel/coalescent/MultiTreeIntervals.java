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
import dr.util.HeapSort;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * Extracts the intervals from a set of trees over a fixed time period
 *
 * @author Andrew Rambaut
 */
public class MultiTreeIntervals extends AbstractModel implements IntervalList {

    public MultiTreeIntervals(Tree tree, double cutoffTime) {
        super("MultiTreeIntervals");
        throw new UnsupportedOperationException("Not implemented yet");
    }


    public MultiTreeIntervals(Collection<Tree> trees, Taxa singletonTaxa, double cuttoffTime) {
        super("MultiTreeIntervals");

        int maxEventCount = 0;
        for (Tree tree : trees) {
            if (tree instanceof Model) {
                addModel((Model) tree);
            }
            // one event for each tip, internal node and one extra for the top of the stem.
            maxEventCount += tree.getNodeCount() + 1;
        }
        // two events each for the singletons (one at the top and bottom).
        maxEventCount += singletonTaxa.getTaxonCount() * 2;

        this.trees = new ArrayList<Tree>(trees);
        this.singletonTaxa = singletonTaxa;
        this.cutoffTime = cuttoffTime;
        this.units = this.trees.get(0).getUnits();

        this.intervals = new Intervals(maxEventCount);
    }


    /**
     * Specifies that the intervals are unknown (i.e., the tree has changed).
     */
    public void setIntervalsUnknown() {
        intervalsKnown = false;
    }

    public int getSampleCount() {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return intervals.getSampleCount();
    }

    /**
     * get number of intervals
     */
    public int getIntervalCount() {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return intervals.getIntervalCount();
    }

    /**
     * Gets an interval.
     */
    public double getInterval(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return intervals.getInterval(i);
    }

    /**
     * Returns the number of uncoalesced lineages within this interval.
     * Required for s-coalescents, where new lineages are added as
     * earlier samples are come across.
     */
    public int getLineageCount(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return intervals.getLineageCount(i);
    }

    /**
     * Returns the number coalescent events in an interval
     */
    public int getCoalescentEvents(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return intervals.getCoalescentEvents(i);
    }

    /**
     * Returns the type of interval observed.
     */
    public IntervalType getIntervalType(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return intervals.getIntervalType(i);
    }

    /**
     * get the total height of the genealogy represented by these
     * intervals.
     */
    public double getTotalDuration() {

        if (!intervalsKnown) {
            calculateIntervals();
        }
        return intervals.getTotalDuration();
    }

    /**
     * Checks whether this set of coalescent intervals is fully resolved
     * (i.e. whether is has exactly one coalescent event in each
     * subsequent interval)
     */
    public boolean isBinaryCoalescent() {
        return intervals.isBinaryCoalescent();
    }

    /**
     * Checks whether this set of coalescent intervals coalescent only
     * (i.e. whether is has exactly one or more coalescent event in each
     * subsequent interval)
     */
    public boolean isCoalescentOnly() {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return intervals.isCoalescentOnly();
    }

    /**
     * Recalculates all the intervals for the given tree.
     */
    private void calculateIntervals() {

        intervals.resetEvents();
        for (Tree tree : trees) {
            for (int i = 0; i < tree.getExternalNodeCount(); i++) {
                NodeRef node = tree.getExternalNode(i);
                intervals.addSampleEvent(tree.getNodeTaxon(node).getHeight());
            }
            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                NodeRef node = tree.getInternalNode(i);
                intervals.addCoalescentEvent(tree.getNodeTaxon(node).getHeight());
            }
            // add the nothing event at the top of the stem of the root of each subtree.
            intervals.addNothingEvent(cutoffTime);
        }
        for (Taxon taxon : singletonTaxa) {
            intervals.addSampleEvent(taxon.getHeight());
            intervals.addNothingEvent(cutoffTime);
        }

        fireModelChanged();

        intervalsKnown = true;
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
        intervalsKnown = false;
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // nothing to do.
    }

    /**
     * Extra functionality to store and restore values for caching
     */
    public void storeState() {
    }

    public void restoreState() {
        intervalsKnown = false;
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
    private final double cutoffTime;
    private final Type units;

    private final Intervals intervals;

    /**
     * are the intervals known?
     */
    private boolean intervalsKnown = false;

}