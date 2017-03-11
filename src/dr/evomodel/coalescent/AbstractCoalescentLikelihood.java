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
public abstract class AbstractCoalescentLikelihood extends AbstractModelLikelihood implements Units, CoalescentIntervalProvider {

    // PUBLIC STUFF

    public AbstractCoalescentLikelihood(
            String name,
            Tree tree,
            TaxonList includeSubtree,
            List<TaxonList> excludeSubtrees) throws TreeUtils.MissingTaxonException {

        super(name);

        this.tree = tree;

        if (includeSubtree != null) {
            includedLeafSet = TreeUtils.getLeavesForTaxa(tree, includeSubtree);
        } else {
            includedLeafSet = null;
        }

        if (excludeSubtrees != null) {
            excludedLeafSets = new Set[excludeSubtrees.size()];
            for (int i = 0; i < excludeSubtrees.size(); i++) {
                excludedLeafSets[i] = TreeUtils.getLeavesForTaxa(tree, excludeSubtrees.get(i));
            }
        } else {
            excludedLeafSets = new Set[0];
        }

        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }

        intervals = new Intervals(tree.getNodeCount());
        storedIntervals = new Intervals(tree.getNodeCount());
        eventsKnown = false;

        this.coalescentEventStatisticValues = new double[getNumberOfCoalescentEvents()];

        addStatistic(new DeltaStatistic());

        likelihoodKnown = false;
    }

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    protected final void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == tree) {
            // treeModel has changed so recalculate the intervals
            eventsKnown = false;
        }

        likelihoodKnown = false;
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
    protected final void storeState() {
        // copy the intervals into the storedIntervals
        storedIntervals.copyIntervals(intervals);

        storedEventsKnown = eventsKnown;
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
    }

    /**
     * Restores the precalculated state: that is the intervals of the tree.
     */
    protected final void restoreState() {
        // swap the intervals back
        Intervals tmp = storedIntervals;
        storedIntervals = intervals;
        intervals = tmp;

        eventsKnown = storedEventsKnown;
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
    }

    protected final void acceptState() {
    } // nothing to do

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    public final Model getModel() {
        return this;
    }

    public final double getLogLikelihood() {
        if (!eventsKnown) {
            setupIntervals();
        }

        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }

        return logLikelihood;
    }

    public final void makeDirty() {
        likelihoodKnown = false;
        eventsKnown = false;
    }

    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given a demographic model.
     */
    public abstract double calculateLogLikelihood();

    /**
     * @return the node ref of the MRCA of this coalescent prior in the given tree.
     */
    protected NodeRef getIncludedMRCA(Tree tree) {
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
    protected Set<NodeRef> getExcludedMRCAs(Tree tree) {

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
    protected final void setupIntervals() {

        intervals.resetEvents();
        collectTimes(tree, getIncludedMRCA(tree), getExcludedMRCAs(tree), intervals);
        // force a calculation of the intervals...
        intervals.getIntervalCount();

        eventsKnown = true;
        likelihoodKnown = false;
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

    public double getCoalescentInterval(int i) {
        if (!eventsKnown) {
            setupIntervals();
        }
        return intervals.getInterval(i);
    }

    public int getCoalescentIntervalDimension() {
        if (!eventsKnown) {
            setupIntervals();
        }
        return intervals.getIntervalCount();
    }

    public int getNumberOfCoalescentEvents() {
        return tree.getExternalNodeCount()-1;
    }

    public int getCoalescentIntervalLineageCount(int i) {
        if (!eventsKnown) {
            setupIntervals();
        }
        return intervals.getLineageCount(i);
    }

    public IntervalType getCoalescentIntervalType(int i) {
        if (!eventsKnown) {
            setupIntervals();
        }
        return intervals.getIntervalType(i);
    }

    public double getCoalescentEventsStatisticValue(int i) {
        if (i == 0) {
            for (int j = 0; j < coalescentEventStatisticValues.length; j++) {
                coalescentEventStatisticValues[j] = 0.0;
            }
            int counter = 0;
            for (int j = 0; j < getCoalescentIntervalDimension(); j++) {
                if (getCoalescentIntervalType(j) == IntervalType.COALESCENT) {
                    this.coalescentEventStatisticValues[counter] += getCoalescentInterval(j) * (getCoalescentIntervalLineageCount(j) * (getCoalescentIntervalLineageCount(j) - 1.0)) / 2.0;
                    counter++;
                } else {
                    this.coalescentEventStatisticValues[counter] += getCoalescentInterval(j) * (getCoalescentIntervalLineageCount(j) * (getCoalescentIntervalLineageCount(j) - 1.0)) / 2.0;
                }
            }
        }
        return coalescentEventStatisticValues[i];
    }

    public String toString() {
        return Double.toString(logLikelihood);

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
    private final Set<String> includedLeafSet;
    private final Set[] excludedLeafSets;

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

    private double logLikelihood;
    private double storedLogLikelihood;
    protected boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;

    private double[] coalescentEventStatisticValues;
}