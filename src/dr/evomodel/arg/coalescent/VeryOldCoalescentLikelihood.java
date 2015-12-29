/*
 * VeryOldCoalescentLikelihood.java
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

package dr.evomodel.arg.coalescent;

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.DemographicModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.math.Binomial;
import dr.util.ComparableDouble;
import dr.util.HeapSort;
import dr.xml.*;

import java.util.ArrayList;


/**
 * A likelihood function for the coalescent. Takes a tree and a demographic model.
 * <p/>
 * Parts of this class were derived from C++ code provided by Oliver Pybus.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: VeryOldCoalescentLikelihood.java,v 1.43 2006/07/28 11:27:32 rambaut Exp $
 */
public class VeryOldCoalescentLikelihood extends AbstractModelLikelihood implements Units {

    // PUBLIC STUFF

    public static final String COALESCENT_LIKELIHOOD = "veryOldCoalescentLikelihood";
    public static final String ANALYTICAL = "analytical";
    public static final String MODEL = "model";
    public static final String POPULATION_TREE = "populationTree";

    /**
     * Denotes an interval after which a coalescent event is observed
     * (i.e. the number of lineages is smaller in the next interval)
     */
    public static final int COALESCENT = 0;

    /**
     * Denotes an interval at the end of which a new sample addition is
     * observed (i.e. the number of lineages is larger in the next interval).
     */
    public static final int NEW_SAMPLE = 1;

    /**
     * Denotes an interval at the end of which nothing is
     * observed (i.e. the number of lineages is the same in the next interval).
     */
    public static final int NOTHING = 2;


    public VeryOldCoalescentLikelihood(Tree tree, DemographicModel demoModel) {
        this(COALESCENT_LIKELIHOOD, tree, demoModel, true);
    }

    public VeryOldCoalescentLikelihood(String name, Tree tree, DemographicModel demoModel, boolean setupIntervals) {

        super(name);

        this.tree = tree;
        this.demoModel = demoModel;
        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }
        if (demoModel != null) {
            addModel(demoModel);
        }
        if (setupIntervals) setupIntervals();

        addStatistic(new DeltaStatistic());
    }

    VeryOldCoalescentLikelihood(String name) {
        super(name);
    }

    // **************************************************************
    // Extendable methods
    // **************************************************************

    /**
     * @return the node ref of the MRCA of this coalescent prior in the given tree.
     */
    public NodeRef getMRCAOfCoalescent(Tree tree) {
        return tree.getRoot();
    }

    /**
     * @return an array of noderefs that represent the MRCAs of subtrees to exclude from coalescent prior.
     *         May return null if no subtrees should be excluded.
     */
    public NodeRef[] getExcludedMRCAs(Tree tree) {
        return null;
    }

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == tree) {
            // treeModel has changed so recalculate the intervals
            intervalsKnown = false;
        } else {
            // demoModel has changed so we don't need to recalculate the intervals
        }

        likelihoodKnown = false;
    }

    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
    } // No parameters to respond to

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    /**
     * Stores the precalculated state: in this case the intervals
     */
    protected void storeState() {
        System.arraycopy(intervals, 0, storedIntervals, 0, intervals.length);
        System.arraycopy(lineageCounts, 0, storedLineageCounts, 0, lineageCounts.length);
        storedIntervalsKnown = intervalsKnown;
        storedIntervalCount = intervalCount;
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
    }

    /**
     * Restores the precalculated state: that is the intervals of the tree.
     */
    protected void restoreState() {
        System.arraycopy(storedIntervals, 0, intervals, 0, storedIntervals.length);
        System.arraycopy(storedLineageCounts, 0, lineageCounts, 0, storedLineageCounts.length);
        intervalsKnown = storedIntervalsKnown;
        intervalCount = storedIntervalCount;
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;

        if (!intervalsKnown) {
            likelihoodKnown = false;
        }
    }

    protected final void acceptState() {
    } // nothing to do

    /**
     * Adopt the state of the model from source.
     */
    protected final void adoptState(Model source) {
        // all we need to do is force a recalculation of intervals
        makeDirty();
    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    public final Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    public final void makeDirty() {
        likelihoodKnown = false;
        intervalsKnown = false;
    }

    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given a demographic model.
     */
    public double calculateLogLikelihood() {

        if ( !intervalsKnown ) setupIntervals();

        if (demoModel == null) return calculateAnalyticalLogLikelihood();

        double logL = 0.0;

        double currentTime = 0.0;

        DemographicFunction demoFunction = demoModel.getDemographicFunction();

        for (int j = 0; j < intervalCount; j++) {


            logL += calculateIntervalLikelihood(demoFunction, intervals[j], currentTime, lineageCounts[j],
                    getIntervalType(j));

            // insert zero-length coalescent intervals
            int diff = getCoalescentEvents(j) - 1;
            for (int k = 0; k < diff; k++) {
                logL += calculateIntervalLikelihood(demoFunction, 0.0, currentTime, lineageCounts[j] - k - 1, COALESCENT);
            }

            currentTime += intervals[j];
        }

        return logL;
    }

    private double calculateAnalyticalLogLikelihood() {

        final double lambda = getLambda();
        final int n = tree.getExternalNodeCount();

        // assumes a 1/theta prior
        //logLikelihood = Math.log(1.0/Math.pow(lambda,n));

        // assumes a flat prior
        double logL = Math.log(1.0 / Math.pow(lambda, n - 1));
        return logL;
    }

    /**
     * Returns the likelihood of a given *coalescent* interval
     */
    public final double calculateIntervalLikelihood(DemographicFunction demoFunction, double width, double timeOfPrevCoal, int lineageCount) {

        return calculateIntervalLikelihood(demoFunction, width, timeOfPrevCoal, lineageCount, COALESCENT);
    }

    /**
     * Returns the likelihood of a given interval,coalescent or otherwise.
     */
    public final double calculateIntervalLikelihood(DemographicFunction demoFunction, double width, double timeOfPrevCoal,
                                                    int lineageCount, int type) {
        //binom.setMax(lineageCount);

        double timeOfThisCoal = width + timeOfPrevCoal;

//        System.err.printf("s: %7.6f   f: %7.6f,  %d, %d\n", timeOfPrevCoal, timeOfThisCoal, lineageCount, type);
        double intervalArea = demoFunction.getIntegral(timeOfPrevCoal, timeOfThisCoal);
        double like = 0;
        switch (type) {
            case COALESCENT:
                like =
                        -Math.log(demoFunction.getDemographic(timeOfThisCoal)) -
                                (Binomial.choose2(lineageCount) * intervalArea);
                break;
            case NEW_SAMPLE:
                like = -(Binomial.choose2(lineageCount) * intervalArea);
                break;
        }

        return like;
    }

    /**
     * Returns a factor lambda such that the likelihood can be expressed as
     * 1/theta^(n-1) * exp(-lambda/theta). This allows theta to be integrated
     * out analytically. :-)
     */
    private double getLambda() {
        double lambda = 0.0;
        for (int i = 0; i < getIntervalCount(); i++) {
            lambda += (intervals[i] * lineageCounts[i]);
        }
        lambda /= 2;

        return lambda;
    }

    /**
     * Recalculates all the intervals from the tree model.
     */
    protected final void setupIntervals() {

        double MULTIFURCATION_LIMIT = 1e-9;

        ArrayList times = new ArrayList();
        ArrayList<Integer> childs = new ArrayList<Integer>();
        collectAllTimes(tree, getMRCAOfCoalescent(tree), getExcludedMRCAs(tree), times, childs);
        int[] indices = new int[times.size()];

        HeapSort.sort(times, indices);

        int maxIntervalCount = tree.getNodeCount();

        if (intervals == null) {
            intervals = new double[maxIntervalCount];
            lineageCounts = new int[maxIntervalCount];
            storedIntervals = new double[maxIntervalCount];
            storedLineageCounts = new int[maxIntervalCount];
        }

        // start is the time of the first tip
        double start = ((ComparableDouble) times.get(indices[0])).doubleValue();
        int numLines = 0;
        int i = 0;
        intervalCount = 0;
        while (i < times.size()) {

            int lineagesRemoved = 0;
            int lineagesAdded = 0;

            double finish = ((ComparableDouble) times.get(indices[i])).doubleValue();
            double next = finish;

            while (Math.abs(next - finish) < MULTIFURCATION_LIMIT) {
                int children = childs.get(indices[i]);
                if (children == 0) {
                    lineagesAdded += 1;
                } else {
                    lineagesRemoved += (children - 1);
                }
                i += 1;
                if (i < times.size()) {
                    next = ((ComparableDouble) times.get(indices[i])).doubleValue();
                } else break;
            }
            //System.out.println("time = " + finish + " removed = " + lineagesRemoved + " added = " + lineagesAdded);
            if (lineagesAdded > 0) {

                if (intervalCount > 0 || ((finish - start) > MULTIFURCATION_LIMIT)) {
                    intervals[intervalCount] = finish - start;
                    lineageCounts[intervalCount] = numLines;
                    intervalCount += 1;
                }

                start = finish;
            }
            // add sample event
            numLines += lineagesAdded;

            if (lineagesRemoved > 0) {

                intervals[intervalCount] = finish - start;
                lineageCounts[intervalCount] = numLines;
                intervalCount += 1;
                start = finish;
            }
            // coalescent event
            numLines -= lineagesRemoved;
        }
        intervalsKnown = true;
    }


    /**
     * extract coalescent times and tip information into ArrayList times from tree.
     *
     * @param node         the node to start from
     * @param excludeBelow an optional array of nodes to exclude (corresponding subtrees) from density.
     */
    private static void collectAllTimes(Tree tree, NodeRef node, NodeRef[] excludeBelow, ArrayList times, ArrayList<Integer> childs) {

        times.add(new ComparableDouble(tree.getNodeHeight(node)));
        childs.add(tree.getChildCount(node));

        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef child = tree.getChild(node, i);
            if (excludeBelow == null) {
                collectAllTimes(tree, child, excludeBelow, times, childs);
            } else {
                // check if this subtree is included in the coalescent density
                boolean include = true;
                for(NodeRef anExcludeBelow : excludeBelow) {
                    if( anExcludeBelow.getNumber() == child.getNumber() ) {
                        include = false;
                        break;
                    }
                }
                if (include) collectAllTimes(tree, child, excludeBelow, times, childs);
            }
        }
    }


    /**
     * get number of intervals
     */
    public final int getIntervalCount() {
        return intervalCount;
    }

    /**
     * Gets an interval.
     */
    public final double getInterval(int i) {
        if (i >= intervalCount) throw new IllegalArgumentException();
        return intervals[i];
    }

    /**
     * Returns the number of uncoalesced lineages within this interval.
     * Required for s-coalescents, where new lineages are added as
     * earlier samples are come across.
     */
    public final int getLineageCount(int i) {
        if (i >= intervalCount) throw new IllegalArgumentException();
        return lineageCounts[i];
    }

    /**
     * Returns the number coalescent events in an interval
     */
    public final int getCoalescentEvents(int i) {

        if (i >= intervalCount) throw new IllegalArgumentException();
        if (i < intervalCount - 1) {
            return lineageCounts[i] - lineageCounts[i + 1];
        } else {
            return lineageCounts[i] - 1;
        }
    }

    /**
     * Returns the type of interval observed.
     */
    public final int getIntervalType(int i) {

        if (i >= intervalCount) throw new IllegalArgumentException();
        int numEvents = getCoalescentEvents(i);

        if (numEvents > 0) return COALESCENT;
        else if (numEvents < 0) return NEW_SAMPLE;
        else return NOTHING;
    }

    /**
     * get the total height of the genealogy represented by these
     * intervals.
     */
    public final double getTotalHeight() {

        double height = 0.0;
        for (int j = 0; j < intervalCount; j++) {
            height += intervals[j];
        }
        return height;
    }

    /**
     * Checks whether this set of coalescent intervals is fully resolved
     * (i.e. whether is has exactly one coalescent event in each
     * subsequent interval)
     */
    public final boolean isBinaryCoalescent() {
        for (int i = 0; i < intervalCount; i++) {
            if (getCoalescentEvents(i) != 1) return false;
        }

        return true;
    }

    /**
     * Checks whether this set of coalescent intervals coalescent only
     * (i.e. whether is has exactly one or more coalescent event in each
     * subsequent interval)
     */
    public final boolean isCoalescentOnly() {
        for (int i = 0; i < intervalCount; i++) {
            if (getCoalescentEvents(i) < 1) return false;
        }
        return true;
    }

    public String toString() {
        return Double.toString(getLogLikelihood());

    }

    // **************************************************************
    // Units IMPLEMENTATION
    // **************************************************************


    public final Type getUnits() {
        return demoModel.getUnits();
    }

    public void setUnits(Type units) {
        demoModel.setUnits(units);

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

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return COALESCENT_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) {

            XMLObject cxo = xo.getChild(MODEL);
            DemographicModel demoModel = (DemographicModel) cxo.getChild(DemographicModel.class);

            cxo = xo.getChild(POPULATION_TREE);
            TreeModel treeModel = (TreeModel) cxo.getChild(TreeModel.class);

            return new VeryOldCoalescentLikelihood(treeModel, demoModel);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents the likelihood of the tree given the demographic function.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(MODEL, new XMLSyntaxRule[]{
                        new ElementRule(DemographicModel.class)
                }),
                new ElementRule(POPULATION_TREE, new XMLSyntaxRule[]{
                        new ElementRule(TreeModel.class)
                }),
        };
    };

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

/*	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return COALESCENT_LIKELIHOOD; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			DemographicModel demoModel = null;
			if (xo.hasAttribute(MODEL)) {
				demoModel = (DemographicModel)xo.getAttribute(MODEL);
			}
			TreeModel treeModel = (TreeModel)xo.getAttribute(TREE);
			return new VeryOldCoalescentLikelihood(treeModel, demoModel);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element represents the likelihood of the tree given the demographic function.";
		}

		public Class getReturnType() { return Likelihood.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new XORRule(
				new EnumAttributeRule(ANALYTICAL, new String[] { "constant" }),
				new AttributeRule(MODEL, DemographicModel.class)
			),
			new AttributeRule(TREE, TreeModel.class)
		};
	};*/

    /**
     * The demographic model.
     */
    DemographicModel demoModel = null;

    /**
     * The tree.
     */
    Tree tree = null;

    /**
     * The widths of the intervals.
     */
    double[] intervals;
    private double[] storedIntervals;

    /**
     * The number of uncoalesced lineages within a particular interval.
     */
    int[] lineageCounts;
    private int[] storedLineageCounts;

    boolean intervalsKnown = false;
    protected boolean storedIntervalsKnown = false;

    double logLikelihood;
    protected double storedLogLikelihood;
    boolean likelihoodKnown = false;
    protected boolean storedLikelihoodKnown = false;

    int intervalCount = 0;
    private int storedIntervalCount = 0;
}