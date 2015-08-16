/*
 * OldAbstractCoalescentLikelihood.java
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

import dr.evolution.coalescent.Coalescent;
import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.ScaledDemographic;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.CoalescentLikelihoodParser;
import dr.inference.model.*;
import dr.math.Binomial;
import dr.util.ComparableDouble;
import dr.util.HeapSort;

import java.util.ArrayList;

/**
 * Forms a base class for a number of coalescent likelihood calculators.
 * <p/>
 * This was the former 'CoalescentLikelihood' which, as of BEAST v1.4.x was replaced as the
 * standard coalescent likelihood by 'NewCoalescentLikelihood'. As this class is used as a base
 * by a number of other classes (i.e., BayesianSkylineLikelihood), I have made this class abstract
 * and removed the parser.
 * <p/>
 * NewCoalescentLikelihood is now CoalesecentLikelihood (it's parser was installed as the default
 * 'coalescentLikelihood' anyway).
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: CoalescentLikelihood.java,v 1.43 2006/07/28 11:27:32 rambaut Exp $
 */
public class OldAbstractCoalescentLikelihood extends AbstractModelLikelihood implements  Units {

    // PUBLIC STUFF

    //public static final String COALESCENT_LIKELIHOOD = "oldcoalescentLikelihood";
    //    public static final String ANALYTICAL = "analytical";
    //    public static final String MODEL = "model";
    //
    //    public static final String POPULATION_TREE = "populationTree";
    //    public static final String POPULATION_FACTOR = "factor";
    protected MultiLociTreeSet treesSet = null;

    public enum CoalescentEventType {
        /**
         * Denotes an interval after which a coalescent event is observed
         * (i.e. the number of lineages is smaller in the next interval)
         */
        COALESCENT,
        /**
         * Denotes an interval at the end of which a new sample addition is
         * observed (i.e. the number of lineages is larger in the next interval).
         */
        NEW_SAMPLE,
        /**
         * Denotes an interval at the end of which nothing is
         * observed (i.e. the number of lineages is the same in the next interval).
         */
        NOTHING
    }

    public OldAbstractCoalescentLikelihood(Tree tree, DemographicModel demoModel) {
        this(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, tree, demoModel, true);
    }

    public OldAbstractCoalescentLikelihood(MultiLociTreeSet treesSet, DemographicModel demoModel) {
        super(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD);
        this.demoModel = demoModel;
        this.tree = null;
        this.treesSet = treesSet;

        if (demoModel != null) {
            addModel(demoModel);
        }

        for (int nt = 0; nt < treesSet.nLoci(); ++nt) {
            final Tree t = treesSet.getTree(nt);
            if (t instanceof Model) {
                addModel((Model) t);
            }
        }
    }

    public OldAbstractCoalescentLikelihood(String name, Tree tree, DemographicModel demoModel, boolean setupIntervals) {
        super(name);

        this.demoModel = demoModel;
        this.tree = tree;

        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }

        if (demoModel != null) {
            addModel(demoModel);
        }

        if (setupIntervals) setupIntervals();

        addStatistic(new DeltaStatistic());
    }

    OldAbstractCoalescentLikelihood(String name) {
        super(name);
    }

    // **************************************************************
    // Extendable methods
    // **************************************************************

    /**
     * @param tree given tree
     * @return the node ref of the MRCA of this coalescent prior in the given tree (i.e. root of tree)
     */
    public NodeRef getMRCAOfCoalescent(Tree tree) {
        return tree.getRoot();
    }

    /**
     * @param tree given tree
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

    // No parameters to respond to

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    /**
     * Stores the precalculated state: in this case the intervals
     */
    protected void storeState() {
        if (tree != null) {
            System.arraycopy(intervals, 0, storedIntervals, 0, intervals.length);
            System.arraycopy(lineageCounts, 0, storedLineageCounts, 0, lineageCounts.length);
            storedIntervalsKnown = intervalsKnown;
            storedIntervalCount = intervalCount;
            storedLikelihoodKnown = likelihoodKnown;
        } else if (treesSet != null) {
            treesSet.storeTheState();
        }
        storedLogLikelihood = logLikelihood;
    }

    /**
     * Restores the precalculated state: that is the intervals of the tree.
     */
    protected void restoreState() {
        if (tree != null) {
            System.arraycopy(storedIntervals, 0, intervals, 0, storedIntervals.length);
            System.arraycopy(storedLineageCounts, 0, lineageCounts, 0, storedLineageCounts.length);
            intervalsKnown = storedIntervalsKnown;
            intervalCount = storedIntervalCount;
        } else if (treesSet != null) {
            treesSet.restoreTheState();
        }

        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;

        if (!intervalsKnown) {
            likelihoodKnown = false;
        }
    }

    protected final void acceptState() {
    } // nothing to do

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
     * @return the log likelihood of this set of coalescent intervals,
     *         given a demographic model
     */
    public double calculateLogLikelihood() {

        if (treesSet != null) {
            final int nTrees = treesSet.nLoci();
            final DemographicFunction demogFunction = demoModel.getDemographicFunction();
            double logLike = 0.0;
            for (int nt = 0; nt < nTrees; ++nt) {
                final double popFactor = treesSet.getPopulationFactor(nt);
                DemographicFunction df = popFactor != 1.0 ?
                        new ScaledDemographic(demogFunction, popFactor) : demogFunction;

                logLike += Coalescent.calculateLogLikelihood(treesSet.getTreeIntervals(nt), df);
            }
            return logLike;
        }

        if (!intervalsKnown) setupIntervals();

        if (demoModel == null) return calculateAnalyticalLogLikelihood();

        double logL = 0.0;

        double currentTime = 0.0;

        DemographicFunction demoFunction = demoModel.getDemographicFunction();

        for (int j = 0; j < intervalCount; j++) {

            logL += calculateIntervalLikelihood(demoFunction, intervals[j], currentTime, lineageCounts[j],
                    getIntervalType(j));

            // insert zero-length coalescent intervals
            final int diff = getCoalescentEvents(j) - 1;
            for (int k = 0; k < diff; k++) {
                logL += calculateIntervalLikelihood(demoFunction, 0.0, currentTime, lineageCounts[j] - k - 1,
                        CoalescentEventType.COALESCENT);
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
        //double logL = Math.log(1.0/Math.pow(lambda,n-1));
        //final double logL = - Math.log(Math.pow(lambda,n-1));
        return -(n - 1) * Math.log(lambda);
    }

    /**
     * Returns the likelihood of a given *coalescent* interval
     */
    public final double calculateIntervalLikelihood(DemographicFunction demoFunction, double width,
                                                    double timeOfPrevCoal, int lineageCount) {

        return calculateIntervalLikelihood(demoFunction, width, timeOfPrevCoal, lineageCount,
                CoalescentEventType.COALESCENT);
    }

    /**
     * k - number of lineages
     * N - population size
     * kingsman coalescent: interval to next coalescent event x ~ exp(lambda), where lambda = C(k,2) / N
     * Like(x ; lambda) = lambda * exp(-lambda * x)
     * so Like(N) = (C(k,2)/N) * exp(- x * C(k,2)/N)
     * lg(Like(N)) = lg(C(k,2)) - lg(N) -C(k,2) * x/N
     * <p/>
     * When N changes over time N = N(t) we have lambda(t) = C(k,2)/N(t) and the likelihood equation is
     * Like(t) = lambda(t) * exp(- integral_0^t(lambda(x) dx) )
     * <p/>
     * lg(Like(t)) = -C(k,2) * integral_0^t(1/N(x) dx) + lg(C(k,2)/N(t))
     * <p/>
     * For a sample event, the likelihood is for no event until time t, and is just the first term of the above.
     *
     * @param demogFunction  the demographic function
     * @param width          the size of the coalescent interval
     * @param timeOfPrevCoal the time of previous coalescent event (going backwards in time)
     * @param lineageCount   the number of lineages spanning this coalescent interval
     * @param type           the type of coalescent event that this interval is terminated by
     * @return likelihood of a given interval,coalescent or otherwise
     */
    public static double calculateIntervalLikelihood(DemographicFunction demogFunction,
                                                     double width, double timeOfPrevCoal, int lineageCount,
                                                     CoalescentEventType type) {
        final double timeOfThisCoal = width + timeOfPrevCoal;

        final double intervalArea = demogFunction.getIntegral(timeOfPrevCoal, timeOfThisCoal);
        final double kchoose2 = Binomial.choose2(lineageCount);
        double like = -kchoose2 * intervalArea;

        switch (type) {
            case COALESCENT:
                final double demographic = demogFunction.getLogDemographic(timeOfThisCoal);
                like += -demographic;

                break;
            case NEW_SAMPLE:
                break;
        }

        return like;
    }

    /**
     * @return the exponent of the population size parameter (shape parameter of Gamma distribution) associated to the
     *         likelihood for this interval, coalescent or otherwise
     */
    public final double calculateIntervalShapeParameter(DemographicFunction demogFunction,
                                                        double width, double timeOfPrevCoal, int lineageCount, CoalescentEventType type) {
        switch (type) {
            case COALESCENT:
                return 1.0;
            case NEW_SAMPLE:
                return 0.0;
        }
        throw new Error("Unknown event found");
    }

    /**
     * @return the intensity of coalescences (rate parameter, or inverse scale parameter, of Gamma distribution)
     *         associated to the likelihood for this interval, coalescent or otherwise
     */
    public final double calculateIntervalRateParameter(DemographicFunction demogFunction, double width,
                                                       double timeOfPrevCoal, int lineageCount, CoalescentEventType type) {
        final double timeOfThisCoal = width + timeOfPrevCoal;
        final double intervalArea = demogFunction.getIntegral(timeOfPrevCoal, timeOfThisCoal);
        return Binomial.choose2(lineageCount) * intervalArea;
    }


    /**
     * @return a factor lambda such that the likelihood can be expressed as
     *         1/theta^(n-1) * exp(-lambda/theta). This allows theta to be integrated
     *         out analytically. :-)
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
     * GL: made public, to give BayesianSkylineGibbsOperator access
     */
    public final void setupIntervals() {

        if (intervals == null) {
            int maxIntervalCount = tree.getNodeCount();

            intervals = new double[maxIntervalCount];
            lineageCounts = new int[maxIntervalCount];
            storedIntervals = new double[maxIntervalCount];
            storedLineageCounts = new int[maxIntervalCount];
        }

        XTreeIntervals ti = new XTreeIntervals(intervals, lineageCounts);
        getTreeIntervals(tree, getMRCAOfCoalescent(tree), getExcludedMRCAs(tree), ti);
        intervalCount = ti.nIntervals;

        intervalsKnown = true;
    }


    /**
     * Extract coalescent times and tip information into ArrayList times from tree.
     * Upon return times contain the time of each node in the subtree below top, and at the corrosponding index
     * of childs is the descendent count for that time.
     *
     * @param top          the node to start from
     * @param excludeBelow an optional array of nodes to exclude (corresponding subtrees) from density.
     * @param tree         given tree
     * @param times        array to fill with times
     * @param childs       array to fill with descendents count
     */
    private static void collectAllTimes(Tree tree, NodeRef top, NodeRef[] excludeBelow,
                                        ArrayList<ComparableDouble> times, ArrayList<Integer> childs) {

        times.add(new ComparableDouble(tree.getNodeHeight(top)));
        childs.add(tree.getChildCount(top));

        for (int i = 0; i < tree.getChildCount(top); i++) {
            NodeRef child = tree.getChild(top, i);
            if (excludeBelow == null) {
                collectAllTimes(tree, child, excludeBelow, times, childs);
            } else {
                // check if this subtree is included in the coalescent density
                boolean include = true;
                for (NodeRef anExcludeBelow : excludeBelow) {
                    if (anExcludeBelow.getNumber() == child.getNumber()) {
                        include = false;
                        break;
                    }
                }
                if (include)
                    collectAllTimes(tree, child, excludeBelow, times, childs);
            }
        }
    }

    private class XTreeIntervals {

        public XTreeIntervals(double[] intervals, int[] lineageCounts) {
            this.intervals = intervals;
            this.lineagesCount = lineageCounts;
        }

        int nIntervals;
        final int[] lineagesCount;
        final double[] intervals;

    }

    private static void getTreeIntervals(Tree tree, NodeRef root, NodeRef[] exclude, XTreeIntervals ti) {
        double MULTIFURCATION_LIMIT = 1e-9;

        ArrayList<ComparableDouble> times = new ArrayList<ComparableDouble>();
        ArrayList<Integer> childs = new ArrayList<Integer>();
        collectAllTimes(tree, root, exclude, times, childs);
        int[] indices = new int[times.size()];

        HeapSort.sort(times, indices);

        final double[] intervals = ti.intervals;
        final int[] lineageCounts = ti.lineagesCount;

        // start is the time of the first tip
        double start = times.get(indices[0]).doubleValue();
        int numLines = 0;
        int i = 0;
        int intervalCount = 0;
        while (i < times.size()) {

            int lineagesRemoved = 0;
            int lineagesAdded = 0;

            final double finish = times.get(indices[i]).doubleValue();
            double next = finish;

            while (Math.abs(next - finish) < MULTIFURCATION_LIMIT) {
                final int children = childs.get(indices[i]);
                if (children == 0) {
                    lineagesAdded += 1;
                } else {
                    lineagesRemoved += (children - 1);
                }
                i += 1;
                if (i == times.size()) break;

                next = times.get(indices[i]).doubleValue();
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

        ti.nIntervals = intervalCount;
    }

    /**
     * @return number of intervals
     */
    public final int getIntervalCount() {
        return intervalCount;
    }

    /**
     * Gets an interval.
     *
     * @param i index of interval
     * @return interval length
     */
    public final double getInterval(int i) {
        if (i >= intervalCount) throw new IllegalArgumentException();
        return intervals[i];
    }

    /**
     * Returns the number of uncoalesced lineages within this interval.
     * Required for s-coalescents, where new lineages are added as
     * earlier samples are come across.
     *
     * @param i lineage index
     * @return number of uncoalesced lineages within this interval.
     */
    public final int getLineageCount(int i) {
        if (i >= intervalCount) throw new IllegalArgumentException();
        return lineageCounts[i];
    }

    /**
     * @param i interval index
     * @return the number coalescent events in an interval
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
     * @param i interval index
     * @return the type of interval observed.
     */
    public final CoalescentEventType getIntervalType(int i) {

        if (i >= intervalCount) throw new IllegalArgumentException();
        int numEvents = getCoalescentEvents(i);

        if (numEvents > 0) return CoalescentEventType.COALESCENT;
        else if (numEvents < 0) return CoalescentEventType.NEW_SAMPLE;
        else return CoalescentEventType.NOTHING;
    }

    /**
     * @return total height of the genealogy represented by these
     *         intervals.
     */
    public final double getTotalHeight() {

        double height = 0.0;
        for (int j = 0; j < intervalCount; j++) {
            height += intervals[j];
        }
        return height;
    }

    /**
     * @return whether this set of coalescent intervals is fully resolved
     *         (i.e. whether is has exactly one coalescent event in each
     *         subsequent interval)
     */
    public final boolean isBinaryCoalescent() {
        for (int i = 0; i < intervalCount; i++) {
            if (getCoalescentEvents(i) != 1) return false;
        }

        return true;
    }

    /**
     * @return whether this set of coalescent intervals coalescent only
     *         (i.e. whether is has exactly one or more coalescent event in each
     *         subsequent interval)
     */
    public final boolean isCoalescentOnly() {
        for (int i = 0; i < intervalCount; i++) {
            if (getCoalescentEvents(i) < 1) return false;
        }
        return true;
    }

    public String toString() {
        return getId(); // Double.toString(getLogLikelihood());
    }

    // **************************************************************
    // Units IMPLEMENTATION
    // **************************************************************

    /**
     * Sets the units these coalescent intervals are
     * measured in.
     */
    public final void setUnits(Type u) {
        demoModel.setUnits(u);
    }

    /**
     * Returns the units these coalescent intervals are
     * measured in.
     */
    public final Type getUnits() {
        return demoModel.getUnits();
    }

    public final boolean getIntervalsKnown() {
        return intervalsKnown;
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

//    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
//
//        public String getParserName() {
//            return COALESCENT_LIKELIHOOD;
//        }
//
//		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
//
//            XMLObject cxo = (XMLObject) xo.getChild(MODEL);
//            DemographicModel demoModel = (DemographicModel) cxo.getChild(DemographicModel.class);
//
//            List<TreeModel> trees = new ArrayList<TreeModel>();
//            List<Double> popFactors = new ArrayList<Double>();
//            MultiLociTreeSet treesSet = demoModel instanceof MultiLociTreeSet ? (MultiLociTreeSet)demoModel : null;
//
//            for(int k = 0; k < xo.getChildCount(); ++k) {
//                final Object child = xo.getChild(k);
//                if( child instanceof XMLObject ) {
//                    cxo = (XMLObject)child;
//                    if( cxo.getName().equals(POPULATION_TREE) ) {
//                        final TreeModel treeModel = (TreeModel) cxo.getChild(TreeModel.class);
//                        if( treeModel == null ) {
//                            // xml check not done yet?
//                            throw new XMLParseException("Expecting a tree model.");
//                        }
//                        trees.add(treeModel);
//
//                        try {
//                            double v = cxo.hasAttribute(POPULATION_FACTOR) ?
//                                    cxo.getDoubleAttribute(POPULATION_FACTOR) : 1.0;
//                            popFactors.add(v);
//                        } catch (XMLParseException e) {
//                            throw new XMLParseException(e.getMessage());
//                        }
//                    }
//                } else if( child instanceof MultiLociTreeSet )  {
//                    treesSet = (MultiLociTreeSet)child;
//                }
//            }
//
//            TreeModel treeModel = null;
//            if( trees.size() == 1 && popFactors.get(0) == 1.0 ) {
//                treeModel = trees.get(0);
//            } else if( trees.size() > 1 ) {
//               treesSet = new MultiLociTreeSet.Default(trees, popFactors);
//            } else if( !(trees.size() == 0 && treesSet != null) ) {
//               throw new XMLParseException("error");
//            }
//
//            if( treeModel != null ) {
//                return new OldAbstractCoalescentLikelihood(treeModel, demoModel);
//            }
//            return new OldAbstractCoalescentLikelihood(treesSet, demoModel);
//        }
//
//
//        //************************************************************************
//		// AbstractXMLObjectParser implementation
//		//************************************************************************
//
//        public String getParserDescription() {
//            return "This element represents the likelihood of the tree given the demographic function.";
//        }
//
//        public Class getReturnType() {
//            return Likelihood.class;
//        }
//
//        public XMLSyntaxRule[] getSyntaxRules() {
//            return rules;
//        }
//
//		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
//			new ElementRule(MODEL, new XMLSyntaxRule[] {
//				new ElementRule(DemographicModel.class)
//			}),
//			new ElementRule(POPULATION_TREE, new XMLSyntaxRule[] {
//                    AttributeRule.newDoubleRule(POPULATION_FACTOR, true),
//                new ElementRule(TreeModel.class)
//			}, 0, Integer.MAX_VALUE),
//		};
//	};

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
			return new CoalescentLikelihood(treeModel, demoModel);
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
    private DemographicModel demoModel = null;

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
    private boolean storedIntervalsKnown = false;

    double logLikelihood;
    private double storedLogLikelihood;
    boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;

    int intervalCount = 0;
    private int storedIntervalCount = 0;
}
