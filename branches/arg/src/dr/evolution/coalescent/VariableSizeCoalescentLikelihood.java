/*
 * CoalescentLikelihood.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */

package dr.evolution.coalescent;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.DemographicModel;
import dr.evomodel.tree.ARGModel;
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
 * @version $Id: CoalescentWithRecombinationLikelihood.java,v 1.1.1.1.2.2 2006/11/06 01:38:30 msuchard Exp $
 */
public class VariableSizeCoalescentLikelihood extends AbstractModel implements Likelihood, Units {

    // PUBLIC STUFF

    public static final String COALESCENT_LIKELIHOOD = "variableSizeCoalescentLikelihood";
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


    public VariableSizeCoalescentLikelihood(Tree tree, DemographicModel demoModel) {
        this(COALESCENT_LIKELIHOOD, tree, demoModel, true);
    }

    public VariableSizeCoalescentLikelihood(String name, Tree tree, DemographicModel demoModel, boolean setupIntervals) {

        super(name);

        this.tree = tree;
        this.demoModel = demoModel;
        if (tree instanceof ARGModel) {
            addModel((ARGModel) tree);
        }
        if (demoModel != null) {
            addModel(demoModel);
        }
        if (setupIntervals) setupIntervals();

        addStatistic(new DeltaStatistic());
    }

    VariableSizeCoalescentLikelihood(String name) {
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

    protected final void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == tree) {
            // treeModel has changed so recalculate the intervals
            intervalsKnown = false;
        } else {
            // demoModel has changed so we don't need to recalculate the intervals
        }

        likelihoodKnown = false;
    }

    // **************************************************************
    // ParameterListener IMPLEMENTATION
    // **************************************************************

    protected final void handleParameterChangedEvent(Parameter parameter, int index) {
    } // No parameters to respond to

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    /**
     * Stores the precalculated state: in this case the intervals
     */
    protected final void storeState() {
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
    protected final void restoreState() {
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
        System.err.print("COALESCENT PRIOR: Start");

        if (true) return calculateAnalyticalLogLikelihood();

        //System.err.println(tree.getExternalNodeCount()+" tips");
        if (intervalsKnown == false) setupIntervals();

        //if (demoModel == null) return calculateAnalyticalLogLikelihood();
        if (true) return calculateAnalyticalLogLikelihood();

        double logL = 0.0;
        //double logL = calculateAnalyticalLogLikelihood();

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

    private final double calculateAnalyticalLogLikelihood() {
        //System.err.println("What the hell?");
        //System.exit(-1);
        double lambda = getLambda();
        int n = ((ARGModel) tree).getReassortmentNodeCount();
        System.err.println("Prior for " + n + " reassortments.");
        double treeHeight = tree.getNodeHeight(tree.getRoot());
        double logL = -lambda * treeHeight;
        double p = 0.00000000000000005;
        logL += n * Math.log(p);

        // assumes a 1/theta prior
        //logLikelihood = Math.log(1.0/Math.pow(lambda,n));

        // assumes a flat prior
        //logL = Math.log(1.0/Math.pow(lambda,n-1));
        System.err.println("COALESCENT PRIOR: end");
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

        double intervalArea = demoFunction.getIntegral(timeOfPrevCoal, timeOfThisCoal);
        double like = 0;
        switch (type) {
            case COALESCENT:
                like = -Math.log(demoFunction.getDemographic(timeOfThisCoal)) -
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
    private final double getLambda() {
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
        ArrayList childs = new ArrayList();
        //System.err.println("Collecting times");
        collectAllTimes(tree, getMRCAOfCoalescent(tree), getExcludedMRCAs(tree), times, childs);
        int[] indices = new int[times.size()];

        //System.err.println("Sorting times");
        HeapSort.sort(times, indices);

        int maxIntervalCount = tree.getNodeCount(); //

//		System.err.println("max l = "+maxIntervalCount);

        if ((intervals == null) || (intervals.length != maxIntervalCount)) {
            intervals = new double[maxIntervalCount];
            lineageCounts = new int[maxIntervalCount];
            storedIntervals = new double[maxIntervalCount];
            storedLineageCounts = new int[maxIntervalCount];
        }

//		System.err.println("Start Processing");
//		System.err.println("yo = "+indices[0]);
        // start is the time of the first tip
        int indicesLength = indices.length - 1;
        double start = ((ComparableDouble) times.get(indices[indicesLength])).doubleValue();
//		System.err.println("start = "+start);
        int numLines = 1; // was 0
//		int i = 0;
        intervalCount = 0;
        while (start > 0) {

            double finish = ((ComparableDouble) times.get(indices[indicesLength - 1])).doubleValue();
//			System.err.println("finish = "+finish);
//			double next = finish;
            int children = ((Integer) childs.get(indices[indicesLength])).intValue();
            if (children == 2)
                numLines++;
            if (children == 1)
                numLines--;
            if (children == 0)
                intervals[intervalCount] = start - finish;
            lineageCounts[intervalCount] = numLines;
            start = finish;
            intervalCount++;
            indicesLength--;
        }
//		System.err.println("intervalCount = "+intervalCount);	
//			int lineagesRemoved = 0;
//			int lineagesAdded = 0;
//
//			double finish = ((ComparableDouble)times.get(indices[i])).doubleValue();
//			double next = finish;
//
//			while (Math.abs(next - finish) < MULTIFURCATION_LIMIT) {
//				int children = ((Integer)childs.get(indices[i])).intValue();
//				if (children == 0) {
//					lineagesAdded += 1;
//				} else {
//					lineagesRemoved += (children - 1);
//				}
//				i += 1;
//				if (i < times.size()) {
//					next = ((ComparableDouble)times.get(indices[i])).doubleValue();
//				} else break;
//			}
//			//System.out.println("time = " + finish + " removed = " + lineagesRemoved + " added = " + lineagesAdded);
//			if (lineagesAdded > 0) {
//
//				if (intervalCount > 0 || ((finish - start) > MULTIFURCATION_LIMIT)) {
//					intervals[intervalCount] = finish - start;
//					lineageCounts[intervalCount] = numLines;
//					intervalCount += 1;
//				}
//
//				start = finish;
//			}
//			// add sample event
//			numLines += lineagesAdded;
//
//			if (lineagesRemoved > 0) {
//
//				intervals[intervalCount] = finish - start;
//				lineageCounts[intervalCount] = numLines;
//				intervalCount += 1;
//				start = finish;
//			}
//			// coalescent event
//			numLines -= lineagesRemoved;
        //}
        intervalsKnown = true;
//		printIntervals();
        //System.exit(-1);
    }

    private void printIntervals() {
        System.err.println("DEBUG: Coalescent intervals");
        int len = intervalCount;
        for (int i = 0; i < len; i++) {
            System.err.println(intervals[i] + " : " + lineageCounts[i]);
        }
    }

    /**
     * extract coalescent times and tip information into ArrayList times from tree.
     *
     * @param node         the node to start from
     * @param excludeBelow an optional array of nodes to exclude (corresponding subtrees) from density.
     */
    private final static void collectAllTimes(Tree tree, NodeRef node, NodeRef[] excludeBelow, ArrayList times, ArrayList childs) {

        for (int i = 0, n = tree.getNodeCount(); i < n; i++) {
            NodeRef nr = tree.getNode(i);
            times.add(new ComparableDouble(tree.getNodeHeight(nr)));
            childs.add(new Integer(tree.getChildCount(nr)));

        }
    }

//		times.add(new ComparableDouble(tree.getNodeHeight(node)));
//		childs.add(new Integer(tree.getChildCount(node)));
//
//		for (int i = 0; i < tree.getChildCount(node); i++) {
//			NodeRef child = tree.getChild(node, i);
//			if (excludeBelow == null) {
//				collectAllTimes(tree, child, excludeBelow, times, childs);
//			} else {
//				// check if this subtree is included in the coalescent density
//				boolean include = true;
//				for (int j =0; j < excludeBelow.length; j++) {
//					if (excludeBelow[j].getNumber() == child.getNumber()) {
//						include = false;
//						break;
//					}
//				}
//				if (include) collectAllTimes(tree, child, excludeBelow, times, childs);
//			}
//		}
    //}


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

    // **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

    /**
     * @return the log columns.
     */
    public final dr.inference.loggers.LogColumn[] getColumns() {
        return new dr.inference.loggers.LogColumn[]{
                new LikelihoodColumn(getId())
        };
    }

    private final class LikelihoodColumn extends dr.inference.loggers.NumberColumn {
        public LikelihoodColumn(String label) {
            super(label);
        }

        public double getDoubleValue() {
            return getLogLikelihood();
        }
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public org.w3c.dom.Element createElement(org.w3c.dom.Document d) {
        throw new RuntimeException("createElement not implemented");
    }

    // **************************************************************
    // Units IMPLEMENTATION
    // **************************************************************

    /**
     * Sets the units these coalescent intervals are
     * measured in.
     */
    public final void setUnits(int u) {
        demoModel.setUnits(u);
    }

    /**
     * Returns the units these coalescent intervals are
     * measured in.
     */
    public final int getUnits() {
        return demoModel.getUnits();
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

            DemographicModel demoModel = null;
            XMLObject cxo;
            try {
                cxo = (XMLObject) xo.getChild(MODEL);
                demoModel = (DemographicModel) cxo.getChild(DemographicModel.class);
            } catch (Exception e) {
            }  // [TODO -- how to do this right?
            cxo = (XMLObject) xo.getChild(POPULATION_TREE);
            ARGModel argModel = (ARGModel) cxo.getChild(ARGModel.class);

            return new VariableSizeCoalescentLikelihood(argModel, demoModel);
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

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                //	new ElementRule(MODEL, new XMLSyntaxRule[] {
                //		new ElementRule(DemographicModel.class)
                //	}),
                new ElementRule(POPULATION_TREE, new XMLSyntaxRule[]{
                        new ElementRule(ARGModel.class)
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
    private boolean storedIntervalsKnown = false;

    double logLikelihood;
    private double storedLogLikelihood;
    boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;

    int intervalCount = 0;
    private int storedIntervalCount = 0;
}