/*
 * CoalescentLikelihood.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */

package dr.evomodel.coalescent;

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.tree.ARGModel;
import dr.inference.model.*;
import dr.math.Binomial;
import dr.util.ComparableDouble;
import dr.util.HeapSort;
import dr.xml.*;

import java.util.ArrayList;


/**
 * A likelihood function for the coalescent with recombination. Takes a tree and a demographic model.
 *
 * @author Marc Suchard
 * @version $Id: CoalescentWithRecombinationLikelihood.java,v 1.1.1.1.2.2 2006/11/06 01:38:30 msuchard Exp $
 */
public class CoalescentWithRecombinationLikelihood extends AbstractModel implements Likelihood, Units {

    // PUBLIC STUFF

    public static final String COALESCENT_LIKELIHOOD = "variableSizeCoalescentLikelihood";
    public static final String ANALYTICAL = "analytical";
    public static final String COALESCENT_MODEL = "coalescentModel";
    public static final String RECOMBINATION_MODEL = "recombinationModel";
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
    public static final int RECOMBINATION = 1;

    /**
     * Denotes an interval at the end of which nothing is
     * observed (i.e. the number of lineages is the same in the next interval).
     */
    public static final int NOTHING = 2;


    public CoalescentWithRecombinationLikelihood(Tree tree,
                                                 DemographicModel coalescentDemoModel,
                                                 DemographicModel recombinationDemoModel) {
        this(COALESCENT_LIKELIHOOD, tree, coalescentDemoModel, recombinationDemoModel, true);
    }

    public CoalescentWithRecombinationLikelihood(String name, Tree tree,
                                                 DemographicModel coalescentDemoModel,
                                                 DemographicModel recombinationDemoModel,
                                                 boolean setupIntervals) {

        super(name);

        this.tree = tree;
        this.coalescentDemoModel = coalescentDemoModel;
        this.recombinationDemoModel = recombinationDemoModel;
        if (tree instanceof ARGModel) {
            addModel((ARGModel) tree);
        }
        if (coalescentDemoModel != null) {
            addModel(coalescentDemoModel);
        }
        if (setupIntervals) setupIntervals();

        addStatistic(new DeltaStatistic());
    }

    CoalescentWithRecombinationLikelihood(String name) {
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

        if (intervalsKnown == false) setupIntervals();


        if (coalescentDemoModel == null || recombinationDemoModel == null) return calculateAnalyticalLogLikelihood();

        double logL = 0.0;

        double currentTime = 0.0;

        DemographicFunction coalescentDemoFunction = coalescentDemoModel.getDemographicFunction();
        DemographicFunction recombinationDemoFunction = recombinationDemoModel.getDemographicFunction();

        for (int j = 0; j < intervalCount; j++) {

            logL += calculateIntervalLikelihood(coalescentDemoFunction, recombinationDemoFunction, intervals[j],
                    currentTime, lineageCounts[j],
                    getIntervalType(j));

            // insert zero-length coalescent intervals
            // assume no zero-length coalescent intervals
            /*           int diff = getCoalescentEvents(j) - 1;
            for (int k = 0; k < diff; k++) {
                logL += calculateIntervalLikelihood(demoFunction, 0.0, currentTime, lineageCounts[j] - k - 1, COALESCENT);
            }*/

            currentTime += intervals[j];
        }

//        printIntervals();
//        ARGModel arg = (ARGModel) tree;
//        if (arg.getReassortmentNodeCount() > 0)   {
//            System.err.println(arg.toGraphString());
//            System.err.println("@@@ = "+intervalCount);
//            System.exit(-1);
//        }


        return logL;
    }

    private final double calculateAnalyticalLogLikelihood() {
        //System.err.println("What the hell?");
        //System.exit(-1);
        double lambda = getLambda();
        int n = ((ARGModel) tree).getReassortmentNodeCount();
//		System.err.println("Prior for " + n + " reassortments.");
        double p = 0.001;
        double logL = n * Math.log(p);

        double height = tree.getNodeHeight(tree.getRoot());
        logL -= height * 1.0;   // Exp(1)
        // assumes a flat prior  over branch lengths given treeHeight
        //logL = Math.log(1.0/Math.pow(lambda,n-1));
//		System.err.println("COALESCENT PRIOR: end yoyoy");
        return logL;
    }

    /**
     * Returns the likelihood of a given *coalescent* interval
     */
    public final double calculateIntervalLikelihood(DemographicFunction coalescentDemoFunction,
                                                    DemographicFunction recombinationDemoFunction,
                                                    double width, double timeOfPrevCoal, int lineageCount) {

        return calculateIntervalLikelihood(coalescentDemoFunction, recombinationDemoFunction,
                width, timeOfPrevCoal, lineageCount, COALESCENT);
    }

    /**
     * Returns the likelihood of a given interval,coalescent or otherwise.
     */
    public final double calculateIntervalLikelihood(DemographicFunction coalescentDemoFunction,
                                                    DemographicFunction recombinationDemoFunction,
                                                    double width, double timeOfPrevEvent,
                                                    int lineageCount, int type) {
        //binom.setMax(lineageCount);

        double recombinationWeight = 1;    // todo get value from recombinationDemoModel

        double timeOfThisEvent = width + timeOfPrevEvent;
//         System.err.printf("s: %7.6f   f: %7.6f,  %d, %d\n",timeOfPrevCoal,timeOfThisCoal,lineageCount, type);

        double intervalAreaCoalescent = coalescentDemoFunction.getIntegral(timeOfPrevEvent, timeOfThisEvent);
        double intervalAreaRecombination = recombinationDemoFunction.getIntegral(timeOfPrevEvent, timeOfThisEvent);
        double like = 0;
        switch (type) {
            case COALESCENT:
                like =  // coalescence occurs at timeOfThisEvent
                        -Math.log(coalescentDemoFunction.getDemographic(timeOfThisEvent))
                                - (Binomial.choose2(lineageCount) * intervalAreaCoalescent)
                                // and recombination did not occur in timeInterval (width)
                                - (lineageCount / 2.0 * intervalAreaRecombination) * recombinationWeight; // intervalAreaRecombination = \effectiveRho * timeInterval
                break;
            case RECOMBINATION:
                like =  // recombination occurs at timeOfThisEvent
                        -Math.log(recombinationDemoFunction.getDemographic(timeOfThisEvent) * recombinationWeight)
                                - (lineageCount / 2.0 * intervalAreaRecombination) * recombinationWeight // todo do I need to rescale?
                                // and coalescence does not occur in timeInterval (width)
                                - (Binomial.choose2(lineageCount) * intervalAreaCoalescent);
                break;
                // todo probably need to add in a relative weight/rate function
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

        ArrayList times = new ArrayList();
        ArrayList childs = new ArrayList();
        collectAllTimes(tree, getMRCAOfCoalescent(tree), getExcludedMRCAs(tree), times, childs);
        int[] indices = new int[times.size()];

        HeapSort.sort(times, indices);

//        for (int i : indices)
//            System.err.println(times.get(i).toString());

        int maxIntervalCount = tree.getNodeCount();

        if (intervals == null || intervals.length != maxIntervalCount) {
            intervals = new double[maxIntervalCount];
            lineageCounts = new int[maxIntervalCount];
            storedIntervals = new double[maxIntervalCount];
            storedLineageCounts = new int[maxIntervalCount];
        }

        // The following only works for strictly bifurcating trees nested in the ARG


        int index = 0;
        double start = ((ComparableDouble) times.get(indices[index])).doubleValue(); // start at most recent time
        double end = tree.getNodeHeight(getMRCAOfCoalescent(tree));                 // end at tMRCA
        int numLines = tree.getExternalNodeCount();
        intervalCount = 0;
        while (start < end) {
            double finish = ((ComparableDouble) times.get(indices[index + 1])).doubleValue();
            if (finish > start) { // Avoid repeated intervals due to reassortment nodes
                int children = ((Integer) childs.get(indices[index]));
                if (children == 2)
                    numLines--;
                else if (children == 1)
                    numLines++;
                /* else if (children == 0)
            ; // what about here*/  // todo
                intervals[intervalCount] = finish - start;
                lineageCounts[intervalCount] = numLines;
                intervalCount++;
            }
            index++;
            start = finish;
        }
        intervalsKnown = true;
    }

    private void printIntervals() {
        System.err.println("DEBUG: Coalescent intervals");
        int len = intervalCount;
        for (int i = 0; i < len; i++) {
            System.err.println(intervals[i] + " : " + lineageCounts[i]);
        }
        System.err.println("END DEBUG");
    }

    /**
     * extract coalescent times and tip information into ArrayList times from tree.
     *
     * @param node         the node to start from
     * @param excludeBelow an optional array of nodes to exclude (corresponding subtrees) from density.
     */
    private final static void collectAllTimes(Tree tree, NodeRef node, NodeRef[] excludeBelow, ArrayList times, ArrayList childs) {

        /*for (int i = 0, n = tree.getNodeCount(); i < n; i++) {
            NodeRef nr = tree.getNode(i);
            times.add(new ComparableDouble(tree.getNodeHeight(nr)));
            childs.add(new Integer(tree.getChildCount(nr)));

        }*/

        //       ARGModel.Node tmp = (ARGModel.Node) node;
        //       if (tmp.taxon != null) {
        //           System.err.println(tmp.taxon.getId()+" "+tree.getChildCount(node));
        //       }


        times.add(new ComparableDouble(tree.getNodeHeight(node)));
        childs.add(new Integer(tree.getChildCount(node)));
//        System.err.println(tree.getChildCount(node)) ;

        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef child = tree.getChild(node, i);
            if (excludeBelow == null) {
                collectAllTimes(tree, child, excludeBelow, times, childs);
            } else {
                // check if this subtree is included in the coalescent density
                boolean include = true;
                for (int j = 0; j < excludeBelow.length; j++) {
                    if (excludeBelow[j].getNumber() == child.getNumber()) {
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
        else if (numEvents < 0) return RECOMBINATION;
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
        coalescentDemoModel.setUnits(u);
    }

    /**
     * Returns the units these coalescent intervals are
     * measured in.
     */
    public final int getUnits() {
        return coalescentDemoModel.getUnits();
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

            DemographicModel coalescentDemoModel = null;
            XMLObject cxo;
            try {
                cxo = (XMLObject) xo.getChild(COALESCENT_MODEL);
                coalescentDemoModel = (DemographicModel) cxo.getChild(DemographicModel.class);
            } catch (Exception e) {
            }  // [TODO -- how to do this right?
            DemographicModel recombinationDemoModel = null;
            try {
                cxo = (XMLObject) xo.getChild(RECOMBINATION_MODEL);
                recombinationDemoModel = (DemographicModel) cxo.getChild(DemographicModel.class);
            } catch (Exception e) {
            }
            cxo = (XMLObject) xo.getChild(POPULATION_TREE);
            ARGModel argModel = (ARGModel) cxo.getChild(ARGModel.class);

            return new CoalescentWithRecombinationLikelihood(argModel, coalescentDemoModel, recombinationDemoModel);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents the likelihood of the ancestral recombination graph given the demographic function.";
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


    /**
     * The demographic models.
     */
    DemographicModel coalescentDemoModel = null;
    DemographicModel recombinationDemoModel = null;

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