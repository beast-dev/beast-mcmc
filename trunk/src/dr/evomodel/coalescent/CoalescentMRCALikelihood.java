/*
 * CoalescentMRCALikelihood.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.math.Binomial;
import dr.util.ComparableDouble;
import dr.util.HeapSort;
import dr.xml.*;

import java.util.ArrayList;
import java.util.Set;


/**
 * A likelihood function for the coalescent. Takes a tree and a demographic model.
 *
 * Parts of this class were derived from C++ code provided by Oliver Pybus.
 *
 * @version $Id: CoalescentMRCALikelihood.java,v 1.13 2005/05/24 20:25:57 rambaut Exp $
 *
 * @author Alexei Drummond
 */
public class CoalescentMRCALikelihood extends AbstractModel implements Likelihood, Units {

	public static final String COALESCENT_MRCA_LIKELIHOOD = "coalescentMRCALikelihood";
	public static final String MODEL = "model";
	public static final String POPULATION_TREE = "populationTree";
	public static final String MRCA = "mrca";
	public static final String EXCLUDE = "exclude";

	/** Denotes an interval after which a coalescent event is observed
	  * (i.e. the number of lineages is smaller in the next interval) */
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


	public CoalescentMRCALikelihood(
		TreeModel treeModel,
		TaxonList ingroup,
		TaxonList[] excludeSubtrees,
		DemographicModel demoModel) throws Tree.MissingTaxonException
	{

		this(COALESCENT_MRCA_LIKELIHOOD, treeModel, ingroup, excludeSubtrees, demoModel);
	}

	public CoalescentMRCALikelihood(
		String name,
		TreeModel treeModel,
		TaxonList ingroup,
		TaxonList[] excludeSubtrees,
		DemographicModel demoModel) throws Tree.MissingTaxonException
	{
		super(name);

		this.treeModel = treeModel;
		this.demoModel = demoModel;

		addModel(treeModel);
		addModel(demoModel);

		ingroupLeafSet = Tree.Utils.getLeavesForTaxa(treeModel, ingroup);
		if (excludeSubtrees != null) {
			excludeLeafSets = new Set[excludeSubtrees.length];
			for (int i =0; i < excludeSubtrees.length; i++) {
				excludeLeafSets[i] = Tree.Utils.getLeavesForTaxa(treeModel, excludeSubtrees[i]);
			}
		}

		setupIntervals();
	}

	/**
	 * @returns the node ref of the MRCA of this coalescent prior in the given tree.
	 */
	public NodeRef getMRCAOfCoalescent(Tree tree) {
		return Tree.Utils.getCommonAncestorNode(tree, ingroupLeafSet);
	}

	/**
	 * @returns an array of noderefs that represent the MRCAs of subtrees to exclude from coalescent prior.
	 * May return null if no subtrees should be excluded.
	 */
	public NodeRef[] getExcludedMRCAs(Tree tree) {

		if (excludeLeafSets == null) return null;
		NodeRef[] excludeBelow = new NodeRef[excludeLeafSets.length];
		for (int i =0; i < excludeLeafSets.length; i++) {
			excludeBelow[i] = Tree.Utils.getCommonAncestorNode(tree, excludeLeafSets[i]);
		}
		return excludeBelow;
	}

	// **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

	protected final void handleModelChangedEvent(Model model, Object object, int index) {
		if (model == treeModel) {
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

	protected final void handleParameterChangedEvent(Parameter parameter, int index) { } // No parameters to respond to

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

	protected final void acceptState() { } // nothing to do

	// **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

	public final Model getModel() { return this; }

	public final double getLogLikelihood() {
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

		if (demoModel == null) return calculateAnalyticalLogLikelihood();

		double logL = 0.0;

		double currentTime = 0.0;

		DemographicFunction demoFunction = demoModel.getDemographicFunction();

		for (int j = 0; j < intervalCount; j++) {

			logL += calculateIntervalLikelihood(demoFunction, intervals[j], currentTime, lineageCounts[j],
			getIntervalType(j));

			// insert zero-length coalescent intervals
			int diff = getCoalescentEvents(j)-1;
			for (int k = 0; k < diff; k++) {
				logL += calculateIntervalLikelihood(demoFunction, 0.0, currentTime, lineageCounts[j]-k-1, COALESCENT);
			}

			currentTime += intervals[j];
		}

		return logL;
	}

	private final double calculateAnalyticalLogLikelihood() {

		double lambda = getLambda();
		int n = treeModel.getExternalNodeCount();

		double logL = 0.0;

		// assumes a 1/theta prior
		//logLikelihood = Math.log(1.0/Math.pow(lambda,n));

		// assumes a flat prior
		logL = Math.log(1.0/Math.pow(lambda,n-1));
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
												int lineageCount, int type)
	{
		//binom.setMax(lineageCount);

		double timeOfThisCoal = width + timeOfPrevCoal;

		double intervalArea = demoFunction.getIntegral(timeOfPrevCoal, timeOfThisCoal);
		double like = 0;
		switch (type) {
			case COALESCENT:
				like = - Math.log(demoFunction.getDemographic(timeOfThisCoal)) -
								(Binomial.choose2(lineageCount)*intervalArea);
				break;
			case NEW_SAMPLE:
				like = -(Binomial.choose2(lineageCount)*intervalArea);
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
		for (int i= 0; i < getIntervalCount(); i++) {
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
		collectAllTimes(treeModel, getMRCAOfCoalescent(treeModel), getExcludedMRCAs(treeModel), times, childs);
		int[] indices = new int[times.size()];

		HeapSort.sort(times, indices);

		int maxIntervalCount = treeModel.getNodeCount();

		if (intervals == null) {
			intervals = new double[maxIntervalCount];
			lineageCounts = new int[maxIntervalCount];
			storedIntervals = new double[maxIntervalCount];
			storedLineageCounts = new int[maxIntervalCount];
		}

		// start is the time of the first tip
		double start = ((ComparableDouble)times.get(indices[0])).doubleValue();
		int numLines = 0;
		int i = 0;
		intervalCount = 0;
		while (i < times.size()) {

			int lineagesRemoved = 0;
			int lineagesAdded = 0;

			double finish = ((ComparableDouble)times.get(indices[i])).doubleValue();
			double next = finish;

			while (Math.abs(next - finish) < MULTIFURCATION_LIMIT) {
				int children = ((Integer)childs.get(indices[i])).intValue();
				if (children == 0) {
					lineagesAdded += 1;
				} else {
					lineagesRemoved += (children - 1);
				}
				i += 1;
				if (i < times.size()) {
					next = ((ComparableDouble)times.get(indices[i])).doubleValue();
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
	 * @param node the node to start from
	 * @param excludeBelow an optional array of nodes to exclude (corresponding subtrees) from density.
	 */
	private final static void collectAllTimes(Tree tree, NodeRef node, NodeRef[] excludeBelow, ArrayList times, ArrayList childs) {

		times.add(new ComparableDouble(tree.getNodeHeight(node)));
		childs.add(new Integer(tree.getChildCount(node)));

		for (int i = 0; i < tree.getChildCount(node); i++) {
			NodeRef child = tree.getChild(node, i);
			if (excludeBelow == null) {
				collectAllTimes(tree, child, excludeBelow, times, childs);
			} else {
				// check if this subtree is included in the coalescent density
				boolean include = true;
				for (int j =0; j < excludeBelow.length; j++) {
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
		if (i < intervalCount-1) {
			return lineageCounts[i]-lineageCounts[i+1];
		} else {
			return lineageCounts[i]-1;
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

		double height=0.0;
		for (int j=0; j < intervalCount; j++) {
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
		return new dr.inference.loggers.LogColumn[] {
			new LikelihoodColumn(getId())
		};
	}

	private final class LikelihoodColumn extends dr.inference.loggers.NumberColumn {
		public LikelihoodColumn(String label) { super(label); }
		public double getDoubleValue() { return getLogLikelihood(); }
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
	public final void setUnits(int u)
	{
		demoModel.setUnits(u);
	}

	/**
	 * Returns the units these coalescent intervals are
	 * measured in.
	 */
	public final int getUnits()
	{
		return demoModel.getUnits();
	}

	// ****************************************************************
	// Private static methods
	// ****************************************************************

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return COALESCENT_MRCA_LIKELIHOOD; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			DemographicModel demoModel = (DemographicModel)xo.getSocketChild("model");
			TreeModel treeModel = (TreeModel)xo.getSocketChild("populationTree");
			TaxonList ingroup = (TaxonList)xo.getSocketChild("mrca");

			TaxonList[] excludeSubtrees = null;
			ArrayList excludeTaxa = new ArrayList();


			// should have one child that is node
			for (int i =0; i < xo.getChildCount(); i++) {
				Object child = xo.getChild(i);
				if ((child instanceof XMLObject) && ((XMLObject)child).getName().equals(EXCLUDE)) {

					XMLObject xchild = (XMLObject)child;

					for (int j = 0; j < xchild.getChildCount(); j++) {
						if  (xchild.getChild(j) instanceof TaxonList) {
							excludeTaxa.add(xchild.getChild(j));
						}
					}
					if (excludeTaxa.size() == 0) {
						throw new XMLParseException(EXCLUDE + " element should contain one or more taxa element");
					} else {
						excludeSubtrees = new TaxonList[excludeTaxa.size()];
						for (int j = 0; j < excludeSubtrees.length; j++) {
							excludeSubtrees[j] = (TaxonList)excludeTaxa.get(j);
						}
					}
				}
			}

			try {
				return new CoalescentMRCALikelihood(treeModel, ingroup, excludeSubtrees, demoModel);
			} catch (Tree.MissingTaxonException mte) {
				throw new XMLParseException("treeModel missing a taxon from taxon list in " + getParserName() + " element");
			}
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "A coalescent likelihood function for a subtree.";
		}

		public Class getReturnType() { return CoalescentMRCALikelihood.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new ElementRule(MODEL, DemographicModel.class),
			new ElementRule(POPULATION_TREE, TreeModel.class),
			new ElementRule(MRCA,
				new XMLSyntaxRule[] { new ElementRule(Taxa.class) }),
			new ElementRule(EXCLUDE,
				new XMLSyntaxRule[] { new ElementRule(Taxa.class, 1, Integer.MAX_VALUE) }, true)
		};
	};

	/** The demographic model. */
	DemographicModel demoModel = null;

	Set ingroupLeafSet = null;
	Set[] excludeLeafSets = null;

	/** The treeModel. */
	TreeModel treeModel = null;

	/** The widths of the intervals. */
	double[] intervals;
	private double[] storedIntervals;

	/** The number of uncoalesced lineages within a particular interval. */
	int[] lineageCounts;
	private int[] storedLineageCounts;

	private boolean intervalsKnown = false;
	private boolean storedIntervalsKnown = false;

	private double logLikelihood;
	private double storedLogLikelihood;
	private boolean likelihoodKnown = false;
	private boolean storedLikelihoodKnown = false;

	private int intervalCount = 0;
	private int storedIntervalCount = 0;

	private Binomial binom = new Binomial();
}