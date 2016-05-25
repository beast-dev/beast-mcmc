/*
 * ARGCoalescentLikelihood.java
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

import dr.evolution.tree.NodeRef;
import dr.evomodel.arg.ARGModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.logging.Logger;

public class ARGCoalescentLikelihood extends VeryOldCoalescentLikelihood {

	public static final String ARG_COALESCENT_MODEL = "argCoalescentLikelihood";
	public static final String RECOMBINATION_RATE = "recombinationRate";
	public static final String POPULATION_SIZE = "populationSize";
	public static final String ARG_MODEL = "argModel";
	public static final String MAX_REASSORTMENTS = "maxReassortments";
	public static final String ANCESTRAL_RESTRICTION = "ancestralRestriction";

	public static final int RECOMBINATION = 3;

	private Parameter popSize;
	private Parameter recomRate;
	protected ARGModel arg;
	private int taxaNumber;
	protected int maxReassortments;
	private boolean ancestralRestriction = false;

	private ArrayList<CoalescentInterval> intervals;
	private ArrayList<CoalescentInterval> storedIntervals;

	public ARGCoalescentLikelihood(String name, ARGModel arg, int max) {
		super(name);
		this.arg = arg;

		intervals = new ArrayList<CoalescentInterval>();

		taxaNumber = arg.getExternalNodeCount();
		this.maxReassortments = max;
	}

	public ARGCoalescentLikelihood(Parameter popSize, Parameter recomRate,
	                               ARGModel arg, boolean setupIntervals, int maxReassort, boolean ancestralRestriction) {
		super(ARG_COALESCENT_MODEL);

		this.popSize = popSize;
		this.recomRate = recomRate;
		this.arg = arg;
		this.ancestralRestriction = ancestralRestriction;
		addVariable(popSize);
		addVariable(recomRate);

		addModel(arg);
		intervals = new ArrayList<CoalescentInterval>(arg.getNodeCount());
		intervalsKnown = false;
		likelihoodKnown = false;

		if (setupIntervals) {
			intervalsKnown = true;
			calculateIntervals();
		}

		taxaNumber = arg.getExternalNodeCount();
		this.maxReassortments = maxReassort;
	}

	public void calculateIntervals() {
		intervals.clear();
		intervals.ensureCapacity(arg.getNodeCount());

		NodeRef x;
		for (int i = 0; i < arg.getInternalNodeCount(); i++) {
			x = arg.getInternalNode(i);
			if (arg.isReassortment(x)) {
				intervals.add(new CoalescentInterval(arg.getNodeHeight(x), RECOMBINATION));
			} else {
				intervals.add(new CoalescentInterval(arg.getNodeHeight(x), COALESCENT));
			}
		}
		for (int i = 0; i < arg.getExternalNodeCount(); i++) {
			x = arg.getExternalNode(i);
			if (arg.getNodeHeight(x) > 0.0) {
				intervals.add(new CoalescentInterval(arg.getNodeHeight(x), NEW_SAMPLE));
			}
		}
		dr.util.HeapSort.sort(intervals);


		double a = 0, b = 0;
		for (int i = 0; i < intervals.size(); i++) {
			b = intervals.get(i).length;
			intervals.get(i).length = intervals.get(i).length - a;
			a = b;
		}

		intervalsKnown = true;
	}

	public void handleModelChangedEvent(Model model, Object object, int index) {
		if (model == arg) {
			intervalsKnown = false;
		}
		likelihoodKnown = false;
	}

	public void handleParameterChangedEvent(Parameter parameter, int index) {
		likelihoodKnown = false;
	}

	public void storeState() {
		storedIntervals = new ArrayList<CoalescentInterval>(intervals.size());
		for (CoalescentInterval interval : intervals) {
			storedIntervals.add(interval.clone());
		}

		intervalsKnown = likelihoodKnown = false;


		storedIntervalsKnown = intervalsKnown;
		storedLikelihoodKnown = likelihoodKnown;
		storedLogLikelihood = logLikelihood;

	}

	public void restoreState() {
		intervals = storedIntervals;
		storedIntervals.clear();
		intervalsKnown = storedIntervalsKnown;
		likelihoodKnown = storedLikelihoodKnown;
		logLikelihood = storedLogLikelihood;

		intervalsKnown = likelihoodKnown = false;
//	     if (!intervalsKnown) {
//	         likelihoodKnown = false;
//	     }
	}

	public boolean currentARGValid(boolean allowDoubleParents) {
		if (!intervalsKnown) {
			calculateIntervals();
		}
		int taxa = taxaNumber;

		for (CoalescentInterval x : intervals) {
			if (taxa == 1)
				return false;
			if (x.type == COALESCENT)
				taxa--;
			else if (x.type == RECOMBINATION)
				taxa++;
			else
				throw new RuntimeException("Not implemented yet");
		}
		if (!allowDoubleParents) {
			for (int i = 0, n = arg.getNodeCount(); i < n; i++) {
				NodeRef x = arg.getNode(i);
				if (arg.isReassortment(x) && arg.getParent(x, 0) == arg.getParent(x, 1)) {
					return false;
				}
			}
		}
		
		return true;
	}

	public double getLogLikelihood() {
		if (likelihoodKnown)
			return logLikelihood;
		if (!intervalsKnown)
			calculateIntervals();

		likelihoodKnown = true;

		if (arg.getReassortmentNodeCount() > maxReassortments){
			logLikelihood = Double.NEGATIVE_INFINITY;
		}else if(ancestralRestriction && !arg.isAncestral()){
			logLikelihood = Double.NEGATIVE_INFINITY;
		}else{
			logLikelihood = calculateLogLikelihood(
					popSize.getParameterValue(0),
					recomRate.getParameterValue(0));
		}


		return logLikelihood;
	}

	private double chooseTwo(int n) {
		return n * (n - 1) / 2.0;
	}

	private double calculateLogLikelihood(double pSize, double rRate) {
		
		double logLike = 0.0;
		int numberOfTaxa = taxaNumber;

		for (CoalescentInterval interval : intervals) {

			if (numberOfTaxa == 1)
				return Double.NEGATIVE_INFINITY;

			double rate = (double) numberOfTaxa *
					(numberOfTaxa - 1 + rRate) / (2.0 * pSize);

			logLike += Math.log(rate) - rate * interval.length;

			if (interval.type == COALESCENT) {
				logLike += Math.log((double) (numberOfTaxa - 1) /
						(numberOfTaxa - 1 + rRate))
						- Math.log(chooseTwo(numberOfTaxa));
				numberOfTaxa--;
			} else if (interval.type == RECOMBINATION) {
				logLike += Math.log(rRate / (numberOfTaxa - 1 + rRate))
						- Math.log((double) numberOfTaxa);

				numberOfTaxa++;
			} else {
				throw new RuntimeException("Not implemented yet");
			}
		}

		assert numberOfTaxa == 1;

		return logLike;
	}

	private class CoalescentInterval implements Comparable<CoalescentInterval>,
			Cloneable {
		public int type;
		public double length;

		public CoalescentInterval(double length, int type) {
			this.length = length;
			this.type = type;
		}

		public int compareTo(CoalescentInterval a) {

			if (a.length > this.length) {
				return -1;
			} else if (a.length == this.length) {
				Logger.getLogger("dr.evomodel.coalescent").severe(
						"The current ARG Model has 2 internal nodes " +
								"at the same height");
				return 0;
			}
			return 1;
		}

		public String toString() {
			if (type == 0) {
				return "(" + length + ", Coalescent)";
			}
			return "(" + length + ", Recombination)";
		}

		public CoalescentInterval clone() {
			return new CoalescentInterval(length, type);
		}

	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserDescription() {
			return "A coalescent likelihood for an ARG model";
		}

		public Class getReturnType() {
			return ARGCoalescentLikelihood.class;
		}

		public String getParserName() {
			return ARG_COALESCENT_MODEL;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private final XMLSyntaxRule[] rules = {
				new ElementRule(POPULATION_SIZE,
						new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
				new ElementRule(RECOMBINATION_RATE,
						new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
				new ElementRule(ARG_MODEL,
						new XMLSyntaxRule[]{new ElementRule(ARGModel.class)}),
				AttributeRule.newBooleanRule(ANCESTRAL_RESTRICTION,true),
		};

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			XMLObject cxo = xo.getChild(RECOMBINATION_RATE);
			Parameter rRate = (Parameter) cxo.getChild(Parameter.class);

			cxo = xo.getChild(POPULATION_SIZE);
			Parameter pSize = (Parameter) cxo.getChild(Parameter.class);

			cxo = xo.getChild(ARG_MODEL);
			ARGModel argModel = (ARGModel) cxo.getChild(ARGModel.class);

			int maxreassort = Integer.MAX_VALUE;
			if (xo.hasAttribute(MAX_REASSORTMENTS)) {
				maxreassort = xo.getIntegerAttribute(MAX_REASSORTMENTS);
			}
			
			boolean ancestral = false;
			if(xo.hasAttribute(ANCESTRAL_RESTRICTION))
				ancestral = xo.getBooleanAttribute(ANCESTRAL_RESTRICTION);

			return new ARGCoalescentLikelihood(pSize, rRate, argModel, false, maxreassort, ancestral);
		}


	};

	public String toString() {
		return getClass().getSimpleName() + " " + super.toString();
	}


}
