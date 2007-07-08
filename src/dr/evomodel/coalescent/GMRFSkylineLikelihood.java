/*
 * SkylineLikelihood.java
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

import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * A likelihood function for the Gaussian Markov random field population trajectory.
 * *
 *
 * @author Vladimir Minin
 * @author Marc Suchard
 * @version $Id: GMRFSkylineLikelihood.java,v 1.3 2007/03/20 22:40:04 msuchard Exp $
 */
public class GMRFSkylineLikelihood extends CoalescentLikelihood {

	// PUBLIC STUFF

	public static final String SKYLINE_LIKELIHOOD = "gmrfSkyLineLikelihood";
	public static final String POPULATION_PARAMETER = "populationSizes";
	public static final String PRECISION_PARAMETER = "precisionParameter";

	// PRIVATE STUFF

	private Parameter popSizeParameter;
	private Parameter precisionParameter;

	public GMRFSkylineLikelihood(Tree tree, Parameter popParameter, Parameter precParameter) {
		super(SKYLINE_LIKELIHOOD);

		this.popSizeParameter = popParameter;
		this.precisionParameter = precParameter;
		int tips = tree.getExternalNodeCount();
		int params = popSizeParameter.getDimension();
		if (tips - params != 1) {
			throw new IllegalArgumentException("Number of tips (" + tips + ") must be one greater than number of pop sizes (" + params + ")");
		}

		this.tree = tree;
		if (tree instanceof TreeModel) {
			addModel((TreeModel) tree);
		}
		addParameter(popSizeParameter);
		addParameter(precisionParameter);
		setupIntervals();

		addStatistic(new DeltaStatistic());
	}

	// **************************************************************
	// Likelihood IMPLEMENTATION
	// **************************************************************

	/**
	 * Calculates the log likelihood of this set of coalescent intervals,
	 * given a demographic model.
	 */
	public double calculateLogLikelihood() {

		if (!intervalsKnown) setupIntervals();

		double logL = 0.0;

		double currentTime = 0.0;

		int popIndex = 0;

		ConstantPopulation cp = new ConstantPopulation(Units.Type.YEARS);

		for (int j = 0; j < intervalCount; j++) {

			cp.setN0(popSizeParameter.getParameterValue(popIndex));
			if (getIntervalType(j) == CoalescentEventType.COALESCENT) {
				popIndex += 1;
			}

			//		logL += calculateIntervalLikelihood(cp, intervals[j], currentTime, lineageCounts[j], getIntervalType(j));

			// insert zero-length coalescent intervals
			int diff = getCoalescentEvents(j) - 1;
			for (int k = 0; k < diff; k++) {
				cp.setN0(popSizeParameter.getParameterValue(popIndex));
				logL += calculateIntervalLikelihood(cp, 0.0, currentTime, lineageCounts[j] - k - 1, CoalescentEventType.COALESCENT);
				popIndex += 1;
			}

			currentTime += intervals[j];


		}

		// Calculate GMRF density; here GMRF = RW(1)
		logL += -0.5 * calculateWeightedSSE();

		return logL;
	}

	public Parameter getPrecisionParameter() {
		return precisionParameter;
	}

	public double calculateWeightedSSE() {
		double weightedSSE = 0;
		double currentPopSize = popSizeParameter.getParameterValue(0); // todo: do we need a prior on N_e(0)?
		double currentInterval = intervals[0];
		for (int j = 1; j < intervalCount; j++) {
			double nextPopSize = popSizeParameter.getParameterValue(j);
			double nextInterval = intervals[j];
			double delta = nextPopSize - currentPopSize;
			double weight = (currentInterval + nextInterval) / 2.0;
			weightedSSE += delta * delta / weight;
			currentPopSize = nextPopSize;
			currentInterval = nextInterval;
		}
		return weightedSSE;

	}

	// ****************************************************************
	// Private and protected stuff
	// ****************************************************************

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return SKYLINE_LIKELIHOOD;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			XMLObject cxo = (XMLObject) xo.getChild(POPULATION_PARAMETER);
			Parameter popParameter = (Parameter) cxo.getChild(Parameter.class);

			cxo = (XMLObject) xo.getChild(PRECISION_PARAMETER);
			Parameter precParameter = (Parameter) cxo.getChild(Parameter.class);

			cxo = (XMLObject) xo.getChild(POPULATION_TREE);
			TreeModel treeModel = (TreeModel) cxo.getChild(TreeModel.class);

			return new GMRFSkylineLikelihood(treeModel, popParameter, precParameter);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element represents the likelihood of the tree given the population size vector.";
		}

		public Class getReturnType() {
			return Likelihood.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				new ElementRule(POPULATION_PARAMETER, new XMLSyntaxRule[]{
						new ElementRule(Parameter.class)
				}),
				new ElementRule(PRECISION_PARAMETER, new XMLSyntaxRule[]{
						new ElementRule(Parameter.class)
				}),
				new ElementRule(POPULATION_TREE, new XMLSyntaxRule[]{
						new ElementRule(TreeModel.class)
				}),
		};
	};


}