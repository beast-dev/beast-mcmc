/*
 * SkylineLikelihood.java
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

import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.CoalescentLikelihoodParser;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.*;

import java.util.Collections;
import java.util.List;

/**
 * A likelihood function for the coalescent. Takes a tree and a exponential markov model.
 * *
 * @version $Id: SkylineLikelihood.java,v 1.3 2004/10/01 22:40:04 alexei Exp $
 *
 * @author Alexei Drummond
 */
public class SkylineLikelihood extends OldAbstractCoalescentLikelihood implements Citable {

	// PUBLIC STUFF

	public static final String SKYLINE_LIKELIHOOD = "skyLineLikelihood";
	public static final String POPULATION_SIZES = "populationSizes";

	public SkylineLikelihood(Tree tree, Parameter popSizeParameter) {
		super(SKYLINE_LIKELIHOOD);

		this.popSizeParameter = popSizeParameter;
		int tips = tree.getExternalNodeCount();
		int params = popSizeParameter.getDimension();
		if (tips - params != 1) {
			throw new IllegalArgumentException("Number of tips (" + tips + ") must be one greater than number of pop sizes (" + params + ")");
		}

		this.tree = tree;
		if (tree instanceof TreeModel) {
			addModel((TreeModel)tree);
		}
		addVariable(popSizeParameter);
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

		int popIndex=0;

		ConstantPopulation cp = new ConstantPopulation(Units.Type.YEARS);

		for (int j = 0; j < intervalCount; j++) {

			cp.setN0(popSizeParameter.getParameterValue(popIndex));
			if (getIntervalType(j) == CoalescentEventType.COALESCENT) {
				popIndex += 1;
			}

			logL += calculateIntervalLikelihood(cp, intervals[j], currentTime, lineageCounts[j], getIntervalType(j));

			// insert zero-length coalescent intervals
			int diff = getCoalescentEvents(j)-1;
			for (int k = 0; k < diff; k++) {
				cp.setN0(popSizeParameter.getParameterValue(popIndex));
				logL += calculateIntervalLikelihood(cp, 0.0, currentTime, lineageCounts[j]-k-1,
                        CoalescentEventType.COALESCENT);
				popIndex += 1;
			}

			currentTime += intervals[j];


		}

		return logL;
	}

	// ****************************************************************
	// Private and protected stuff
	// ****************************************************************

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return SKYLINE_LIKELIHOOD; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			XMLObject cxo = (XMLObject)xo.getChild(POPULATION_SIZES);
			Parameter param = (Parameter)cxo.getChild(Parameter.class);

			cxo = (XMLObject)xo.getChild(CoalescentLikelihoodParser.POPULATION_TREE);
			TreeModel treeModel = (TreeModel)cxo.getChild(TreeModel.class);

			return new SkylineLikelihood(treeModel, param);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element represents the likelihood of the tree given the population size vector.";
		}

		public Class getReturnType() { return Likelihood.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new ElementRule(POPULATION_SIZES, new XMLSyntaxRule[] {
				new ElementRule(Parameter.class)
			}),
			new ElementRule(CoalescentLikelihoodParser.POPULATION_TREE, new XMLSyntaxRule[] {
				new ElementRule(TreeModel.class)
			}),
		};
	};

	/** The demographic model. */
	Parameter popSizeParameter = null;

	@Override
	public Citation.Category getCategory() {
		return Citation.Category.TREE_PRIORS;
	}

	@Override
	public String getDescription() {
		return "Bayesian Skyline Coalescent";
	}

	@Override
	public List<Citation> getCitations() {
		return Collections.singletonList(CITATION);
	}

	public static Citation CITATION = new Citation(
			new Author[]{
					new Author("AJ", "Drummond"),
					new Author("A", "Rambaut"),
					new Author("B", "Shapiro"),
					new Author("OG", "Pybus")
			},
			"Bayesian coalescent inference of past population dynamics from molecular sequences",
			2005,
			"Mol Biol Evol",
			22, 1185, 1192
	);
}