/*
 * NewCoalescentLikelihood.java
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

import dr.evolution.coalescent.Coalescent;
import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.tree.Tree;
import dr.evolution.util.*;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.xml.*;

import java.util.*;
import java.util.logging.Logger;


/**
 * A likelihood function for the coalescent. Takes a tree and a demographic model.
 *
 * Parts of this class were derived from C++ code provided by Oliver Pybus.
 *
 * @version $Id: NewCoalescentLikelihood.java,v 1.6 2005/05/24 20:25:57 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public final class CoalescentLikelihood extends AbstractCoalescentLikelihood implements Likelihood, Units {

	// PUBLIC STUFF

	public static final String COALESCENT_LIKELIHOOD = "coalescentLikelihood";
	public static final String MODEL = "model";
	public static final String POPULATION_TREE = "populationTree";

	public static final String INCLUDE = "include";
	public static final String EXCLUDE = "exclude";

	public CoalescentLikelihood(Tree tree,
	                            TaxonList includeSubtree,
	                            List<TaxonList> excludeSubtrees,
	                            DemographicModel demoModel) throws Tree.MissingTaxonException {

		super(COALESCENT_LIKELIHOOD, tree, includeSubtree, excludeSubtrees);

		this.demoModel = demoModel;

		addModel(demoModel);
	}

	// **************************************************************
	// Likelihood IMPLEMENTATION
	// **************************************************************

	/**
	 * Calculates the log likelihood of this set of coalescent intervals,
	 * given a demographic model.
	 */
	public double calculateLogLikelihood() {

		DemographicFunction demoFunction = demoModel.getDemographicFunction();

		double lnL =  Coalescent.calculateLogLikelihood(getIntervals(), demoFunction);

		if (Double.isNaN(lnL) || Double.isInfinite(lnL)) {
			Logger.getLogger("error").severe("CoalescentLikelihood is " + Double.toString(lnL));
		}

		return lnL;
	}

	// **************************************************************
	// Units IMPLEMENTATION
	// **************************************************************

	/**
	 * Sets the units these coalescent intervals are
	 * measured in.
	 */
	public final void setUnits(Type u)
	{
		demoModel.setUnits(u);
	}

	/**
	 * Returns the units these coalescent intervals are
	 * measured in.
	 */
	public final Type getUnits()
	{
		return demoModel.getUnits();
	}

	// ****************************************************************
	// Private and protected stuff
	// ****************************************************************

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return COALESCENT_LIKELIHOOD; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			XMLObject cxo = (XMLObject)xo.getChild(MODEL);
			DemographicModel demoModel = (DemographicModel)cxo.getChild(DemographicModel.class);

			cxo = (XMLObject)xo.getChild(POPULATION_TREE);
			TreeModel treeModel = (TreeModel)cxo.getChild(TreeModel.class);

			TaxonList includeSubtree = null;

			if (xo.hasSocket(INCLUDE)) {
				includeSubtree = (TaxonList)xo.getElementFirstChild(INCLUDE);
			}

			List<TaxonList> excludeSubtrees = new ArrayList<TaxonList>();

			if (xo.hasSocket(EXCLUDE)) {
				cxo = (XMLObject)xo.getChild(EXCLUDE);
				for (int i =0; i < cxo.getChildCount(); i++) {
					excludeSubtrees.add((TaxonList)cxo.getChild(i));
				}
			}

			try {
				return new CoalescentLikelihood(treeModel, includeSubtree, excludeSubtrees, demoModel);
			} catch (Tree.MissingTaxonException mte) {
				throw new XMLParseException("treeModel missing a taxon from taxon list in " + getParserName() + " element");
			}
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
				new ElementRule(MODEL, new XMLSyntaxRule[] {
						new ElementRule(DemographicModel.class)
				}, "The demographic model which describes the coalescent rate over time"),
				new ElementRule(POPULATION_TREE, new XMLSyntaxRule[] {
						new ElementRule(TreeModel.class)
				}, "The treeModel"),
				new ElementRule(INCLUDE, new XMLSyntaxRule[] {
						new ElementRule(Taxa.class)
				}, "An optional subset of taxa on which to calculate the likelihood (should be monophyletic)", true),
				new ElementRule(EXCLUDE, new XMLSyntaxRule[] {
						new ElementRule(Taxa.class, 1, Integer.MAX_VALUE)
				}, "One or more subsets of taxa which should be excluded from calculate the likelihood (should be monophyletic)", true)
		};
	};



	// ****************************************************************
	// Private and protected stuff
	// ****************************************************************


	/** The demographic model. */
	private DemographicModel demoModel = null;
}