/*
 * CoalescentLikelihood.java
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
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodelxml.coalescent.CoalescentLikelihoodParser;

import java.util.List;
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
public final class CoalescentLikelihood extends AbstractCoalescentLikelihood implements Units {

	// PUBLIC STUFF
	public CoalescentLikelihood(Tree tree,
	                            TaxonList includeSubtree,
	                            List<TaxonList> excludeSubtrees,
	                            DemographicModel demoModel) throws TreeUtils.MissingTaxonException {

		super(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, tree, includeSubtree, excludeSubtrees);

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

		//double lnL =  Coalescent.calculateLogLikelihood(getIntervals(), demoFunction);
        double lnL =  Coalescent.calculateLogLikelihood(getIntervals(), demoFunction, demoFunction.getThreshold());

		if (Double.isNaN(lnL) || Double.isInfinite(lnL)) {
			Logger.getLogger("warning").severe("CoalescentLikelihood for " + demoModel.getId() + " is " + Double.toString(lnL));
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

	/** The demographic model. */
	private DemographicModel demoModel = null;
}