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

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.IntervalList;
import dr.evolution.coalescent.IntervalType;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.demographicmodel.DemographicModel;
import dr.evomodelxml.coalescent.CoalescentLikelihoodParser;
import dr.math.Binomial;

import java.util.logging.Logger;


/**
 * A likelihood function for the coalescent. Takes an intervalList and a demographic model.
 * If the interval list is a model then it will listen for changes.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public final class CoalescentLikelihood extends AbstractCoalescentLikelihood implements Units {

	// PUBLIC STUFF

	/**
	 * A constructor that takes an IntervalList. This is uses the older DemographicModel which
	 * is no deprecated but left here for backwards compatibility
	 * @param intervalList the interval list
	 * @param demographicModel a demographic model
	 */
	public CoalescentLikelihood(IntervalList intervalList,
								DemographicModel demographicModel) {

		super(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, intervalList);

		this.populationSizeModel = null;
		this.demographicModel = demographicModel;

		addModel(demographicModel);
	}

	/**
	 * A constructor that takes an IntervalList
	 * @param intervalList the interval list
	 * @param populationSizeModel a population size model
	 */
	public CoalescentLikelihood(IntervalList intervalList,
								PopulationSizeModel populationSizeModel) {

		super(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, intervalList);

		this.populationSizeModel = populationSizeModel;
		this.demographicModel = null;

		addModel(populationSizeModel);
	}

	// **************************************************************
	// Likelihood IMPLEMENTATION
	// **************************************************************

	/**
	 * Calculates the log likelihood of this set of coalescent intervals,
	 * given a demographic model.
	 */
	protected double calculateLogLikelihood() {

		double lnL;

		if (populationSizeModel != null) {
			PopulationSizeFunction popFunction = populationSizeModel.getPopulationSizeFunction();

			lnL = calculateLogLikelihood(popFunction);

		} else {
			DemographicFunction demoFunction = demographicModel.getDemographicFunction();

			lnL = calculateLogLikelihood(demoFunction);
		}
		if (Double.isNaN(lnL) || Double.isInfinite(lnL)) {
			Logger.getLogger("warning").severe("CoalescentLikelihood for " + demographicModel.getId() + " is " + Double.toString(lnL));
		}

		return lnL;
	}

	/**
	 * Calculates the log likelihood of this set of coalescent intervals,
	 * given a demographic model.
	 */
	protected double calculateLogLikelihood(DemographicFunction demographicFunction) {

		double logL = 0.0;

		IntervalList intervals = getIntervalList();

		final int n = intervals.getIntervalCount();

		if (n == 0) {
			return 0.0;
		}

		double startTime = intervals.getStartTime();

		for (int i = 0; i < n; i++) {

			final double duration = intervals.getInterval(i);
			final double finishTime = startTime + duration;

			final double intervalArea = demographicFunction.getIntegral(startTime, finishTime);
			if( intervalArea == 0 && duration != 0 ) {
				return Double.NEGATIVE_INFINITY;
			}
			final int lineageCount = intervals.getLineageCount(i);


			final double kChoose2 = Binomial.choose2(lineageCount);
			// common part
			logL += -kChoose2 * intervalArea;

			if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {

				final double demographicAtCoalPoint = demographicFunction.getDemographic(finishTime);

				// if value at end is many orders of magnitude different than mean over interval reject the interval
				// This is protection against cases where ridiculous infinitesimal population size at the end of a
				// linear interval drive coalescent values to infinity.

				if( duration == 0.0 || demographicAtCoalPoint * (intervalArea/duration) >= demographicFunction.getThreshold() ) {
					//                if( duration == 0.0 || demographicAtCoalPoint >= threshold * (duration/intervalArea) ) {
					logL -= Math.log(demographicAtCoalPoint);
				} else {
					// remove this at some stage
					//  System.err.println("Warning: " + i + " " + demographicAtCoalPoint + " " + (intervalArea/duration) );
					return Double.NEGATIVE_INFINITY;
				}

			}

			startTime = finishTime;
		}

		return logL;
	}

	protected double calculateLogLikelihood(PopulationSizeFunction populationSizeFunction) {

		double logL = 0.0;

		IntervalList intervals = getIntervalList();

		final int n = intervals.getIntervalCount();

		if (n == 0) {
			return 0.0;
		}

		double startTime = intervals.getStartTime();

		for (int i = 0; i < n; i++) {

			final double duration = intervals.getInterval(i);
			final double finishTime = startTime + duration;

			final double intervalArea = populationSizeFunction.getIntegral(startTime, finishTime);
			if( intervalArea == 0 && duration != 0 ) {
				return Double.NEGATIVE_INFINITY;
			}
			final int lineageCount = intervals.getLineageCount(i);


			final double kChoose2 = Binomial.choose2(lineageCount);

			// common part
			logL += -kChoose2 * intervalArea;

			if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {

				logL -= populationSizeFunction.getLogDemographic(finishTime);
			}

			startTime = finishTime;
		}

		return logL;
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
		demographicModel.setUnits(u);
	}

	/**
	 * Returns the units these coalescent intervals are
	 * measured in.
	 */
	public final Type getUnits()
	{
		return demographicModel.getUnits();
	}

	// ****************************************************************
	// Private and protected stuff
	// ****************************************************************

	/** the population size model */
	private final PopulationSizeModel populationSizeModel;

	/** The demographic model. */
	private final DemographicModel demographicModel;
}