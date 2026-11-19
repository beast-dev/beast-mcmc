/*
 * ExponentialProductLikelihood.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.coalescent.workingprior;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;

/**
 * Calculates a product of exponential densities and exponential tail probabilities: an exponential
 * density for each interval that ends in a coalescent event and a tail probability for each that
 * ends in a sampling event. The log population size is either the same throughout or changes at
 * each coalescent event (the posterior means of the coalescentEventsStatistic).
 *
 * @author Guy Baele
 */
public class ExponentialProductLikelihood extends Likelihood.Abstract {

	private final TreeModel treeModel;
	// the log population size before each coalescent event, or one value for them all
	private final double[] logPopSizes;

	/**
	 * @param logPopSize the log population size, used for every interval
	 */
	public ExponentialProductLikelihood(TreeModel treeModel, double logPopSize) {
		this(treeModel, new double[] { logPopSize });
	}

	/**
	 * @param logPopSizes the log population size before each coalescent event (use in combination
	 *                    with the coalescentEventsStatistic), or a single value
	 */
	public ExponentialProductLikelihood(TreeModel treeModel, double[] logPopSizes) {
		super(treeModel);
		this.treeModel = treeModel;
		this.logPopSizes = logPopSizes;
	}

	public double calculateLogLikelihood() {

		double logPDF = 0.0;

		CoalescentTreeIntervalStatistic ctis = new CoalescentTreeIntervalStatistic(treeModel);
		int coalescentIntervalCounter = 0;

		for (int i = 0; i < ctis.getDimension(); i++) {

			int combinations = (int)ctis.getLineageCount(i)*((int)ctis.getLineageCount(i)-1)/2;
			double branchLength = ctis.getStatisticValue(i);

			//single-lineage intervals are not counted
			if (ctis.getLineageCount(i) != 1) {

				final double logPopSize = logPopSizes.length == 1 ? logPopSizes[0] : logPopSizes[coalescentIntervalCounter];

				if (i == ctis.getDimension()-1) {
					//coalescent event at root: exponential density
					logPDF += -logPopSize - combinations*branchLength*Math.exp(-logPopSize);
				} else if (ctis.getLineageCount(i) > ctis.getLineageCount(i+1)) {
					//coalescent event: exponential density
					logPDF += -logPopSize - combinations*branchLength*Math.exp(-logPopSize);
					coalescentIntervalCounter++;
				} else {
					//sampling event: exponential tail probability
					logPDF += -combinations*branchLength*Math.exp(-logPopSize);
				}

			}

		}

		return logPDF;

	}

	/**
	 * Overridden to always return false.
	 */
	protected boolean getLikelihoodKnown() {
		return false;
	}

}
