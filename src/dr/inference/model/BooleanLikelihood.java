/*
 * BooleanLikelihood.java
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

package dr.inference.model;

import java.util.ArrayList;

/**
 * A class that returns a log likelihood of a set of boolean statistics.
 * If all the statistics are true then it returns 0.0 otherwise -INF.
 *
 * @author Andrew Rambaut
 *
 * @version $Id: BooleanLikelihood.java,v 1.8 2005/05/24 20:25:59 rambaut Exp $
 */
public class BooleanLikelihood extends Likelihood.Abstract {

	public BooleanLikelihood() {
		super(null);
	}

	/**
	 * Adds a statistic, this is the data for which the likelihood is calculated.
	 */
	public void addData(BooleanStatistic data) { dataList.add(data); }

	protected ArrayList<BooleanStatistic> dataList = new ArrayList<BooleanStatistic>();

	// **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

	/**
	 * Overridden to always return false.
	 */
	protected boolean getLikelihoodKnown() {
		return false;
	}

	/**
     * Calculate the log likelihood of the current state.
	 * If all the statistics are true then it returns 0.0 otherwise -INF.
     * @return the log likelihood.
     */
	public double calculateLogLikelihood() {

		if (getBooleanState()) {
			return Double.NEGATIVE_INFINITY;
		} else {
			return 0.0;
		}
	}

	public boolean getBooleanState() {

        for (BooleanStatistic statistic : dataList) {
            for (int j = 0; j < statistic.getDimension(); j++) {
                if (!statistic.getBoolean(j)) {
                    return true;
                }
            }
        }
        return false;
	}

    /**
     * Boolean likelihoods would generally want to be evaluated early to allow for a quick
     * termination of the evaluation in the case of a zero likelihood.
     * @return
     */
    public boolean evaluateEarly() {
        return true;
    }

}

