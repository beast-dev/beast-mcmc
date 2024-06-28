/*
 * DummyLikelihood.java
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

package dr.inference.model;

/**
 * A class that always returns a log likelihood of 0 but contains models that would otherwise
 * be unregistered with the MCMC. This is an ugly solution to a rare problem.
 *
 * @author Andrew Rambaut
 *
 */
public class DummyLikelihood extends Likelihood.Abstract {

	public DummyLikelihood(Model model) {
		super(model);
	}

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
		return 0.0;
	}

}

