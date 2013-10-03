/*
 * PriorParsers.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import jebl.math.Binomial;
import dr.evolution.coalescent.TreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;

/**
 * Calculates the coalescent constant for a given tree(Model).
 *
 * @author Guy Baele
 */
public final class CoalescentConstantLikelihood extends Likelihood.Abstract {

	private TreeModel treeModel;
	private TreeIntervals intervals;
	
	// PUBLIC STUFF
	public CoalescentConstantLikelihood(TreeModel treeModel) {

		super(treeModel);
		this.treeModel = treeModel;
		this.intervals = new TreeIntervals(treeModel);

	}

    // **************************************************************
	// Likelihood IMPLEMENTATION
	// **************************************************************

	/**
	 * Calculates the coalescent constant given the treeModel.
	 */
	public double calculateLogLikelihood() {

		intervals.setIntervalsUnknown();
		final int nIntervals = intervals.getIntervalCount();
		//System.err.println("Interval count: " + nIntervals);
		double logPDF = 0.0;
        for (int i = 0; i < nIntervals; i++) {
        	//System.err.println("Lineage count " + i + ": " + intervals.getLineageCount(i));
        	if (intervals.getLineageCount(i) > 2) {
        		logPDF += Math.log(Binomial.choose2(intervals.getLineageCount(i)));
        	}
        }
        //System.err.println("logPDF = " + (-logPDF));
        return -logPDF;
        
	}

}