/*
 * CovariateGMRFSkylineLikelihood.java
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.logging.Logger;

/**
 * @author Erik Bloomquist
 * @author Marc A. Suchard
 */
public class CovariateGMRFSkylineLikelihood extends GMRFSkyrideLikelihood {


	private Parameter covariateData;
	private Parameter covariateTimes;

	private ArrayList<CoalescentIntervalWithData> intervals;
	private ArrayList<CoalescentIntervalWithData> storedIntervals;


	public CovariateGMRFSkylineLikelihood(Tree tree, Parameter popParameter, Parameter precParameter,
	                                      Parameter lambda, Parameter beta, MatrixParameter dMatrix,
	                                      Parameter data, Parameter times) {
		super(tree, popParameter, null, precParameter, lambda, beta, dMatrix, false, true);

		covariateData = data;
		covariateTimes = times;

		fieldLength += covariateData.getDimension();

		addVariable(covariateData); // this can have missing values for imputation

	}

	//	@Override
	public void sSetupIntervals() {

		intervals.clear();
		intervals.ensureCapacity(fieldLength);

		NodeRef x;
		for (int i = 0; i < tree.getInternalNodeCount(); i++) {
			x = tree.getInternalNode(i);
			intervals.add(new CoalescentIntervalWithData(tree.getNodeHeight(x), Double.NaN, CoalescentEventType.COALESCENT));
		}
		for (int i = 0; i < tree.getExternalNodeCount(); i++) {
			x = tree.getExternalNode(i);
			if (tree.getNodeHeight(x) > 0.0) {
				intervals.add(new CoalescentIntervalWithData(tree.getNodeHeight(x), Double.NaN, CoalescentEventType.NEW_SAMPLE));
			}
		}
		for (int i = 0; i < covariateTimes.getDimension(); i++) {
			intervals.add(new CoalescentIntervalWithData(covariateTimes.getParameterValue(i),
					covariateData.getParameterValue(i), CoalescentEventType.NOTHING));
		}
		dr.util.HeapSort.sort(intervals);


		double a;
		for (int i = 0; i < intervals.size() - 1; i++) {
			a = intervals.get(i).length;
			intervals.get(i).length = a - intervals.get(i + 1).length;
		}
		intervals.remove(intervals.size() - 1);
		intervalsKnown = true;

	}

	public void setupGMRFWeights() {
		super.setupGMRFWeights();
	}

	public void storeState() {
		storedIntervals = new ArrayList<CoalescentIntervalWithData>(intervals.size());
		for (CoalescentIntervalWithData interval : intervals) {
			storedIntervals.add(interval.clone());
		}
	}

	public void restoreState() {
		intervals = storedIntervals;
		storedIntervals.clear();
	}

	private class CoalescentIntervalWithData implements Comparable<CoalescentIntervalWithData>, Cloneable {
		public final CoalescentEventType type;
		public double length;
		public final double datum;

		public CoalescentIntervalWithData(double length, double datum, CoalescentEventType type) {
			this.length = length;
			this.type = type;
			this.datum = datum;
		}

		public int compareTo(CoalescentIntervalWithData a) {

			if (a.length < this.length) {
				return -1;
			} else if (a.length == this.length) {
				Logger.getLogger("dr.evomodel.coalescent").severe("The current model " +
						"has 2 internal nodes or 1 node and 1 covariate at the same height");
				return 0;
			}
			return 1;
		}

		public String toString() {
			return "(" + length + "," + type + "," + datum + ")";
		}

		public CoalescentIntervalWithData clone() {
			return new CoalescentIntervalWithData(length, datum, type);
		}

	}

}
