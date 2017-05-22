/*
 * SiteLogLikelihoodLogger.java
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

package dr.app.beagle.tools;

import dr.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.xml.Reportable;

public class SiteLogLikelihoodLogger implements Loggable, Reportable {

	BeagleTreeLikelihood beagleTreeLikelihood;
	int patternCount;
	private SiteLogLikelihoodColumn[] columns = null;

	public SiteLogLikelihoodLogger(BeagleTreeLikelihood beagleTreeLikelihood) {
		this.beagleTreeLikelihood = beagleTreeLikelihood;

		patternCount = beagleTreeLikelihood.getPatternCount();

	}// END: Constructor

	@Override
	public LogColumn[] getColumns() {

		if (columns == null) {
			columns = new SiteLogLikelihoodColumn[patternCount];
			for (int site = 0; site < patternCount; site++) {

				columns[site] = new SiteLogLikelihoodColumn(site);

			}
		}

		return columns;
	}// END: getColumns

	private double getSiteLogLikelihood(int site) {
		double[] siteLikelihoods = beagleTreeLikelihood.getSiteLogLikelihoods();
		return siteLikelihoods[site];
	}// END: getSiteLogLikelihoods

	private class SiteLogLikelihoodColumn extends NumberColumn {

		final int site;

		public SiteLogLikelihoodColumn(int site) {
			super("SiteLogLikelihoodColumn");
			this.site = site;
		}

		@Override
		public double getDoubleValue() {
			return getSiteLogLikelihood(site);
		}

	}// END: SiteLogLikelihoodColumn class

	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (int site = 0; site < patternCount; ++site) {
			if (site > 0) {
				sb.append(", ");
			}

			sb.append(getColumns()[site].getFormatted());
		}

		return sb.toString();
	}// END: toString

	@Override
	public String getReport() {
		return toString();
	}// END: getReport

}// END: class
