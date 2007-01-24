/*
 * EmpiricalDistributionLikelihood.java
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

package dr.inference.distribution;

import dr.inference.model.Likelihood;
import dr.inference.model.Statistic;
import dr.inference.trace.Trace;
import dr.util.HeapSort;
import dr.xml.*;

import java.io.FileReader;
import java.io.IOException;

/**
 * A class that returns the log likelihood of a set of data (statistics)
 * being distributed according to the given empirical distribution.
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: EmpiricalDistributionLikelihood.java,v 1.8 2006/08/08 08:43:11 rambaut Exp $
 */

public class EmpiricalDistributionLikelihood extends AbstractDistributionLikelihood {

	public static final String EMPIRICAL_DISTRIBUTION_LIKELIHOOD = "empiricalDistributionLikelihood";

	public static final String LOG_FILE = "logFile";

	public static final String FILE_NAME = "fileName";
	public static final String STATISTIC = "statistic";
	public static final String BURN_IN = "burnIn";
	public static final String DELTA = "delta";

	public static final String HISTOGRAM = "histogram";

	public static final String DATA = "data";

	Trace trace = null;
	double[] values = null;
	int[] indices = null;
	double halfDx = 0.1;

	/**
	 * Constructor that tries to guess a value of dx
	 */
	public EmpiricalDistributionLikelihood(String fileName, String statName, int burnIn) {

		super(null);

		try {
			FileReader fileReader = new FileReader(fileName);

			trace = Trace.Utils.loadTrace(fileReader, statName);

			if (trace == null) throw new IllegalArgumentException("Statistic " + statName + " does not exist in file " + fileName);

			values = trace.getValues(burnIn);
			indices = new int[values.length];
			HeapSort.sort(values, indices);

			if (values.length < 100) {
				throw new RuntimeException("Not enough values to construct reliable empirical distribution!");
			}

			double proportion = 0.01;
			double suggestedDx = getDxForMaxPValue(proportion);

			System.out.println("Estimated dx = " + suggestedDx);
			System.out.println("   * This would allow a maximum of " + (int)Math.round(values.length * proportion) + " probability density classes.");
			System.out.println("   * The number of classes is also the largest finite probability ratio.");

			checkForZeroDensityRegions();

			this.halfDx = suggestedDx / 2.0;
		} catch (IOException ioe) {
			throw new RuntimeException("File not found or somesuch!");
		}
	}

	public EmpiricalDistributionLikelihood(String fileName, String statName, int burnIn, double dx) {

		super(null);

		try {
			FileReader fileReader = new FileReader(fileName);

			trace = Trace.Utils.loadTrace(fileReader, statName);

			if (trace == null) throw new IllegalArgumentException("Statistic " + statName + " does not exist in file " + fileName);

			values = trace.getValues(burnIn);
			indices = new int[values.length];
			HeapSort.sort(values, indices);

			if (values.length < 100) {
				throw new RuntimeException("Not enough values to construct reliable empirical distribution!");
			}

			double proportion = 0.01;
			double suggestedDx = getDxForMaxPValue(proportion);
			System.out.println("Suggested dx = " + suggestedDx);
			System.out.println("   * This would allow a maximum of " + (int)Math.round(values.length * proportion) + " probability density classes.");
			System.out.println("   * The number of classes is also the largest finite probability ratio.");
			System.out.println("Your dx = " + dx);
			checkForZeroDensityRegions();

			this.halfDx = dx / 2.0;
		} catch (IOException ioe) {
			throw new RuntimeException("File not found or somesuch!");
		}
	}

	private void checkForZeroDensityRegions() {

		double dx = 2.0 * halfDx;
		for (int i =0; i < (values.length-1); i++) {
			if (values[i+1] - values[i] > dx) {
				System.out.println("Warning: Zero probability density region between " + values[i] + " and " + values[i+1]);
			}
		}
	}

	/**
	 * @return a value for dx that will contain a maximum of the given proportion of points.
	 */
	private double getDxForMaxPValue(double maxPValue) {

		double minRange = Double.MAX_VALUE;
		int diff = (int)Math.round(maxPValue * (double)values.length);
		for (int i =0; i <= (values.length - diff); i++) {
			double minValue = values[indices[i]];
			double maxValue = values[indices[i+diff-1]];
			double range = Math.abs(maxValue - minValue);
			if (range < minRange) {
				minRange = range;
			}
		}
		return minRange;
	}

	/**
	 * @return the probability density in a small area.
	 */
	private double getDensity(double minValue, double maxValue) {
		double density = (((double)getCount(minValue, maxValue)) / values.length) / (maxValue-minValue);
		return density;
	}

	/**
	 * @return the number of values in the given range
	 */
	private int getCount(double minValue, double maxValue) {

		//WARNING! EXTREMELY INEFFICIENT
		// SHOULD USE BINARY SEARCH!

		int count = 0;
		for (int i =0; i < values.length; i++) {
			if (values[i] >= minValue && values[i] <= maxValue) count += 1;
		}
		return count;
	}

	// **************************************************************
	// Likelihood IMPLEMENTATION
	// **************************************************************

	/**
	 * Calculate the log likelihood of the current state.
	 * @return the log likelihood.
	 */
	public double calculateLogLikelihood() {

		double value, logL = 0.0;

        for (Statistic statistic : dataList) {
            for (int j = 0; j < statistic.getDimension(); j++) {
                value = statistic.getStatisticValue(j);
                logL += Math.log(getDensity(value - halfDx, value + halfDx));
            }
        }
        return logL;
	}

	/**
	 * Reads a distribution likelihood from a DOM Document element.
	 */
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return EMPIRICAL_DISTRIBUTION_LIKELIHOOD; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			AbstractDistributionLikelihood likelihood = null;

			if (xo.getChild(LOG_FILE) != null) {
				XMLObject cxo = (XMLObject)xo.getChild(LOG_FILE);
				String fileName = xo.getStringAttribute(FILE_NAME);
				String statName = xo.getStringAttribute(STATISTIC);
				int burnIn = xo.getIntegerAttribute(BURN_IN);
				if (xo.hasAttribute(DELTA)) {

					double delta = xo.getDoubleAttribute(DELTA);
					likelihood = new EmpiricalDistributionLikelihood(fileName, statName, burnIn, delta);
				} else {
					likelihood = new EmpiricalDistributionLikelihood(fileName, statName, burnIn);
				}
			} else {
				XMLObject cxo = (XMLObject)xo.getChild(HISTOGRAM);

				throw new RuntimeException("Still to be implemented....");
			}

			XMLObject cxo = (XMLObject)xo.getChild(DATA);

			for (int i = 0; i < xo.getChildCount(); i++) {

				if (cxo.getChildCount() == 0) {
					throw new XMLParseException("no statistic element in " + cxo.getName() + " element");
				}

				for (int j = 0; j < cxo.getChildCount(); j++) {
					if (cxo.getChild(j) instanceof Statistic) {

						likelihood.addData( (Statistic)cxo.getChild(j));
					} else {
						throw new XMLParseException("illegal element in " + cxo.getName() + " element");
					}
				}
			}

			return likelihood;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
				new XORRule(
						new ElementRule(LOG_FILE, new XMLSyntaxRule[] {
								new StringAttributeRule(FILE_NAME, "The file name of an BEAST log file"),
								new StringAttributeRule(STATISTIC, "The name of the column in the log file that will be used as an empirical distribution"),
								AttributeRule.newIntegerRule(BURN_IN),
								AttributeRule.newBooleanRule(DELTA, true)
						}),
						new ElementRule(HISTOGRAM, new XMLSyntaxRule[] {
								new ContentRule("X,Y data describing the distribution" )})
				),
				new ElementRule(DATA, new XMLSyntaxRule[] { new ElementRule(Statistic.class, 1, Integer.MAX_VALUE )})
		};

		public String getParserDescription() {
			return "Calculates the likelihood of some data given some parametric or empirical distribution.";
		}

		public Class getReturnType() { return Likelihood.class; }
	};
}

