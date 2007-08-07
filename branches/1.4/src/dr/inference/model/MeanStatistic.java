/*
 * MeanStatistic.java
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

package dr.inference.model;

import dr.xml.*;

import java.util.Vector;

/**
 * @version $Id: MeanStatistic.java,v 1.9 2005/05/24 20:26:00 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class MeanStatistic extends Statistic.Abstract {

	public static String MEAN_STATISTIC = "meanStatistic";

	public MeanStatistic(String name) {
		super(name);
	}

	public void addStatistic(Statistic statistic) {
		statistics.add(statistic);
	}

	public int getDimension() { return 1; }

	/** @return mean of contained statistics */
	public double getStatisticValue(int dim) {
		double sum = 0.0;
		int dimensionCount = 0;
		int n;
		Statistic statistic;

		for (int i = 0; i < statistics.size(); i++) {
			statistic = (Statistic)statistics.get(i);
			n = statistic.getDimension();

			for (int j = 0; j < n; j++) {
				sum += statistic.getStatisticValue(j);
			}
			dimensionCount += n;
		}

		return sum / dimensionCount;
	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String[] getParserNames() { return new String[] { getParserName(), "mean" }; }
		public String getParserName() { return MEAN_STATISTIC; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			MeanStatistic meanStatistic = new MeanStatistic(MEAN_STATISTIC);

			for (int i =0; i < xo.getChildCount(); i++) {
				Object child = xo.getChild(i);
				if (child instanceof Statistic) {
					meanStatistic.addStatistic((Statistic)child);
				} else {
					throw new XMLParseException("Unknown element found in " + getParserName() + " element:" + child);
				}
			}

			return meanStatistic;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element returns a statistic that is the mean of the child statistics.";
		}

		public Class getReturnType() { return MeanStatistic.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new ElementRule(Statistic.class, 1, Integer.MAX_VALUE )
		};
	};


	// ****************************************************************
	// Private and protected stuff
	// ****************************************************************

	private Vector statistics = new Vector();
}
