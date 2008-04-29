/*
 * LogarithmStatistic.java
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

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: ExponentialStatistic.java,v 1.3 2005/05/24 20:26:00 rambaut Exp $
 */
public class LogarithmStatistic extends Statistic.Abstract {

	public static String LOGARITHM_STATISTIC = "logarithmStatistic";
	public static String BASE = "base";

	private final Statistic statistic;
	private final double base;

	public LogarithmStatistic(String name, Statistic statistic, double base) {
		super(name);
		this.statistic = statistic;
		this.base = base;
	}

	public int getDimension() {
		return statistic.getDimension();
	}

	/**
	 * @return mean of contained statistics
	 */
	public double getStatisticValue(int dim) {
		if (base <= 1.0) {
			return Math.log(statistic.getStatisticValue(dim));
		} else {
			return Math.log(statistic.getStatisticValue(dim)) / Math.log(base);
		}
	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String[] getParserNames() {
			return new String[]{getParserName(), "logarithm"};
		}

		public String getParserName() {
			return LOGARITHM_STATISTIC;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			LogarithmStatistic logStatistic;

			double base = 0.0; // base 0.0 means natural logarithm
			if (xo.hasAttribute(BASE)) {
				base = xo.getDoubleAttribute(BASE);
			}

			if (base <= 1.0 && base != 0.0) {
				throw new XMLParseException("Error parsing " + getParserName() + " element: base attribute should be > 1");
			}

			Object child = xo.getChild(0);
			if (child instanceof Statistic) {
				logStatistic = new LogarithmStatistic(LOGARITHM_STATISTIC, (Statistic) child, base);
			} else {
				throw new XMLParseException("Unknown element found in " + getParserName() + " element:" + child);
			}

			return logStatistic;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element returns a statistic that is the element-wise natural logarithm of the child statistic.";
		}

		public Class getReturnType() {
			return ExponentialStatistic.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				AttributeRule.newDoubleRule(BASE, true, "An optional base for the logarithm (default is the natural logarithm, base e)"),
				new ElementRule(Statistic.class, 1, 1)
		};
	};
}