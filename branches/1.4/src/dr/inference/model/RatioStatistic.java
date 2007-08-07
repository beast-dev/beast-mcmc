/*
 * ProductStatistic.java
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
 * @version $Id: ProductStatistic.java,v 1.2 2005/05/24 20:26:00 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class RatioStatistic extends Statistic.Abstract {

	public static String RATIO_STATISTIC = "ratioStatistic";

	public RatioStatistic(String name, Statistic numerator, Statistic denominator) {
		super(name);

        this.numerator = numerator;
        this.denominator = denominator;
        
        if (denominator.getDimension() != 1 &&
                numerator.getDimension() != 1 &&
                denominator.getDimension() != numerator.getDimension()) {
            throw new IllegalArgumentException();
        }

        if (denominator.getDimension() == 1) {
            dimension = numerator.getDimension();
        } else {
            dimension = denominator.getDimension();
        }
    }

    public int getDimension() { return dimension; }

	/** @return mean of contained statistics */
	public double getStatisticValue(int dim) {

        if (numerator.getDimension() == 1) {
            return numerator.getStatisticValue(0) / denominator.getStatisticValue(dim);
        } else if (denominator.getDimension() == 1) {
            return numerator.getStatisticValue(dim) / denominator.getStatisticValue(0);
        } else {
            return numerator.getStatisticValue(dim) / denominator.getStatisticValue(dim);
        }
    }

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String[] getParserNames() { return new String[] { getParserName(), "ratio" }; }
		public String getParserName() { return RATIO_STATISTIC; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Statistic numerator = (Statistic)xo.getChild(0);
            Statistic denominator = (Statistic)xo.getChild(1);

            RatioStatistic ratioStatistic;
            try {
                ratioStatistic = new RatioStatistic(RATIO_STATISTIC, numerator, denominator);
            } catch (IllegalArgumentException iae) {
                throw new XMLParseException("Error parsing " + getParserName() + " the numerator and denominator statistics " +
                        "should be of the same dimension or of dimension 1");
            }
			return ratioStatistic;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element returns a statistic that is the ratio of the 2 child statistics.";
		}

		public Class getReturnType() { return ProductStatistic.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
                new ElementRule(Statistic.class, "The numerator statistic" ),
                new ElementRule(Statistic.class, "The denominator statistic" )
		};
	};


	// ****************************************************************
	// Private and protected stuff
	// ****************************************************************

    private int dimension = 0;
    private Statistic numerator = null;
	private Statistic denominator = null;
}