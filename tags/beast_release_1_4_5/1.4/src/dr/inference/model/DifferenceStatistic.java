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
public class DifferenceStatistic extends Statistic.Abstract {

	public static String DIFFERENCE_STATISTIC = "differenceStatistic";

	public DifferenceStatistic(String name, Statistic term1, Statistic term2) {
		super(name);

        this.term1 = term1;
        this.term2 = term2;

        if (term1.getDimension() != 1 &&
                term2.getDimension() != 1 &&
                term1.getDimension() != term2.getDimension()) {
            throw new IllegalArgumentException();
        }

        if (term2.getDimension() == 1) {
            dimension = term1.getDimension();
        } else {
            dimension = term2.getDimension();
        }
    }

    public int getDimension() { return dimension; }

	/** @return mean of contained statistics */
	public double getStatisticValue(int dim) {

        if (term1.getDimension() == 1) {
            return term1.getStatisticValue(0) - term2.getStatisticValue(dim);
        } else if (term2.getDimension() == 1) {
            return term1.getStatisticValue(dim) - term1.getStatisticValue(0);
        } else {
            return term1.getStatisticValue(dim) - term2.getStatisticValue(dim);
        }
    }

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String[] getParserNames() { return new String[] { getParserName(), "difference" }; }
		public String getParserName() { return DIFFERENCE_STATISTIC; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Statistic term1 = (Statistic)xo.getChild(0);
            Statistic term2 = (Statistic)xo.getChild(1);

            DifferenceStatistic differenceStatistic;
            try {
                differenceStatistic = new DifferenceStatistic(DIFFERENCE_STATISTIC, term1, term2);
            } catch (IllegalArgumentException iae) {
                throw new XMLParseException("Error parsing " + getParserName() + " the left and right statistics " +
                        "should be of the same dimension or of dimension 1");
            }
			return differenceStatistic;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element returns a statistic that is the difference of the 2 child statistics.";
		}

		public Class getReturnType() { return ProductStatistic.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
                new ElementRule(Statistic.class, "The first statistic" ),
                new ElementRule(Statistic.class, "The second statistic" )
		};
	};


	// ****************************************************************
	// Private and protected stuff
	// ****************************************************************

    private int dimension = 0;
    private Statistic term1 = null;
	private Statistic term2 = null;
}