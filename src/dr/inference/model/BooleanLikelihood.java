/*
 * BooleanLikelihood.java
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

	public static final String BOOLEAN_LIKELIHOOD = "booleanLikelihood";

	public static final String DATA = "data";

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
	 * Reads a distribution likelihood from a DOM Document element.
	 */
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return BOOLEAN_LIKELIHOOD; }

		public Object parseXMLObject(XMLObject xo) {

			BooleanLikelihood likelihood = new BooleanLikelihood();

			for (int i = 0; i < xo.getChildCount(); i++) {
				if (xo.getChild(i) instanceof BooleanStatistic) {
					likelihood.addData( (BooleanStatistic)xo.getChild(i));
				}
			}

			return likelihood;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "A function that log likelihood of a set of boolean statistics. "+
					"If all the statistics are true then it returns 0.0 otherwise -infinity.";
		}

		public Class getReturnType() { return BooleanLikelihood.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new ElementRule(BooleanStatistic.class, 1, Integer.MAX_VALUE )
		};

	};
}

