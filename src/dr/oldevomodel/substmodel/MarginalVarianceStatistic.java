/*
 * MarginalVarianceStatistic.java
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

package dr.oldevomodel.substmodel;

import dr.inferencexml.distribution.MultivariateOUModel;
import dr.xml.*;
import dr.inference.model.Statistic;

/**
 * @author Marc A. Suchard
 */
public class MarginalVarianceStatistic extends Statistic.Abstract {

	public static final String VARIANCE_STATISTIC = "marginalVariance";

	public MarginalVarianceStatistic(MultivariateOUModel mvou) {
		this.mvou = mvou;
	}

	public int getDimension() {
		return mvou.getDimension();
	}

	public double getStatisticValue(int dim) {
		return mvou.getStatisticValue(dim);
	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return VARIANCE_STATISTIC;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			MultivariateOUModel mvou = (MultivariateOUModel) xo.getChild(MultivariateOUModel.class);
			return new MarginalVarianceStatistic(mvou);

		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element returns a statistic that is the matrix inverse of the child statistic.";
		}

		public Class getReturnType() {
			return MarginalVarianceStatistic.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				new ElementRule(MultivariateOUModel.class)
		};
	};

	private MultivariateOUModel mvou;
}
