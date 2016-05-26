/*
 * CorrelationStatistic.java
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

package dr.inference.model;

import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.xml.*;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jul 7, 2007
 * Time: 8:15:31 AM
 * To change this template use File | Settings | File Templates.
 */
public class CorrelationStatistic extends Statistic.Abstract {

//	public class CorrelationStatistic extends Statistic.Abstract {

	public static final String CORRELATION_STATISTIC = "correlation";
	public static final String DIMENSION1 = "dimension1";
	public static final String DIMENSION2 = "dimension2";

	private MatrixParameter precision = null;
	private int dim1;
	private int dim2;

	public CorrelationStatistic(String name, MatrixParameter precision,
	                            int dim1, int dim2) {
		super(name);
		this.precision = precision;
		this.dim1 = dim1 - 1;
		this.dim2 = dim2 - 1;
//			System.err.println("MAKE!");
//			System.exit(0);
	}

	public int getDimension() {
		return 1;
	}

	public double getStatisticValue(int dim) {

		double[][] variance = new SymmetricMatrix(
				precision.getParameterAsMatrix()).inverse().toComponents();
		return variance[dim1][dim2] /
				Math.sqrt(variance[dim1][dim1] * variance[dim2][dim2]);

	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return CORRELATION_STATISTIC;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

//				ReciprocalStatistic recipStatistic = null;
			MatrixParameter matrix = (MatrixParameter) xo.getChild(MatrixParameter.class);
			int dim1 = xo.getIntegerAttribute(DIMENSION1);
			int dim2 = xo.getIntegerAttribute(DIMENSION2);

			if (dim1 < 1 || dim1 > matrix.getRowDimension() || dim2 < 1 || dim2 > matrix.getColumnDimension())
				throw new XMLParseException("Invalid dimensions in " + getParserName() + " element");

			return new CorrelationStatistic(CORRELATION_STATISTIC, matrix, dim1, dim2);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element returns a precision that is the element-wise reciprocal of the child precision.";
		}

		public Class getReturnType() {
			return ReciprocalStatistic.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				new ElementRule(MatrixParameter.class, 1, 1),
				AttributeRule.newIntegerRule(DIMENSION1),
				AttributeRule.newIntegerRule(DIMENSION2)

		};
	};
//	}

}
