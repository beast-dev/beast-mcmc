/*
 * UniformDistributionModel.java
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

import dr.math.*;
import dr.xml.*;
import dr.inference.model.*;

/**
 * A class that acts as a model for normally distributed data.
 * @author Alexei Drummond
 * $Id$
 */

public class UniformDistributionModel extends AbstractModel implements ParametricDistributionModel {


	public static final String UNIFORM_DISTRIBUTION_MODEL = "uniformDistributionModel";
	public static final String LOWER = "lower";
	public static final String UPPER = "upper";

	/**
	 * Constructor.
	 */
	public UniformDistributionModel(Parameter lowerParameter, Parameter upperParameter) {

		super(UNIFORM_DISTRIBUTION_MODEL);

		this.lowerParameter = lowerParameter;
		addParameter(lowerParameter);
		lowerParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
		this.upperParameter = upperParameter;
		addParameter(upperParameter);
		upperParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
	}

	public double getLower() {
		return lowerParameter.getParameterValue(0);
	}

	public double getUpper() {
		return upperParameter.getParameterValue(0);
	}

	// *****************************************************************
	// Interface Distribution
	// *****************************************************************

	public double pdf(double x) { return UniformDistribution.pdf(x, getLower(), getUpper()); }
	public double logPdf(double x) { return UniformDistribution.logPdf(x, getLower(), getUpper()); }
	public double cdf(double x) { return UniformDistribution.cdf(x, getLower(), getUpper()); }
	public double quantile(double y) { return UniformDistribution.quantile(y, getLower(), getUpper()); }
	public double mean() { return UniformDistribution.mean(getLower(), getUpper()); }
	public double variance() { return UniformDistribution.variance(getLower(), getUpper()); }

	public final UnivariateFunction getProbabilityDensityFunction() { return pdfFunction; }

	private UnivariateFunction pdfFunction = new UnivariateFunction() {
		public final double evaluate(double x) { return 1.0; }
		public final double getLowerBound() { return getLower(); }
		public final double getUpperBound() { return getUpper(); }
	};

	// *****************************************************************
	// Interface Model
	// *****************************************************************

	public void handleModelChangedEvent(Model model, Object object, int index) {
		// no intermediates need to be recalculated...
	}

	public void handleParameterChangedEvent(Parameter parameter, int index) {
		// no intermediates need to be recalculated...
	}

	protected void storeState() {} // no additional state needs storing
	protected void restoreState() {} // no additional state needs restoring
	protected void acceptState() {} // no additional state needs accepting
	protected void adoptState(Model source) {} // no additional state needs adopting

	/**
	 * Reads a normal distribution model from a DOM Document element.
	 */
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return UNIFORM_DISTRIBUTION_MODEL; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			Parameter upperParam;
			Parameter lowerParam;

			XMLObject cxo = (XMLObject)xo.getChild(LOWER);
			if (cxo.getChild(0) instanceof Parameter) {
				lowerParam = (Parameter)cxo.getChild(Parameter.class);
			} else {
				lowerParam = new Parameter.Default(cxo.getDoubleChild(0));
			}

			cxo = (XMLObject)xo.getChild(UPPER);
			if (cxo.getChild(0) instanceof Parameter) {
				upperParam = (Parameter)cxo.getChild(Parameter.class);
			} else {
				upperParam = new Parameter.Default(cxo.getDoubleChild(0));
			}

			return new UniformDistributionModel(lowerParam, upperParam);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new ElementRule(LOWER,
				new XMLSyntaxRule[] {
					new XORRule(
						new ElementRule(Parameter.class),
						new ElementRule(Double.class)
					)}
				),
			new ElementRule(UPPER,
				new XMLSyntaxRule[] {
					new XORRule(
						new ElementRule(Parameter.class),
						new ElementRule(Double.class)
					)}
				)
		};

		public String getParserDescription() {
			return "Describes a normal distribution with a given mean and standard deviation " +
				"that can be used in a distributionLikelihood element";
		}

		public Class getReturnType() { return UniformDistributionModel.class; }
	};

	// **************************************************************
    // Private instance variables
    // **************************************************************

	private Parameter lowerParameter;
	private Parameter upperParameter;

}
