/*
 * BitFlipOperator.java
 *
 * Copyright (C) 2002-2007 Alexei Drummond and Andrew Rambaut
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

package dr.evomodel.operators;

import dr.evomodel.substmodel.SVSGeneralSubstitutionModel;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * A generic operator that flips bits.
 *
 * @author Marc Suchard
 * @version $Id$
 */
public class BitFlipInSubstitutionModelOperator extends SimpleMCMCOperator implements CoercableMCMCOperator {

	public static final String BIT_FLIP_OPERATOR = "bitFlipInSubstitutionModelOperator";
	public static final String SCALE_FACTOR = "scaleFactor";

	public BitFlipInSubstitutionModelOperator(SVSGeneralSubstitutionModel subModel, Parameter rateParameter, double weight, double scaleFactor, int mode) {
		this.subModel = subModel;
		this.indicatorParameter = subModel.getRateIndicators();
		this.rateParameter = rateParameter;
		this.mode = mode;
		this.scaleFactor = scaleFactor;
		setWeight(weight);

	}

	/**
	 * @return the indicatorParameter this operator acts on.
	 */
	public Parameter getIndicatorParameter() {
		return indicatorParameter;
	}

	/**
	 * Change the indicatorParameter and return the hastings ratio.
	 * Flip (Switch a 0 to 1 or 1 to 0) for a random bit in a bit vector.
	 * Return the hastings ratio which makes all subsets of vectors with the same number of 1 bits
	 * equiprobable.
	 */
	public final double doOperation() throws OperatorFailedException {
		final int dim = indicatorParameter.getDimension();
		double sum = 0.0;

		for (int i = 0; i < dim; i++) {
			sum += indicatorParameter.getParameterValue(i);
		}

		final int pos = MathUtils.nextInt(dim);

		final int value = (int) indicatorParameter.getParameterValue(pos);

		double rand = //0.5 *
				MathUtils.nextDouble();

		double logq;
		if (value == 0) {
			indicatorParameter.setParameterValue(pos, 1.0);

			logq = -Math.log((dim - sum) / (sum + 1));
//	        rand = 0.5 - rand;

		} else if (value == 1) {
			indicatorParameter.setParameterValue(pos, 0.0);

			logq = -Math.log(sum / (dim - sum + 1));
//	        rand = 0.5 + rand;
			rand *= -1;

		} else {
			throw new RuntimeException("expected 1 or 0");
		}

//	    final double scale = (scaleFactor + (rand * ((1.0/scaleFactor) - scaleFactor)));
//      logq -= Math.log(scale);
		final double scale = Math.exp((rand) * scaleFactor);
		logq += Math.log(scale);

		final double oldValue = rateParameter.getParameterValue(0);
		final double newValue = scale * oldValue;

		rateParameter.setParameterValue(0, newValue);

//        if (!subModel.myIsValid()) {
//		    System.err.println("invalid model");
//            throw new OperatorFailedException("Out of bounds");
//        } //                  else System.err.println("valid");

		// hastings ratio is designed to make move symmetric on sum of 1's
		return logq;
	}

	// Interface MCMCOperator
	public final String getOperatorName() {
		return "bitflip(" + indicatorParameter.getParameterName() + ")";
	}

	public double getCoercableParameter() {
//	     return Math.log(1.0 / scaleFactor - 1.0);
		return Math.log(scaleFactor);
	}

	public void setCoercableParameter(double value) {
//	     scaleFactor = 1.0 / (Math.exp(value) + 1.0);
		scaleFactor = Math.exp(value);
	}

	public double getRawParameter() {
		return scaleFactor;
	}

	public int getMode() {
		return mode;
	}

	public double getScaleFactor() {
		return scaleFactor;
	}

	public double getTargetAcceptanceProbability() {
		return 0.234;
	}

	public double getMinimumAcceptanceLevel() {
		return 0.1;
	}

	public double getMaximumAcceptanceLevel() {
		return 0.4;
	}

	public double getMinimumGoodAcceptanceLevel() {
		return 0.20;
	}

	public double getMaximumGoodAcceptanceLevel() {
		return 0.30;
	}

	public final String getPerformanceSuggestion() {

		double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
		double targetProb = getTargetAcceptanceProbability();
		dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
		double sf = OperatorUtils.optimizeWindowSize(scaleFactor, prob, targetProb);
		if (prob < getMinimumGoodAcceptanceLevel()) {
			return "Try setting scaleFactor to about " + formatter.format(sf);
		} else if (prob > getMaximumGoodAcceptanceLevel()) {
			return "Try setting scaleFactor to about " + formatter.format(sf);
		} else return "";
	}


	public String toString() {
		return getOperatorName();
	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return BIT_FLIP_OPERATOR;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			double weight = xo.getDoubleAttribute(WEIGHT);

			int mode = CoercableMCMCOperator.DEFAULT;

			if (xo.hasAttribute(AUTO_OPTIMIZE)) {
				mode = xo.getBooleanAttribute(AUTO_OPTIMIZE) ?
						CoercableMCMCOperator.COERCION_ON : CoercableMCMCOperator.COERCION_OFF;
			}

			final double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);

			if (scaleFactor <= 0.0 || scaleFactor >= 1.0) {
				throw new XMLParseException("scaleFactor must be between 0.0 and 1.0");
			}


			Parameter rateParameter = (Parameter) xo.getChild(Parameter.class);
			SVSGeneralSubstitutionModel subModel = (SVSGeneralSubstitutionModel) xo.getChild(SVSGeneralSubstitutionModel.class);

//	        if (xo.hasAttribute(MIN) && xo.hasAttribute(MAX)) {
//		        int min = xo.getIntegerAttribute(MIN);
//		        int max = xo.getIntegerAttribute(MAX);
//		        return new BitFlipOperator(indicatorParameter, weight, min, max);
//	        }

			return new BitFlipInSubstitutionModelOperator(subModel, rateParameter, weight, scaleFactor, mode);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element returns a bit-flip operator on a given indicatorParameter.";
		}

		public Class getReturnType() {
			return MCMCOperator.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				AttributeRule.newDoubleRule(WEIGHT),
				AttributeRule.newDoubleRule(SCALE_FACTOR),
				AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
//		        AttributeRule.newIntegerRule(MIN,true),
//		        AttributeRule.newIntegerRule(MAX,true),
				new ElementRule(Parameter.class),
				new ElementRule(SVSGeneralSubstitutionModel.class)
		};

	};
	// Private instance variables

	private Parameter indicatorParameter = null;
	private Parameter rateParameter = null;
	private SVSGeneralSubstitutionModel subModel;
	private double scaleFactor;
	private int mode;
}
