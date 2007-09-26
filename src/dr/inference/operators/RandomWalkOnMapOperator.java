/*
 * RandomWalkOnMapOperator.java
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

package dr.inference.operators;

import dr.evomodel.continuous.MapDiffusionModel;
import dr.evomodel.continuous.TopographicalMap;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.*;


/**
 * @author Marc Suchard
 */
public class RandomWalkOnMapOperator extends SimpleMCMCOperator implements CoercableMCMCOperator {

	public static final String OPERATOR_NAME = "randomWalkOnMapOperator";
	public static final String WINDOW_SIZE = "windowSize";
	public static final String UPDATE_INDEX = "updateIndex";
	public static final String GRID_X_DIMENSION = "xGridDimension";
	public static final String GRID_Y_DIMENSION = "yGridDimension";

	public RandomWalkOnMapOperator(Parameter parameter,
//	                               int xDim, int yDim,
MapDiffusionModel mapModel,
double windowSize,
double weight, int mode) {
		this.parameter = parameter;
		this.map = mapModel.getMap();
		this.windowSize = windowSize;
		setWeight(weight);
		this.mode = mode;
		this.numberPoints = parameter.getDimension() / 2;

	}


	/**
	 * @return the parameter this operator acts on.
	 */
	public Parameter getParameter() {
		return parameter;
	}

	public final double getWindowSize() {
		return windowSize;
	}

	/**
	 * change the parameter and return the hastings ratio.
	 */
	public final double doOperation() throws OperatorFailedException {

		// a random dimension to perturb

		int index;
		index = MathUtils.nextInt(numberPoints) * 2;

		int deltaX = 0;
		int deltaY = 0;

		while (deltaX == 0 && deltaY == 0) {  // eight uniform choices
			deltaX = MathUtils.nextInt(3) - 1;
			deltaY = MathUtils.nextInt(3) - 1;
		}

		// a random point around old value within window [-1,1]  todo expand windowsize

		int newX = (int) parameter.getParameterValue(index) + deltaX;
		int newY = (int) parameter.getParameterValue(index + 1) + deltaY;

		if (!map.isValidPoint(newX, newY)) {
			throw new OperatorFailedException("proposed value outside boundaries");
		}

		parameter.setParameterValue(index, newX);
		parameter.setParameterValue(index + 1, newY);

		return 0.0;
	}

	//MCMCOperator INTERFACE
	public final String getOperatorName() {
		return parameter.getParameterName();
	}


	public double getCoercableParameter() {
		return Math.log(windowSize);
	}

	public void setCoercableParameter(double value) {
		windowSize = Math.exp(value);
	}

	public double getRawParameter() {
		return windowSize;
	}

	public int getMode() {
		return mode;
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

		double ws = OperatorUtils.optimizeWindowSize(windowSize, parameter.getParameterValue(0) * 2.0, prob, targetProb);

		if (prob < getMinimumGoodAcceptanceLevel()) {
			return "Try decreasing windowSize to about " + ws;
		} else if (prob > getMaximumGoodAcceptanceLevel()) {
			return "Try increasing windowSize to about " + ws;
		} else return "";
	}

	public static dr.xml.XMLObjectParser PARSER = new AbstractXMLObjectParser() {


		public String getParserName() {
			return OPERATOR_NAME;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			int mode = CoercableMCMCOperator.DEFAULT;
			if (xo.hasAttribute(AUTO_OPTIMIZE)) {
				if (xo.getBooleanAttribute(AUTO_OPTIMIZE)) {
					mode = CoercableMCMCOperator.COERCION_ON;
				} else {
					mode = CoercableMCMCOperator.COERCION_OFF;
				}
			}

			double weight = xo.getDoubleAttribute(WEIGHT);
			double windowSize = xo.getDoubleAttribute(WINDOW_SIZE);
			Parameter parameter = (Parameter) xo.getChild(Parameter.class);

			MapDiffusionModel mapModel = (MapDiffusionModel) xo.getChild(MapDiffusionModel.class);


			return new RandomWalkOnMapOperator(parameter, mapModel, windowSize, weight, mode);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element returns a random walk operator on a given map.";
		}

		public Class getReturnType() {
			return MCMCOperator.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				AttributeRule.newDoubleRule(WINDOW_SIZE),
				AttributeRule.newDoubleRule(WEIGHT),
				AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
				new ElementRule(MapDiffusionModel.class),
				new ElementRule(Parameter.class)
		};

	};

	public String toString() {
		return OPERATOR_NAME + "(" + parameter.getParameterName() + ")";
	}

	//PRIVATE STUFF

	private Parameter parameter = null;
	private double windowSize = 0.01;
	private int mode = CoercableMCMCOperator.DEFAULT;
	//	private int weight = 1;
	private TopographicalMap map;
	private int numberPoints;

}

