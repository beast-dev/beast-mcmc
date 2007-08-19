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

import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.*;


/**
 * @author Marc Suchard
 */
public class RandomWalkOnMapOperator extends SimpleMCMCOperator implements CoercableMCMCOperator {

    public static final String WINDOW_SIZE = "windowSize";
    public static final String UPDATE_INDEX = "updateIndex";
    public static final String GRID_X_DIMENSION = "xGridDimension";
    public static final String GRID_Y_DIMENSION = "yGridDimension";

    public RandomWalkOnMapOperator(Parameter parameter,
                                   int xDim, int yDim,
                                   double windowSize,
                                   double weight, int mode) {
        this.parameter = parameter;
        this.windowSize = windowSize;
        setWeight(weight);
        this.mode = mode;
        dim = new int[2];
        dim[0] = xDim;
        dim[1] = yDim;
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
        // todo (x,y) simultaneous update
        int index;
        index = MathUtils.nextInt(parameter.getDimension());

        // a random point around old value within windowSize * 2
        int newValue = (int) (parameter.getParameterValue(index) + ((2.0 * MathUtils.nextDouble() - 1.0) * windowSize));

        // check boundary
//			if (newValue < parameter.getBounds().getLowerLimit(index) || newValue > parameter.getBounds().getUpperLimit(index)) {
//				throw new OperatorFailedException("proposed value outside boundaries");
//			}

        // todo Reflect about graph edges
        if (newValue < 0 || newValue >= dim[0]) {    // todo fix boundaries
            throw new OperatorFailedException("proposed value outside boundaries");
        }


        parameter.setParameterValue(index, newValue);

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
            return "randomWalkOnMapOperator";
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

            int xDim = xo.getIntegerAttribute(GRID_X_DIMENSION);
            int yDim = xo.getIntegerAttribute(GRID_Y_DIMENSION);

            return new RandomWalkOnMapOperator(parameter, xDim, yDim, windowSize, weight, mode);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a random walk operator on a given parameter.";
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
                AttributeRule.newIntegerRule(GRID_X_DIMENSION),
                AttributeRule.newIntegerRule(GRID_Y_DIMENSION),
                new ElementRule(Parameter.class)
        };

    };

    public String toString() {
        return "randomWalkOnMapOperator(" + parameter.getParameterName() + ", " + windowSize + ", " + getWeight() + ")";
    }

    //PRIVATE STUFF

    private Parameter parameter = null;
    private double windowSize = 0.01;
    private int mode = CoercableMCMCOperator.DEFAULT;
    private int[] dim;
}

