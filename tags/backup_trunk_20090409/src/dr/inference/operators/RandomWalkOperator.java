/*
 * RandomWalkOperator.java
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

import java.util.ArrayList;
import java.util.List;

/**
 * A generic random walk operator for use with a multi-dimensional parameters.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: RandomWalkOperator.java,v 1.16 2005/06/14 10:40:34 rambaut Exp $
 */
public class RandomWalkOperator extends AbstractCoercableOperator {

    public static final String WINDOW_SIZE = "windowSize";
    public static final String UPDATE_INDEX = "updateIndex";
    public static final String UPPER = "upper";
    public static final String LOWER = "lower";

    public static final String BOUNDARY_CONDITION = "boundaryCondition";

    enum BoundaryCondition {
        reflecting,
        absorbing
    }

    public RandomWalkOperator(Parameter parameter, double windowSize, BoundaryCondition bc, double weight, CoercionMode mode) {

        super(mode);
        this.parameter = parameter;
        this.windowSize = windowSize;
        this.condition = bc;

        setWeight(weight);
    }

    public RandomWalkOperator(Parameter parameter, Parameter updateIndex, double windowSize, BoundaryCondition bc, double weight, CoercionMode mode) {
        super(mode);
        this.parameter = parameter;
        this.windowSize = windowSize;
        this.condition = bc;

        setWeight(weight);
        updateMap = new ArrayList<Integer>();
        for (int i = 0; i < updateIndex.getDimension(); i++) {
            if (updateIndex.getParameterValue(i) == 1.0)
                updateMap.add(i);
        }
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
        if (updateMap == null)
            index = MathUtils.nextInt(parameter.getDimension());
        else
            index = updateMap.get(MathUtils.nextInt(updateMap.size()));

        // a random point around old value within windowSize * 2
        double newValue = parameter.getParameterValue(index) + ((2.0 * MathUtils.nextDouble() - 1.0) * windowSize);

        double lower = parameter.getBounds().getLowerLimit(index);
        double upper = parameter.getBounds().getUpperLimit(index);

        while (newValue < lower || newValue > upper) {
            if (newValue < lower) {
                if (condition == BoundaryCondition.reflecting) {
                    newValue = lower + (lower - newValue);
                } else {
                    throw new OperatorFailedException("proposed value outside boundaries");
                }

            }
            if (newValue > upper) {
                if (condition == BoundaryCondition.reflecting) {
                    newValue = upper - (newValue - upper);
                } else {
                    throw new OperatorFailedException("proposed value outside boundaries");
                }

            }
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
            return "randomWalkOperator";
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            CoercionMode mode = CoercionMode.parseMode(xo);

            double weight = xo.getDoubleAttribute(WEIGHT);
            double windowSize = xo.getDoubleAttribute(WINDOW_SIZE);
            Parameter parameter = (Parameter) xo.getChild(Parameter.class);

            BoundaryCondition condition = BoundaryCondition.valueOf(
                    xo.getAttribute(BOUNDARY_CONDITION, BoundaryCondition.reflecting.name()));

            if (xo.hasChildNamed(UPDATE_INDEX)) {
                XMLObject cxo = (XMLObject) xo.getChild(UPDATE_INDEX);
                Parameter updateIndex = (Parameter) cxo.getChild(Parameter.class);
                if (updateIndex.getDimension() != parameter.getDimension())
                    throw new RuntimeException("Parameter to update and missing indices must have the same dimension");
                return new RandomWalkOperator(parameter, updateIndex, windowSize, condition,
                        weight, mode);
            }

            return new RandomWalkOperator(parameter, windowSize, condition, weight, mode);
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
                new StringAttributeRule(BOUNDARY_CONDITION, null, BoundaryCondition.values(), true),
                new ElementRule(Parameter.class)
        };

    };

    public String toString() {
        return "randomWalkOperator(" + parameter.getParameterName() + ", " + windowSize + ", " + getWeight() + ")";
    }

    //PRIVATE STUFF

    private Parameter parameter = null;
    private double windowSize = 0.01;
    private List<Integer> updateMap = null;
    private BoundaryCondition condition;
}