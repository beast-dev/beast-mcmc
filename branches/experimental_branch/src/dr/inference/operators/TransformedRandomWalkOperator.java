/*
 * RandomWalkOperator.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import java.util.ArrayList;
import java.util.List;

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inferencexml.operators.RandomWalkOperatorParser;
import dr.math.MathUtils;
import dr.util.Transform;

/**
 * A random walk operator that takes a (list of) parameter(s) and corresponding transformations.
 *
 * @author Guy Baele
 */
public class TransformedRandomWalkOperator extends AbstractCoercableOperator {

    private Parameter parameter = null;
    private final Transform[] transformations;
    private double windowSize = 0.01;
    private List<Integer> updateMap = null;
    private int updateMapSize;
    private final BoundaryCondition condition;

    private final Double lowerOperatorBound;
    private final Double upperOperatorBound;

    public enum BoundaryCondition {
        reflecting,
        absorbing
    }

    public TransformedRandomWalkOperator(Parameter parameter, Transform[] transformations, double windowSize, BoundaryCondition bc, double weight, CoercionMode mode) {
        this(parameter, transformations, null, windowSize, bc, weight, mode);
    }

    public TransformedRandomWalkOperator(Parameter parameter, Transform[] transformations, Parameter updateIndex, double windowSize, BoundaryCondition bc,
            double weight, CoercionMode mode) {
        this(parameter, transformations, updateIndex, windowSize, bc, weight, mode, null, null);
    }

    public TransformedRandomWalkOperator(Parameter parameter, Transform[] transformations, Parameter updateIndex, double windowSize, BoundaryCondition bc,
            double weight, CoercionMode mode, Double lowerOperatorBound, Double upperOperatorBound) {
        super(mode);
        this.parameter = parameter;
        this.transformations = transformations;
        this.windowSize = windowSize;
        this.condition = bc;

        setWeight(weight);
        if (updateIndex != null) {
            updateMap = new ArrayList<Integer>();
            for (int i = 0; i < updateIndex.getDimension(); i++) {
                if (updateIndex.getParameterValue(i) == 1.0)
                    updateMap.add(i);
            }
            updateMapSize=updateMap.size();
        }

        this.lowerOperatorBound = lowerOperatorBound;
        this.upperOperatorBound = upperOperatorBound;
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

        //store MH-ratio in logq
        double logJacobian = 0.0;

        // a random dimension to perturb
        int index;
        if (updateMap == null) {
            index = MathUtils.nextInt(parameter.getDimension());
        } else {
            index = updateMap.get(MathUtils.nextInt(updateMapSize));
        }

        //System.err.println("index: " + index);

        // a random point around old value within windowSize * 2
        double draw = (2.0 * MathUtils.nextDouble() - 1.0) * windowSize;

        //System.err.println("draw: " + draw);

        //transform parameter values first
        double[] x = parameter.getParameterValues();
        int dim = parameter.getDimension();

        //System.err.println("parameter " + index + ": " + x[index]);

        double[] transformedX = new double[dim];
        for (int i = 0; i < dim; i++) {
            transformedX[i] = transformations[i].transform(x[i]);
        }

        //System.err.println("transformed parameter " + index + ": " + transformedX[index]);

        //double newValue = parameter.getParameterValue(index) + draw;
        double newValue = transformedX[index] + draw;

        //System.err.println("new value: " + newValue);

        final Bounds<Double> bounds = parameter.getBounds();
        final double lower = (lowerOperatorBound == null ? bounds.getLowerLimit(index) : Math.max(bounds.getLowerLimit(index), lowerOperatorBound));
        final double upper = (upperOperatorBound == null ? bounds.getUpperLimit(index) : Math.min(bounds.getUpperLimit(index), upperOperatorBound));

        if (condition == BoundaryCondition.reflecting) {
            newValue = reflectValue(newValue, lower, upper);
        } else if (newValue < lower || newValue > upper) {
            throw new OperatorFailedException("proposed value outside boundaries");
        }

        //parameter.setParameterValue(index, newValue);
        parameter.setParameterValue(index, transformations[index].inverse(newValue));

        //System.err.println("set parameter to: " + parameter.getValue(index));

        //this should be correct
        //logJacobian += transformations[index].getLogJacobian(parameter.getParameterValue(index)) - transformations[index].getLogJacobian(x[index]);

        logJacobian += transformations[index].getLogJacobian(x[index]) - transformations[index].getLogJacobian(parameter.getParameterValue(index));  

        //return 0.0;
        return logJacobian;
    }

    public double reflectValue(double value, double lower, double upper) {

        double newValue = value;

        if (value < lower) {
            if (Double.isInfinite(upper)) {
                // we are only going to reflect once as the upper bound is at infinity...
                newValue = lower + (lower - value);
            } else {
                double remainder = lower - value;

                double widths = Math.floor(remainder / (upper - lower));
                remainder -= (upper - lower) * widths;

                // even reflections
                if (widths % 2 == 0) {
                    newValue = lower + remainder;
                    // odd reflections
                } else {
                    newValue = upper - remainder;
                }
            }
        } else if (value > upper) {
            if (Double.isInfinite(lower)) {
                // we are only going to reflect once as the lower bound is at -infinity...
                newValue = upper - (newValue - upper);
            } else {

                double remainder = value - upper;

                double widths = Math.floor(remainder / (upper - lower));
                remainder -= (upper - lower) * widths;

                // even reflections
                if (widths % 2 == 0) {
                    newValue = upper - remainder;
                    // odd reflections
                } else {
                    newValue = lower + remainder;
                }
            }
        }

        return newValue;
    }

    public double reflectValueLoop(double value, double lower, double upper) {
        double newValue = value;

        while (newValue < lower || newValue > upper) {
            if (newValue < lower) {
                newValue = lower + (lower - newValue);
            }
            if (newValue > upper) {
                newValue = upper - (newValue - upper);

            }
        }

        return newValue;
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

    public String toString() {
        return RandomWalkOperatorParser.RANDOM_WALK_OPERATOR + "(" + parameter.getParameterName() + ", " + windowSize + ", " + getWeight() + ")";
    }

}
