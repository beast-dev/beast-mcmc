/*
 * RandomWalkOperator.java
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

package dr.inference.operators;

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inferencexml.operators.RandomWalkOperatorParser;
import dr.math.MathUtils;
import dr.util.Transform;

import java.util.ArrayList;
import java.util.List;

/**
 * A generic random walk operator for use with a multi-dimensional parameters.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: RandomWalkOperator.java,v 1.16 2005/06/14 10:40:34 rambaut Exp $
 */
public class RandomWalkOperator extends AbstractAdaptableOperator {

    public enum BoundaryCondition {
        rejecting,
        reflecting,
        absorbing,
        log,
        logit
    }

    public RandomWalkOperator(Parameter parameter, double windowSize, BoundaryCondition bc, double weight, AdaptationMode mode) {
        this(parameter, null, windowSize, bc, weight, mode);
    }

    public RandomWalkOperator(Parameter parameter, Parameter updateIndex, double windowSize, BoundaryCondition boundaryCondition,
                              double weight, AdaptationMode mode) {

        this(parameter, windowSize, boundaryCondition, weight, mode, makeUpdateMap(updateIndex));
    }

    public RandomWalkOperator(Parameter parameter, double windowSize, BoundaryCondition boundaryCondition,
                              double weight, AdaptationMode mode, List<Integer> updateMap) {
        super(mode);

        setWeight(weight);
        this.parameter = parameter;
        this.windowSize = windowSize;
        this.boundaryCondition = boundaryCondition;
        this.updateMap = updateMap;
        if (updateMap != null) {
            updateMapSize = updateMap.size();
        }
    }

    private static ArrayList<Integer> makeUpdateMap(Parameter updateIndex) {
        ArrayList<Integer> updateMap = null;
        if (updateIndex != null) {
            updateMap = new ArrayList<Integer>();
            for (int i = 0; i < updateIndex.getDimension(); i++) {
                if (updateIndex.getParameterValue(i) == 1.0)
                    updateMap.add(i);
            }
        }
        return updateMap;
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

    public final BoundaryCondition getBoundaryCondition() {
        return boundaryCondition;
    }

    public final List<Integer> getUpdateMap() {
        return updateMap;
    }

    /**
     * change the parameter and return the hastings ratio.
     */
    public double doOperation() {

        // a random dimension to perturb
        if (parameter.getDimension() <= 0) {
            throw new RuntimeException("Illegal Dimension");
        }

        int dim;
        if (updateMap == null) {
            dim = MathUtils.nextInt(parameter.getDimension());
        } else {
            dim = updateMap.get(MathUtils.nextInt(updateMapSize));
        }

        // a random point around old value within windowSize * 2
        double draw = (2.0 * MathUtils.nextDouble() - 1.0) * windowSize;

        double oldValue = parameter.getParameterValue(dim);

        final Bounds<Double> bounds = parameter.getBounds();
        final double lower = bounds.getLowerLimit(dim);
        final double upper = bounds.getUpperLimit(dim);

        if (boundaryCondition == BoundaryCondition.logit) {
            // scale oldValue to [0,1]
            double x1 = (oldValue - lower) / (upper - lower);
            // logit transform it, add the draw, inverse transform it
            double x2 = Transform.LOGIT.inverse(Transform.LOGIT.transform(x1) + draw);

            // parameter takes new value scaled back into interval [lower, upper]
            parameter.setParameterValue(dim, (x2 * (upper - lower)) + lower);
            
            // HR is the ratio of Jacobians for the before and after values in interval [0,1]
            return Transform.LOGIT.getLogJacobian(x1) - Transform.LOGIT.getLogJacobian(x2);

        } else if (boundaryCondition == BoundaryCondition.log) {
            // offset oldValue to [0,+Inf]
            double x1 = oldValue - lower;
            // logit transform it, add the draw, inverse transform it
            double x2 = Transform.LOG.inverse(Transform.LOG.transform(x1) + draw);

            // parameter takes new value tranlated back into interval [lower, +Inf]
            parameter.setParameterValue(dim, x2 + lower);

            // HR is the ratio of Jacobians for the before and after values
            return Transform.LOG.getLogJacobian(x1) - Transform.LOG.getLogJacobian(x2);

        } else {
            double newValue = oldValue + draw;
            if (boundaryCondition == BoundaryCondition.reflecting) {
                newValue = reflectValue(newValue, lower, upper);
            } else if (boundaryCondition == BoundaryCondition.absorbing && (newValue < lower || newValue > upper)) {
                return 0.0;
            } else if (boundaryCondition == BoundaryCondition.rejecting && (newValue < lower || newValue > upper)) {
                return Double.NEGATIVE_INFINITY;
            }

            parameter.setParameterValue(dim, newValue);

            if (parameter.check()) {
                return 0.0;
            } else {
                return Double.NEGATIVE_INFINITY;
            }

        }

}

    public static double reflectValue(double value, double lower, double upper) {

        double newValue = value;

        if (upper == lower) {
            newValue = upper;
        } else if (value < lower) {
            if (Double.isInfinite(upper)) {
                // we are only going to reflect once as the upper bound is at infinity...
                newValue = lower + (lower - value);
            } else {
//                double remainder = lower - value;
//
//                double widths = Math.floor(remainder / (upper - lower));
//                remainder -= (upper - lower) * widths;

                final double ratio = (lower - value) / (upper - lower);
                final double widths = Math.floor(ratio);
                final double remainder = (ratio - widths) * (upper - lower);

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

//                double remainder = value - upper;
//
//                double widths = Math.floor(remainder / (upper - lower));
//                remainder -= (upper - lower) * widths;

                final double ratio = (value - upper) / (upper - lower);
                final double widths = Math.floor(ratio);
                final double remainder = (ratio - widths) * (upper - lower);

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
        return RandomWalkOperatorParser.RANDOM_WALK_OPERATOR + "(" + parameter.getParameterName() + ")";
    }

    @Override
    protected double getAdaptableParameterValue() {
        return Math.log(windowSize);
    }

    @Override
    protected void setAdaptableParameterValue(double value) {
        windowSize = Math.exp(value);
    }

    @Override
    public double getRawParameter() {
        return windowSize;
    }

    @Override
    public String getAdaptableParameterName() {
        return "windowSize";
    }

    public String toString() {
        return RandomWalkOperatorParser.RANDOM_WALK_OPERATOR + "(" + parameter.getParameterName() + ", " + windowSize + ", " + getWeight() + ")";
    }

    //PRIVATE STUFF

    protected Parameter parameter = null;
    private double windowSize = 0.01;
    private List<Integer> updateMap;
    private int updateMapSize;
    private final BoundaryCondition boundaryCondition;
}
