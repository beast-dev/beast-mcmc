/*
 * LogRandomWalkOperator.java
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
import dr.inferencexml.operators.LogRandomWalkOperatorParser;
import dr.math.MathUtils;

/**
 * A generic scale operator for use with a multi-dimensional parameters.
 * Either scale all dimentions at once or scale one dimention at a time.
 *
 * @author Joseph Heled
 * @version $Id$
 */
@Deprecated // Doing a random walk on a log scaled parameter is now supported in RandomWalkOperator
public class LogRandomWalkOperator extends AbstractCoercableOperator {

    private double size;
    private Parameter parameter = null;
    private final boolean scaleAll;
    private final boolean scaleAllInd;


    public LogRandomWalkOperator(Parameter parameter, double size,
                                 CoercionMode mode, double weight, boolean scaleAll, boolean scaleAllInd) {

        super(mode);
        setWeight(weight);

        this.parameter = parameter;
        assert parameter.getDimension() == 1;
        this.size = size;

        this.scaleAll = scaleAll;
        this.scaleAllInd = scaleAllInd;
    }


    private double scaleOne(int j) {
        final double w = size * (MathUtils.nextDouble() - 0.5);
        final double newValue = parameter.getParameterValue(j) * Math.exp(w);
        parameter.setParameterValue(j, newValue);

        return w;
    }

    /**
     * change the parameter and return the hastings ratio.
     */
    public final double doOperation() {
        final int dim = parameter.getDimension();
        final Bounds<Double> bounds = parameter.getBounds();
        double hastingsRatio = 0;

        int checkStart = 0, checkEnd = dim;

        // Must first set all parameters first and check for boundaries later for the operator to work
        // correctly with dependent parameters such as tree node heights.

        if( scaleAllInd ) {
            for(int i = 0; i < dim; i++) {
                hastingsRatio += scaleOne(i);
            }
        } else if( scaleAll ) {
            final double w = size * (MathUtils.nextDouble() - 0.5);
            final double f = Math.exp(w);
            for(int i = 0; i < dim; i++) {
                parameter.setParameterValue(i, parameter.getParameterValue(i) * f);
            }
            hastingsRatio += dim * w;
        } else {
            int j = MathUtils.nextInt(dim);
            hastingsRatio += scaleOne(j);
            checkStart = j;
            checkEnd = j + 1;
        }

        for(int i = checkStart; i < checkEnd; i++) {
            final double value = parameter.getParameterValue(i);
            if( value < bounds.getLowerLimit(i) || value > bounds.getUpperLimit(i) ) {
//                throw new OperatorFailedException("proposed value outside boundaries");
                return Double.NEGATIVE_INFINITY;
            }
        }
        return hastingsRatio;
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return LogRandomWalkOperatorParser.LOGRANDOMWALK_OPERATOR + (scaleAllInd ? "-all" : "") +
                (scaleAllInd ? "-independently" : "") +
                "(" + parameter.getParameterName() + ")";
    }

    public double getCoercableParameter() {
        return size;
    }

    public void setCoercableParameter(double value) {
        if( getMode() != CoercionMode.COERCION_OFF ) {
            size = value;
        }
    }

    public double getRawParameter() {
        return size;
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    public final String getPerformanceSuggestion() {

        double prob = Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();
        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
        double sf = size;
        if( prob < getMinimumGoodAcceptanceLevel() ) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if( prob > getMaximumGoodAcceptanceLevel() ) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else {
            return "";
        }
    }

    public String toString() {
        return getOperatorName();
    }
}