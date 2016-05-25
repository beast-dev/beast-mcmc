/*
 * UniformIntegerOperator.java
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

import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inferencexml.operators.UniformIntegerOperatorParser;
import dr.math.MathUtils;

/**
 * A generic uniform sampler/operator for use with a multi-dimensional parameter.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: UniformOperator.java,v 1.16 2005/06/14 10:40:34 rambaut Exp $
 */
public class UniformIntegerOperator extends SimpleMCMCOperator {
    private final int howMany;

    public UniformIntegerOperator(Parameter parameter, int lower, int upper, double weight, int howMany) {
        this(parameter, weight, howMany);
        this.lower = lower;
        this.upper = upper;
    }

    public UniformIntegerOperator(Parameter parameter, int lower, int upper, double weight) {
        this(parameter, lower, upper, weight, 1);
    }

    public UniformIntegerOperator(Variable parameter, double weight, int howMany) { // Bounds.Staircase
        this.parameter = parameter;
        this.howMany = howMany;
        setWeight(weight);
    }

    /**
     * @return the parameter this operator acts on.
     */
    public Parameter getParameter() {
        return (Parameter) parameter;
    }

    /**
     * @return the Variable this operator acts on.
     */
    public Variable getVariable() {
        return parameter;
    }

    /**
     * change the parameter and return the hastings ratio.
     */
    public final double doOperation() {

        for (int n = 0; n < howMany; ++n) {
            // do not worry about duplication, does not matter
            int index = MathUtils.nextInt(parameter.getSize());

            if (parameter instanceof Parameter) {
                int newValue = MathUtils.nextInt(upper - lower + 1) + lower; // from 0 to n-1, n must > 0,
                ((Parameter) parameter).setParameterValue(index, newValue);
            } else { // Variable<Integer>, Bounds.Staircase

                int upper = ((Variable<Integer>) parameter).getBounds().getUpperLimit(index);
                int lower = ((Variable<Integer>) parameter).getBounds().getLowerLimit(index);
                int newValue = MathUtils.nextInt(upper - lower + 1) + lower; // from 0 to n-1, n must > 0,
                ((Variable<Integer>) parameter).setValue(index, newValue);

            }

        }

        return 0.0;
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return "uniformInteger(" + parameter.getId() + ")";
    }

    public final void optimize(double targetProb) {

        throw new RuntimeException("This operator cannot be optimized!");
    }

    public boolean isOptimizing() {
        return false;
    }

    public void setOptimizing(boolean opt) {
        throw new RuntimeException("This operator cannot be optimized!");
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

    public String getPerformanceSuggestion() {
        if (Utils.getAcceptanceProbability(this) < getMinimumAcceptanceLevel()) {
            return "";
        } else if (Utils.getAcceptanceProbability(this) > getMaximumAcceptanceLevel()) {
            return "";
        } else {
            return "";
        }
    }

    public String toString() {
        return UniformIntegerOperatorParser.UNIFORM_INTEGER_OPERATOR + "(" + parameter.getId() + ")";
    }

    //PRIVATE STUFF

    private Variable parameter = null;
    private int upper;
    private int lower;
}
