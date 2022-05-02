/*
 * FireParameterOperator.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.inferencexml.operators.FireParameterOperatorParser;

/**
 * @author Marc Suchard
 * @author Gabriel Hassler
 */
public class FireParameterOperator extends SimpleMCMCOperator implements GibbsOperator {

    public FireParameterOperator(Parameter parameter, double[] values, double weight) {
        this.parameter = parameter;
        this.values = values;
        setWeight(weight);
    }

    public String getPerformanceSuggestion() {
        return null;
    }

    public String getOperatorName() {
        return FireParameterOperatorParser.FIRE_PARAMETER_OPERATOR;
    }

    public double doOperation() {
        if (values != null) {

            for (int i = 0; i < parameter.getDimension(); i++) {
                parameter.setParameterValueQuietly(i, values[i]);
            }
            parameter.fireParameterChangedEvent();

        } else {
            parameter.setParameterValue(0, parameter.getParameterValue(0));
        }
        return 0;
    }

    public int getStepCount() {
        return 1;
    }

    private Parameter parameter;
    private final double[] values;
}
