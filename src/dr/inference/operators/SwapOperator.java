/*
 * SwapOperator.java
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
import dr.inferencexml.operators.SwapOperatorParser;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * A generic operator swapping a number of pairs in a multi-dimensional parameter.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class SwapOperator extends SimpleMCMCOperator {
    private int size = 1;

    public SwapOperator(Parameter parameter, int size) {
        this.parameter = parameter;
        this.size = size;
        if (parameter.getDimension() < 2 * size) {
            throw new IllegalArgumentException();
        }

        int dimension = parameter.getDimension();
        List<Integer> list = new ArrayList<Integer>();
        for (int i = 0; i < dimension; i++) {
            list.add(i);
        }
        masterList = Collections.unmodifiableList(list);
    }

    /**
     * @return the parameter this operator acts on.
     */
    public Parameter getParameter() {
        return parameter;
    }

    /**
     * swap the values in two random parameter slots.
     */
    public final double doOperation() {

        List<Integer> allIndices = new ArrayList<Integer>(masterList);
        int left, right;

        for (int i = 0; i < size; i++) {
            left = allIndices.remove(MathUtils.nextInt(allIndices.size()));
            right = allIndices.remove(MathUtils.nextInt(allIndices.size()));
            double value1 = parameter.getParameterValue(left);
            double value2 = parameter.getParameterValue(right);
            parameter.setParameterValue(left, value2);
            parameter.setParameterValue(right, value1);
        }

        return 0.0;
    }

    public String getOperatorName() {
        return SwapOperatorParser.SWAP_OPERATOR + "(" + parameter.getParameterName() + ")";
    }

    public String getPerformanceSuggestion() {
        return "No suggestions";
    }

    //PRIVATE STUFF

    private Parameter parameter = null;
    private List<Integer> masterList = null;
}
