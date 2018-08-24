/*
 * BitMoveOperator.java
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
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * A generic operator that moves k 1 bits to k zero locations.
 *
 * @author Alexei Drummond
 * @version $Id$
 */
public class BitMoveOperator extends SimpleMCMCOperator {

    public BitMoveOperator(Parameter bitsParameter, Parameter valuesParameter, int numBitsToMove, double weight) {
        this.bitsParameter = bitsParameter;
        this.valuesParameter = valuesParameter;

        if (valuesParameter != null && bitsParameter.getDimension() != valuesParameter.getDimension()) {
            throw new IllegalArgumentException("bits parameter must be same length as values parameter");
        }

        this.numBitsToMove = numBitsToMove;
        setWeight(weight);
    }

    /**
     * Pick a random k ones in the vector and move them to a random k zero positions.
     */
    public final double doOperation() {

        final int dim = bitsParameter.getDimension();
        List<Integer> ones = new ArrayList<Integer>();
        List<Integer> zeros = new ArrayList<Integer>();

        for (int i = 0; i < dim; i++) {
            if (bitsParameter.getParameterValue(i) == 1.0) {
                ones.add(i);
            } else {
                zeros.add(i);
            }
        }

        if (ones.size() >= numBitsToMove && zeros.size() >= numBitsToMove) {

            for (int i = 0; i < numBitsToMove; i++) {

                int myOne = ones.remove(MathUtils.nextInt(ones.size()));
                int myZero = zeros.remove(MathUtils.nextInt(zeros.size()));

                bitsParameter.setParameterValue(myOne, 0.0);
                bitsParameter.setParameterValue(myZero, 1.0);

                if (valuesParameter != null) {
                    double value1 = valuesParameter.getParameterValue(myOne);
                    double value2 = valuesParameter.getParameterValue(myZero);
                    valuesParameter.setParameterValue(myOne, value2);
                    valuesParameter.setParameterValue(myZero, value1);
                }

            }
        } else {
            // might be better to return negative infinity and reject the move - there is likely
            // to be a bit flip operator on which will sort the problem out.
            //throw new RuntimeException("Not enough bits to move!");

            return Double.NEGATIVE_INFINITY;
        }

        return 0.0;
    }

    // Interface MCMCOperator
    public final String getOperatorName() {
        StringBuilder builder = new StringBuilder();
        builder.append("bitMove(");
        builder.append(bitsParameter.getParameterName());

        if (valuesParameter != null) {
            builder.append(", ").append(valuesParameter.getParameterName());
        }
        builder.append(", ").append(numBitsToMove).append(")");

        return builder.toString();
    }

    public final String getPerformanceSuggestion() {
        return "no performance suggestion";
    }

    public String toString() {
        return getOperatorName();
    }

    // Private instance variables

    private Parameter bitsParameter = null;
    private Parameter valuesParameter = null;
    private int numBitsToMove = 1;
}
