/*
 * RandomWalkIntegerNodeHeightWeightedOperator.java
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

package dr.evomodel.operators;

import dr.inference.model.Parameter;
import dr.evomodelxml.operators.RandomWalkIntegerNodeHeightWeightedOperatorParser;
import dr.inference.operators.RandomWalkIntegerOperator;
import dr.math.MathUtils;

/**
 * @author Chieh-Hsi Wu
 *
 * The probability an internal node is picked to have its state changed depends on the node height.
 */
public class RandomWalkIntegerNodeHeightWeightedOperator extends RandomWalkIntegerOperator {

    private Parameter internalNodeHeights;

    public RandomWalkIntegerNodeHeightWeightedOperator(
            Parameter parameter, int windowSize, double weight, Parameter internalNodeHeights){
        super(parameter, windowSize, weight);
        this.internalNodeHeights = internalNodeHeights;
    }

    public double doOperation() {

        // a random dimension to perturb
        int index = MathUtils.randomChoicePDF(internalNodeHeights.getParameterValues());
      
        int newValue = calculateNewValue(index);
        parameter.setValue(index, newValue);

        return 0.0;
    }

    //MCMCOperator INTERFACE
    public String getOperatorName() {
        return "randomWalkIntegerNodeHeightWeighted(" + parameter.getId() + ")";
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
   
    public String toString() {
        return RandomWalkIntegerNodeHeightWeightedOperatorParser.RANDOM_WALK_INT_NODE_HEIGHT_WGT_OP +
                "(" + parameter.getId() + ", " + windowSize + ", " + getWeight() + ")";
    }
}
