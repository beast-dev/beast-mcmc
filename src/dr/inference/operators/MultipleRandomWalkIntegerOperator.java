/*
 * MultipleRandomWalkIntegerOperator.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.math.MathUtils;

/**
 * @author Chieh-Hsi Wu
 * Operator that performs independent random walks on k elements in a parameter
 */
public class MultipleRandomWalkIntegerOperator extends RandomWalkIntegerOperator{
    private int sampleSize = 1;
    public MultipleRandomWalkIntegerOperator(Parameter parameter, int windowSize, int sampleSize, double weight) {
        super(parameter,windowSize, weight);
        this.sampleSize = sampleSize;
    }

    public double doOperation() {

        int[] shuffledInd = MathUtils.shuffled(parameter.getSize()); // use getSize(), which = getDimension()
        int index;
        for(int i = 0; i < sampleSize; i++){
            index = shuffledInd[i];

            int newValue = calculateNewValue(index);
            parameter.setValue(index, newValue);
        }
        return 0.0;
    }

}
