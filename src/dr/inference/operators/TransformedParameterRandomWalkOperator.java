/*
 * TransformedParameterRandomWalkOperator.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.inference.model.TransformedParameter;

/**
 * @author Paul Bastide
 */

public class TransformedParameterRandomWalkOperator extends RandomWalkOperator {

    public TransformedParameterRandomWalkOperator(TransformedParameter parameter, double windowSize, BoundaryCondition bc, double weight, CoercionMode mode) {
        super(parameter, windowSize, bc, weight, mode);
    }

    public TransformedParameterRandomWalkOperator(RandomWalkOperator randomWalkOperator) {
        super((TransformedParameter) randomWalkOperator.getParameter(),
                randomWalkOperator.getWindowSize(),
                randomWalkOperator.getBoundaryCondition(),
                randomWalkOperator.getWeight(),
                randomWalkOperator.getMode());
    }

    @Override
    public double doOperation() {
        // Store old states
        double[] oldValues = ((TransformedParameter) parameter).getParameterUntransformedValues();
        // Do operation
        double ratio = super.doOperation();
        // New states
        double[] newValues = ((TransformedParameter) parameter).getParameterUntransformedValues();
        // Compute Jacobians
        ratio += ((TransformedParameter) parameter).diffLogJacobian(oldValues, newValues);
        return ratio;
    }

//    public class TransformedMultivariateParameterRandomWalkOperator extends TransformedParameterRandomWalkOperator {
//        public TransformedMultivariateParameterRandomWalkOperator(TransformedMultivariateParameter parameter, double windowSize, BoundaryCondition bc, double weight, CoercionMode mode) {
//            super(parameter, windowSize, bc, weight, mode);
//        }
//
//        @Override
//        public final double doOperation() {
//            double[] oldValues = parameter.getParameterValues();
//            double ratio = super.doOperation();
//            double[] newValues = parameter.getParameterValues();
//            ratio += (TransformedMultivariateParameter)parameter.diffLogJacobian(oldValues, newValues);
//            return ratio;
//        }
//    }
}
