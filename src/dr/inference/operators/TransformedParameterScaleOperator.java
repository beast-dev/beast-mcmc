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

import dr.inference.model.Parameter;
import dr.inference.model.TransformedParameter;
import dr.inference.model.Variable;
import dr.math.matrixAlgebra.Matrix;
import figtree.treeviewer.ScaleAxis;

/**
 * @author Paul Bastide
 */

public class TransformedParameterScaleOperator extends ScaleOperator {

    private static boolean DEBUG = false;

    public TransformedParameterScaleOperator(ScaleOperator scaleOperator) {
        super((TransformedParameter) scaleOperator.getVariable(),
                scaleOperator.getScaleFactor(),
                scaleOperator.getMode(),
                scaleOperator.getWeight());
    }

    public TransformedParameterScaleOperator(Variable variable, double scale) {

        super(variable, scale);
    }

    public TransformedParameterScaleOperator(Variable<Double> variable, double scale, AdaptationMode mode, double weight) {

        super(variable, scale, mode, weight);
    }

    public TransformedParameterScaleOperator(Variable<Double> variable, boolean scaleAll, int degreesOfFreedom, double scale,
                         AdaptationMode mode, Parameter indicator, double indicatorOnProb, boolean scaleAllInd) {

        super(variable, scaleAll, degreesOfFreedom,scale,mode, indicator, indicatorOnProb, scaleAllInd);

    }
    @Override
    public double doOperation() {
        // Store old states
        Parameter par = (Parameter)getVariable();
        double[] oldValues = par.getParameterValues();

        // Do operation
        double ratio = super.doOperation();
        // New states
        double[] newValues = par.getParameterValues();

        // Compute Jacobians
        ratio += ((TransformedParameter) par).diffLogJacobian(oldValues, newValues);

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
