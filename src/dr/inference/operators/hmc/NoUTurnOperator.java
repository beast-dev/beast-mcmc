/*
 * NoUTurnOperator.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.operators.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.GeneralOperator;

/**
 * @author Marc A. Suchard
 */

public class NoUTurnOperator extends HamiltonianMonteCarloOperator implements GeneralOperator {

    public NoUTurnOperator(CoercionMode mode, double weight, GradientWrtParameterProvider gradientProvider,
                           Parameter parameter, Parameter mask, double stepSize, int nSteps, double drawVariance) {
        super(mode, weight, gradientProvider, parameter, null, mask, stepSize, nSteps, drawVariance);
    }

    @Override
    public double doOperation(Likelihood likelihood) {
        return leafFrog(); // TODO Implement NUTS algorithm
    }

    @Override
    public double doOperation() {
        throw new IllegalStateException("Should not be called.");
    }

    private void buildTree() {

    }

    private double findReasonableEpsilon() {
        return 0.0;
    }

    private double computeStopCriterion() {
        return 0.0;
    }
}
