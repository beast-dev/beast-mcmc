/*
 * StochasticDifferentialEquationModel.java
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

package dr.evomodel.continuous;

/**
 * @author Marc Suchard
 *         <p/>
 *         Provides a numerical approximation to a SDE solution for a diffusion model.
 *         Approximation is via the Euler-Maruyama method.
 *         <p/>
 *         Projected uses:
 *         - Diffusion on a sphere
 *         - Incorporating geographic features
 */
public class StochasticDifferentialEquationModel extends MultivariateDiffusionModel {

    private static double maxTimeIncrement = 1E-3;
    private static int defaultNumberSteps = 100;

    protected double calculateLogDensity(double[] start, double[] stop, double time) {

        int numSteps = defaultNumberSteps;
        if (time / numSteps > maxTimeIncrement) {
            numSteps = (int) (time / maxTimeIncrement);
        }

        return 0;
    }

    private double[] getDriftVector(double[] X) {
        return null;
    }

    private double[] getVarianceMatrix(double[] X) {
        return null;
    }

}
