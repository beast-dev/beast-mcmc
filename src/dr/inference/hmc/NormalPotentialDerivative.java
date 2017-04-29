/*
 * NormalPotentialDerivative.java
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

package dr.inference.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

/**
 * @author Max Tolkoff
 */
@Deprecated // TODO Should be implemented in NormalDistribution, etc.
public class NormalPotentialDerivative implements GradientWrtParameterProvider {
    double mean;
    double stdev;
    Parameter parameter;

    public NormalPotentialDerivative(double mean, double stdev, Parameter parameter){
        this.mean = mean;
        this.stdev = stdev;
        this.parameter = parameter;
    }

    @Override
    public Likelihood getLikelihood() {
        return null;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return parameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] derivative = new double[parameter.getDimension()];

        for (int i = 0; i < derivative.length; i++) {
            derivative[i] -= (parameter.getParameterValue(i) - mean) / (stdev * stdev);
            /* Sign change */
        }

        return derivative;
    }
}
