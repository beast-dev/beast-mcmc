/*
 * TransformedGradientWrtParameter.java
 *
 * Copyright (c) 2002-2023 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.util.Transform;
import dr.xml.Reportable;

/**
 * @author Marc Suchard
 * @author Xiang Ji
 */

public class TransformedGradientWrtParameter implements GradientWrtParameterProvider, Reportable {

    private final GradientWrtParameterProvider gradient;
    private final Transform transform;

    public TransformedGradientWrtParameter(GradientWrtParameterProvider gradient,
                                           Transform transform) {
        this.gradient = gradient;
        this.transform = transform;
    }
    @Override
    public Likelihood getLikelihood() {
        return null;
    }

    @Override
    public Parameter getParameter() {
        return null;
    }

    @Override
    public int getDimension() {
        return 0;
    }

    @Override
    public double[] getGradientLogDensity() {
        return new double[0];
    }

    @Override
    public String getReport() {
        return null;
    }
}
