/*
 * MaskedGradient.java
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

import dr.inference.model.*;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */

public class MaskedGradient implements GradientWrtParameterProvider, VariableListener {

    private final int dimension;
    private final GradientWrtParameterProvider gradient;
    private final Parameter parameter;

    private final int[] map;

    public MaskedGradient(GradientWrtParameterProvider gradient, Parameter mask) {
        this.gradient = gradient;
        mask.addVariableListener(this);

        final Parameter originalParameter = gradient.getParameter();
        this.parameter = new MaskedParameter(originalParameter, mask, true);

        this.map = new int[originalParameter.getDimension()];
        int[] inverseMap = new int[originalParameter.getDimension()];

        this.dimension = MaskedParameter.updateMask(mask, map, inverseMap, 1);
    }

    @Override
    public Likelihood getLikelihood() {
        return gradient.getLikelihood();
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return dimension;
    }

    @Override
    public double[] getGradientLogDensity() {

        double[] originalGradient = gradient.getGradientLogDensity();
        double[] result = new double[dimension];

        for (int i = 0; i < dimension; ++i) {
            result[i] = originalGradient[map[i]];
        }

        return result;
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        throw new RuntimeException("Changing mask is not implemented");
    }
}
