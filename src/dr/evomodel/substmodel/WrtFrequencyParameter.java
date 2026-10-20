/*
 * WrtFrequencyParameter.java
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

package dr.evomodel.substmodel;

import dr.evomodel.substmodel.DifferentialMassProvider.DifferentialWrapper.WrtParameter;
import dr.inference.model.Parameter;

import java.util.Arrays;

/**
 * Differential with respect to one (unconstrained) stationary frequency pi_k of a model whose
 * unnormalized generator is linear in pi.  Built by {@link DifferentiableSubstitutionModelUtil#frequencyFactory}.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class WrtFrequencyParameter implements WrtParameter {

    private final FrequencyModel frequencyModel;
    private final Parameter parameter;
    private final int dim;

    WrtFrequencyParameter(FrequencyModel frequencyModel, Parameter parameter, int dim) {
        this.frequencyModel = frequencyModel;
        this.parameter = parameter;
        this.dim = dim;
    }

    public FrequencyModel getFrequencyModel() {
        return frequencyModel;
    }

    public Parameter getParameter() {
        return parameter;
    }

    public int getDim() {
        return dim;
    }

    public double[] getDifferentialFrequencies() {
        double[] differential = new double[frequencyModel.getFrequencyCount()];
        differential[dim] = 1.0;
        return differential;
    }

    @Override
    public double getRate(int switchCase) {
        throw new UnsupportedOperationException("Frequency differentials do not use rate switch cases");
    }

    @Override
    public double getNormalizationDifferential() {
        return 1.0;
    }

    @Override
    public void setupDifferentialFrequencies(double[] differentialFrequencies, double[] frequencies) {
        Arrays.fill(differentialFrequencies, 0.0);
        differentialFrequencies[dim] = 1.0;
    }

    @Override
    public void setupDifferentialRates(double[] differentialRates, double[] relativeRates, double normalizingConstant) {
        for (int i = 0; i < relativeRates.length; i++) {
            differentialRates[i] = relativeRates[i] / normalizingConstant;
        }
    }
}
