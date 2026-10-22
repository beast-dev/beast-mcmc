/*
 * WrtRateParameter.java
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

/**
 * Differential with respect to one model parameter that sets a single exchangeability (relative rate) s_m of a
 * reversible generator Q_ij = s_ij pi_j, e.g. TN93 kappa1/kappa2 or a GTR rate.
 * Built by {@link DifferentiableSubstitutionModelUtil#rateFactory}.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class WrtRateParameter implements WrtParameter {

    private final Parameter parameter;
    private final int dim;
    private final int rateIndex;
    private final int rateCount;

    WrtRateParameter(Parameter parameter, int dim, int rateIndex, int rateCount) {
        this.parameter = parameter;
        this.dim = dim;
        this.rateIndex = rateIndex;
        this.rateCount = rateCount;
    }

    public Parameter getParameter() {
        return parameter;
    }

    public int getDim() {
        return dim;
    }

    public int getRateIndex() {
        return rateIndex;
    }

    // d relativeRates / d theta: the parameter value is the relative rate itself
    public double[] getDifferentialRates() {
        double[] differential = new double[rateCount];
        differential[rateIndex] = 1.0;
        return differential;
    }

    // The switch-case interface below belongs to the per-model enum differentials and is not used for rates

    @Override
    public double getRate(int switchCase) {
        throw new UnsupportedOperationException("Rate differentials do not use rate switch cases");
    }

    @Override
    public double getNormalizationDifferential() {
        throw new UnsupportedOperationException("Rate differentials are built by getRateDifferentialMatrix");
    }

    @Override
    public void setupDifferentialFrequencies(double[] differentialFrequencies, double[] frequencies) {
        throw new UnsupportedOperationException("Rate differentials are built by getRateDifferentialMatrix");
    }

    @Override
    public void setupDifferentialRates(double[] differentialRates, double[] relativeRates, double normalizingConstant) {
        throw new UnsupportedOperationException("Rate differentials are built by getRateDifferentialMatrix");
    }
}
