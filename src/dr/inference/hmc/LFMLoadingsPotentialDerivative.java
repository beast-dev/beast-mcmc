/*
 * LFMLoadingsPotentialDerivative.java
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
import dr.inference.model.LatentFactorModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

/**
 * @author Max Tolkoff
 */
public class LFMLoadingsPotentialDerivative implements GradientWrtParameterProvider {
    LatentFactorModel lfm;

    public LFMLoadingsPotentialDerivative(LatentFactorModel lfm){
        this.lfm = lfm;
    }

    @Override
    public Likelihood getLikelihood() {
        return lfm;
    }

    @Override
    public Parameter getParameter() {
        return lfm.getLoadings();
    }

    @Override
    public int getDimension() {
        return lfm.getLoadings().getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] derivative = new double[lfm.getLoadings().getDimension()];
        Parameter missingIndicator = lfm.getMissingIndicator();
        int ntaxa = lfm.getFactors().getColumnDimension();
        int ntraits = lfm.getLoadings().getRowDimension();
        int nfac = lfm.getLoadings().getColumnDimension();
        double[] residual = lfm.getResidual();

        for (int i = 0; i < nfac; i++) {
            for (int j = 0; j < ntraits; j++) {
                for (int k = 0; k < ntaxa; k++) {
                    if(missingIndicator == null || missingIndicator.getParameterValue(k * ntraits + j) != 1){
                        derivative[i * ntraits + j] += lfm.getFactors().getParameterValue(i, k) * lfm.getColumnPrecision().getParameterValue(j, j) *
                                residual[k * ntraits + j];
                        /* Sign change */
                    }
                }
            }
        }

        return derivative;
    }
}
