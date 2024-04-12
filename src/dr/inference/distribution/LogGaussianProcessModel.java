/*
 * LogLinearModel.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.distribution;

import dr.evomodel.substmodel.LogAdditiveCtmcRateProvider;
import dr.inference.model.DesignMatrix;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */

// TODO this class is not really a LogLinearModel nor GeneralizedLinearModel; need to disassociate.
// TODO disassocation requires refactoring
public class LogGaussianProcessModel extends LogLinearModel implements LogAdditiveCtmcRateProvider {

    public LogGaussianProcessModel(Parameter dependentParameter) {
        super(dependentParameter);
    }

    @Override
    protected double calculateLogLikelihood() {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    public double[] getXBeta() {

        final int fieldDim = independentParam.get(0).getDimension();

        double[] rates = new double[fieldDim];

        for (int k = 0; k < numIndependentVariables; ++k) {
            Parameter field = independentParam.get(k);
            DesignMatrix X = designMatrix.get(k);

            for (int i = 0; i < fieldDim; ++i) {
                rates[i] += field.getParameterValue(i);
            }
        }

        for(int i = 0; i < rates.length; i++) {
            rates[i] = Math.exp(rates[i]);
        }
        return rates;
    }
}
