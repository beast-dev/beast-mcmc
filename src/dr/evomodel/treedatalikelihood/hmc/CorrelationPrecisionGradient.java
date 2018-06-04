/*
 * CorrelationPrecisionGradient.java
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

package dr.evomodel.treedatalikelihood.hmc;

import dr.inference.model.CompoundSymmetricMatrix;
import dr.inference.model.Likelihood;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;
import dr.math.matrixAlgebra.SymmetricMatrix;

/**
 * @author Paul Bastide
 * @author Marc A. Suchard
 */
public class CorrelationPrecisionGradient extends AbstractPrecisionGradient {

    public CorrelationPrecisionGradient(ConjugateWishartStatisticsProvider wishartStatistics,
                                        Likelihood likelihood,
                                        CompoundSymmetricMatrix parameter) {

        super(wishartStatistics, likelihood, parameter);
    }


    @Override
    public int getDimension() {
        return getDimensionCorrelation();
    }

    @Override
    double[] getGradientParameter(SymmetricMatrix weightedSumOfSquares,
                                  int numberTips,
                                  SymmetricMatrix correlationPrecision,
                                  double[] precisionDiagonal) {

        return getGradientCorrelation(weightedSumOfSquares, numberTips,
                correlationPrecision, precisionDiagonal);
    }

}
