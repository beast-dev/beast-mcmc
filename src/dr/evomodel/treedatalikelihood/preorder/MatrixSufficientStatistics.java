/*
 * MatrixSufficientStatistics.java
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

package dr.evomodel.treedatalikelihood.preorder;

import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.DenseMatrix64F;

/**
 * @author Marc A. Suchard
 * @author Paul Bastide
 */

public class MatrixSufficientStatistics extends NormalSufficientStatistics {

    private final DenseMatrix64F actualization;

    MatrixSufficientStatistics(double[] displacement,
                               double[] precision,
                               double[] actualization,
                               int index,
                               int dim,
                               DenseMatrix64F Pd,
                               PrecisionType precisionType) {

        super(displacement, precision, index, dim, Pd, precisionType);

        int actualizationOffset = (dim * dim) * index;
        this.actualization = MissingOps.wrap(actualization, actualizationOffset, dim, dim);

    }

    public double getDisplacement(int row) {
        return getMean(row);
    }

    public double getActualization(int row, int col) {
        return actualization.unsafe_get(row, col);
    }

    public DenseMatrix64F getRawDisplacement() {
        return getRawMean();
    }

    public DenseMatrix64F getRawActualization() {
        return actualization;
    }

}
