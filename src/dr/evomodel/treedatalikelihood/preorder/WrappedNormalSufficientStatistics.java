/*
 * WrappedNormalSufficientStatistics.java
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

package dr.evomodel.treedatalikelihood.preorder;

import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.math.matrixAlgebra.ReadableMatrix;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedMatrix;
import dr.math.matrixAlgebra.WrappedVector;
import org.ejml.data.DenseMatrix64F;

/**
 * @author Marc A. Suchard
 */
public class WrappedNormalSufficientStatistics {

    private final WrappedVector mean;
    private final WrappedMatrix precision;
    private final WrappedMatrix variance;
    private final double precisionScalar;

    public WrappedNormalSufficientStatistics(WrappedVector mean,
                                             WrappedMatrix precision,
                                             WrappedMatrix variance) {
        this.mean = mean;
        this.precision = precision;
        this.variance = variance;
        this.precisionScalar = 1.0;
    }

    public WrappedNormalSufficientStatistics(double[] buffer,
                                             int index,
                                             int dim,
                                             DenseMatrix64F Pd,
                                             PrecisionType precisionType) {

        int partialOffset = (precisionType.getPartialsDimension(dim)) * index;
        this.mean = new WrappedVector.Raw(buffer, partialOffset, dim);
        if (precisionType == PrecisionType.SCALAR) {
            this.precision = new WrappedMatrix.Raw(Pd.getData(), 0, dim, dim);
            this.precisionScalar = buffer[partialOffset + dim];
            this.variance = null;
        } else {
            this.precisionScalar = 1.0;
            this.precision = new WrappedMatrix.Raw(buffer, partialOffset + dim, dim, dim);
            this.variance = new WrappedMatrix.Raw(buffer, partialOffset + dim + dim * dim, dim, dim);
        }
    }

    public WrappedVector getMean() { return mean; }

    public WrappedMatrix getPrecision() { return precision; }

    public WrappedMatrix getVariance() { return variance; }

    public double getMean(int row) {
        return mean.get(row);
    }

    public double getPrecision(int row, int col) {
        return precision.get(row, col);
    }

    public double getPrecisionScalar() { return precisionScalar; }

    public String toString() {
        return mean + " " + precision;
    }
}
