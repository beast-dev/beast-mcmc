/*
 * MaximumEigenvalue.java
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

package dr.math;

import dr.math.matrixAlgebra.*;

import static dr.math.matrixAlgebra.ReadableMatrix.Utils.product;
import static dr.math.matrixAlgebra.ReadableVector.Utils.innerProduct;
import static dr.math.matrixAlgebra.ReadableVector.Utils.norm;

/**
 * @author Marc A. Suchard
 */
public interface MaximumEigenvalue {

    double find(double[][] matrix);

    class PowerMethod implements  MaximumEigenvalue {

        private final int numIterations;
        private final double err;

        public PowerMethod(int numIterations, double err) {
            this.numIterations = numIterations;
            this.err = err;
        }

        @Override
        public double find(double[][] matrix) {

            final ReadableMatrix mat = new WrappedMatrix.ArrayOfArray(matrix);

            WrappedVector y = getInitialGuess(matrix.length);
            double maxEigenvalue = Double.NEGATIVE_INFINITY;

            for (int i = 0; i < numIterations; ++i) {

                ReadableVector v = new ReadableVector.Scale(1 / norm(y), y);
                y = product(mat, v);
                maxEigenvalue = innerProduct(v, y);

                ReadableVector diff = new ReadableVector.Sum(y,
                        new ReadableVector.Scale(-maxEigenvalue, v));

                if (ReadableVector.Utils.norm(diff) < err) {
                    break;
                }
            }

            return maxEigenvalue;
        }

        private static WrappedVector getInitialGuess(int dim) {

            double[] y = new double[dim];
            for (int i = 0; i < dim; ++i) {
                y[i] = MathUtils.nextDouble();
            }

            return new WrappedVector.Raw(y);
        }
    }
}
