/*
 * GradientWrtParameterProvider.java
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

import java.util.HashMap;
import java.util.Map;

/**
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */
public interface PrecisionColumnProvider {

    double[] getColumn(int index);

    @SuppressWarnings("unused")
    class Generic extends AbstractModel implements PrecisionColumnProvider {

        private final MatrixParameterInterface matrix;
        private final Map<Integer, double[]> cache = new HashMap<>();

        public Generic(MatrixParameterInterface matrix) {
            super("precisionColumnProvider.Generic");
            this.matrix = matrix;

            addVariable(matrix);
        }

        @Override
        public double[] getColumn(int index) {

            if (!cache.containsKey(index)) {

                final int dim = matrix.getRowDimension();

                double[] column = new double[dim];
                for (int row = 0; row < dim; ++row) {
                    column[row] = matrix.getParameterValue(row, index);
                }

                cache.put(index, column);
            }

            return cache.get(index);
        }

        @Override
        protected void handleModelChangedEvent(Model model, Object object, int index) { }

        @Override
        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            if (variable == matrix) {
                cache.clear();
            }
        }

        @Override
        protected void storeState() { } // TODO

        @Override
        protected void restoreState() { } // TODO

        @Override
        protected void acceptState() { }
    }
}