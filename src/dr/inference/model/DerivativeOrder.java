/*
 * GradientProvider2.java
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

package dr.inference.model;

/**
 * @author Marc A. Suchard
 */
public enum DerivativeOrder {

    ZEROTH("zeroth", 0) {
        @Override
        public int getDerivativeDimension(int parameterDimension) {
            return 1;
        }
    },

    GRADIENT("gradient", 1) {
        @Override
        public int getDerivativeDimension(int parameterDimension) {
            return parameterDimension;
        }
    },

    DIAGONAL_HESSIAN("diagonalHessian", 2) {
        @Override
        public int getDerivativeDimension(int parameterDimension) {
            return parameterDimension;
        }
    },

    FULL_HESSIAN("Hessian", 3) {
        @Override
        public int getDerivativeDimension(int parameterDimension) {
            return parameterDimension * parameterDimension;
        }
    };

    DerivativeOrder(String type, int order) {
        this.type = type;
        this.order = order;
    }

    abstract public int getDerivativeDimension(int parameterDimension);

    final private String type;
    final int order;

    final public String toString() {
        return type;
    }

    final public int getValue() {
        return order;
    }
}
