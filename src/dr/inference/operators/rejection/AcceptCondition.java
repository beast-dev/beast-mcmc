/*
 * AcceptCondition.java
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

package dr.inference.operators.rejection;

import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;

public interface AcceptCondition {

    boolean satisfiesCondition(double[] values);

    enum SimpleAcceptCondition implements AcceptCondition {
        DescendingAbsoluteValue("descendingAbsoluteValue") {
            @Override
            public boolean satisfiesCondition(double[] values) {
                for (int i = 1; i < values.length; i++) {
                    if (Math.abs(values[i - 1]) < Math.abs(values[i])) {
                        return false;
                    }
                }
                return true;
            }
        },

//        DescendingAbsoluteValueSpaced("descendingAbsoluteValueSpaced") {
//            @Override
//            public boolean satisfiesCondition(double[] values) {
//                for (int i = 1; i < values.length; i++) {
//                    if (0.9 * Math.abs(values[i - 1]) < Math.abs(values[i])) {
//                        return false;
//                    }
//                }
//                return true;
//            }
//        },

        AlternatingSigns("descendingAlternatingSigns") {
            @Override
            public boolean satisfiesCondition(double[] values) {
                for (int i = 1; i < values.length; i++) {
                    Boolean signa = (values[i] > 0);
                    Boolean signb = (values[i - 1] > 0);
                    if (Math.abs(values[i - 1]) < Math.abs(values[i]) || signa == signb) {
                        return false;
                    }
                }
                return true;
            }
        },

        PositiveDefinite("positiveDefinite") {
            @Override
            public boolean satisfiesCondition(double[] values) {
                int n = (int) Math.sqrt(values.length);
                Matrix M = new Matrix(values, n, n);
                CholeskyDecomposition chol;
                try {
                    chol = new CholeskyDecomposition(M);
                } catch (IllegalDimension illegalDimension) {
                    throw new RuntimeException("Matrix must be square");
                }
                return chol.isSPD();
            }
        };

        private final String name;

        SimpleAcceptCondition(String name) {
            this.name = name;
        }

        public String getName() {
            return name;
        }

        @Override
        public abstract boolean satisfiesCondition(double[] values);
    }

}


