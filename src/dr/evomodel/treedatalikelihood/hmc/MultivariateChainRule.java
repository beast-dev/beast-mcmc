/*
 * MultivariateChainRule.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

public interface MultivariateChainRule {

    double[] chainGradient(double[] lhs);

    void chainGradient(DenseMatrix64F gradient);

    class Chain implements MultivariateChainRule {

        private final MultivariateChainRule[] rules;

        Chain(MultivariateChainRule[] rules) {
            this.rules = rules;
        }

        @Override
        public double[] chainGradient(double[] gradient) {

            for (MultivariateChainRule rule : rules) {
                gradient = rule.chainGradient(gradient);
            }
            return gradient;
        }

        @Override
        public void chainGradient(DenseMatrix64F gradient) {
            for (MultivariateChainRule rule : rules) {
                rule.chainGradient(gradient);
            }
        }
    }

    class Inverse implements MultivariateChainRule {

        private final double[] vecP;
        private final double[] vecV;
        private final int dim;

        Inverse(double[] vecP, double[] vecV) {
            this.vecP = vecP;
            this.vecV = vecV;
            this.dim = (int) Math.sqrt(vecP.length);
        }

        @Override
        public double[] chainGradient(double[] lhs) {

            assert lhs.length == dim * dim;

            double[] gradient = new double[dim * dim];

            for (int i = 0; i < dim * dim; ++i) {

                if (vecV[i] == 0 || Double.isNaN(vecV[i])) {
                    throw new RuntimeException("0 or NaN value in variance. check start value or use smaller step size for hmc");
                }
                gradient[i] = -lhs[i] * vecP[i] / vecV[i];
            }

            return gradient;
        }

        @Override
        public void chainGradient(DenseMatrix64F gradient) {
            throw new RuntimeException("not yet implemented");
        }
    }

    class InverseGeneral implements MultivariateChainRule {

        private final DenseMatrix64F Mat;
        private final DenseMatrix64F temp;
        private final int dim;


        public InverseGeneral(double[] vecMat) {
            this.dim = (int) Math.sqrt(vecMat.length);
            this.Mat = DenseMatrix64F.wrap(dim, dim, vecMat);
            this.temp = new DenseMatrix64F(dim, dim);
        }

        public InverseGeneral(DenseMatrix64F Mat) {
            this.dim = Mat.getNumCols();
            assert dim == Mat.getNumRows() : "Inverse rule is only valid for square matrices.";
            this.Mat = Mat;
            this.temp = new DenseMatrix64F(dim, dim);
        }

        @Override
        public double[] chainGradient(double[] lhs) {

            assert lhs.length == dim * dim;

            DenseMatrix64F gradient = new DenseMatrix64F(dim, dim);

            DenseMatrix64F LHS = DenseMatrix64F.wrap(dim, dim, lhs);
            CommonOps.mult(Mat, LHS, temp);
            CommonOps.mult(-1, temp, Mat, gradient);

            return gradient.getData();
        }

        @Override
        public void chainGradient(DenseMatrix64F gradient) {
            CommonOps.mult(Mat, gradient, temp);
            CommonOps.mult(-1, temp, Mat, gradient);
        }
    }
}
