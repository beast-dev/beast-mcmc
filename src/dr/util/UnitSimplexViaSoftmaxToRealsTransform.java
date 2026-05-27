/*
 * UnitSimplexViaSoftmaxToRealsTransform.java
 *
 * Copyright © 2002-2025 the BEAST Development Team
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

package dr.util;

import dr.xml.*;

/**
 * @author Marc Suchard
 * @author Xinghua Tao
 *
 * @notes <a href="https://mc-stan.org/docs/reference-manual/transforms.html#simplex-transform.section"/>
 */

public class UnitSimplexViaSoftmaxToRealsTransform extends Transform.MultivariateTransform {

    private final SumToZeroTransform sumToZeroTransform;

    public UnitSimplexViaSoftmaxToRealsTransform(int dim) {
        super(dim);

        this.sumToZeroTransform = new SumToZeroTransform(dim);
    }

    @Override
    protected double[] transform(double[] valuesOnSimplex) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected double[] inverse(double[] valuesOnReals) {

        double[] z = sumToZeroTransform.inverse(valuesOnReals);

        double total = sumExp(z, dim);
        double[] valuesOnSimplex = new double[dim];
        for (int i = 0; i < dim; ++i) {
            valuesOnSimplex[i] = Math.exp(z[i]) / total;
        }

        return valuesOnSimplex;
    }

    @Override
    public double[] inverse(double[] values, int from, int to, double sum) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public boolean isInInteriorDomain(double[] values) {
        double total = 0.0;
        for (double v : values) {
            if (v < 0.0 || v > 1.0) {
                return false;
            }
            total += v;
        }
        return Math.abs(total - 1.0) < 1E-6;
    }

    public String getTransformName() {
        return "UnitSimplexTransform";
    }

    @Override
    public double[] gradient(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] gradientInverse(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double getLogJacobian(double[] valuesOnSimplex) {

        double logDet = sumToZeroTransform.getLogJacobian(null);

        for (double value : valuesOnSimplex) {
            logDet += Math.log(value);
        }

        return logDet;
    }

    @Override
    public double[] getGradientLogJacobianInverse(double[] valuesOnReals) {

        double[] valuesOnSimplex = inverse(valuesOnReals);
//        double[][] computeJacobianMatrixInverse

        double[][] partial = new double[1][dim];
        for (int i = 0; i < dim; ++i) {
            partial[0][i] = 1 / valuesOnSimplex[i];
        }

        double[][] chain = sumToZeroTransform.computeJacobianMatrixInverse(valuesOnReals);
//        return multiply(jacobian, false, chain, true, dim, dim,  dim -1 ); // TODO compute everything in place!
        double[][] product = multiply(partial, false, chain, true, 1, dim, dim - 1);

        return product[0];

//        return chain;
    }

    // ************************************************************************* //
    // Computation of the (transposed) Jacobian matrix
    // ************************************************************************* //

    private static final boolean K_PLUS_ONE = true;

    public double[][] computeJacobianMatrixInverse(double[] valuesOnReals) {

        double[] valuesOnSimplex = inverse(valuesOnReals); // also computes z[]

        double[][] jacobian;
        if (K_PLUS_ONE) {
            jacobian = new double[dim][dim];
        } else {
            jacobian = new double[dim - 1][dim - 1];
        }

        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                if (i == j) {
                    jacobian[i][i] = valuesOnSimplex[i] * (1.0 - valuesOnSimplex[i]);
                } else {
                    jacobian[i][j] = -valuesOnSimplex[i] * valuesOnSimplex[j];
                }
            }
        }
//        if (K_PLUS_ONE) {
//            jacobian[dim - 1][dim - 1] = 1.0;
//        }

//        return jacobian;
        double[][] chain = sumToZeroTransform.computeJacobianMatrixInverse(valuesOnReals);
//        return multiply(jacobian, false, chain, true, dim, dim,  dim -1 ); // TODO compute everything in place!
        return multiply(chain, false, jacobian, false, dim - 1, dim, dim);
//        return multiply(chain, true, jacobian, false, dim, dim - 1, dim - 1);
    }

    private static double[][] multiply(double[][] a, boolean transposeA, double[][] b, boolean transposeB,
                                       int N, int K, int M) {
        double[][] result = new double[N][M];
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                for (int k = 0; k < K; ++k) {
                    double aValue = (!transposeA) ? a[i][k] : a[k][i];
                    double bValue = (!transposeB) ? b[k][j] : b[j][k];
                    result[i][j] += aValue * bValue; // TODO B is symmetric; remove transpose here
                }
            }
        }
        return result;
    }

    private static double sumExp(double[] values, int end) {
        double total = 0.0;
        for (int i = 0; i < end; ++i) {
            total += Math.exp(values[i]);
        }
        return total;
    }

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        private static final String DIMENSION = "dimension";
        private static final String UNIT_SIMPLEX_VIA_SOFTMAX_TRANSFORM = "UnitSimplexViaSoftmaxTransform";

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            int dim = xo.getIntegerAttribute(DIMENSION);
            return new UnitSimplexViaSoftmaxToRealsTransform(dim);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    AttributeRule.newIntegerRule(DIMENSION)
            };
        }

        @Override
        public String getParserDescription() {
            return "transforms values on the simplex to values on the real-line";
        }

        @Override
        public Class getReturnType() {
            return UnitSimplexViaSoftmaxToRealsTransform.class;
        }

        @Override
        public String getParserName() {
            return UNIT_SIMPLEX_VIA_SOFTMAX_TRANSFORM;
        }
    };
}

