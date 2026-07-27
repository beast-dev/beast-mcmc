/*
 * UnitSimplexToRealsTransform.java
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

import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.xml.*;

/**
 * @author Marc Suchard
 * @author Xinghua Tao
 *
 * @notes <a href="https://mc-stan.org/docs/reference-manual/transforms.html#sum-to-zero-transforms"/>
 */

public class SumToZeroTransform extends Transform.MultivariateTransform {

    private double[][] constantJacobian;
    private double logDeterminant = Double.NaN;

    public SumToZeroTransform(int sumToZeroDim) {
        super(sumToZeroDim, sumToZeroDim - 1);
    }

    @Override
    protected double[] transform(double[] valuesOnSimplex) {
        // x \in \Re^{N + 1} --> y \in \Re^{N} (unconstrained)

        assert valuesOnSimplex.length == getInputDimension();


        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected double[] inverse(double[] unconstrainedValues) {
        // y \in \Re^{N} --> x \in \Re^{N + 1}

        assert unconstrainedValues.length == getOutputDimension();

        double[] sumToZeroValues = new double[getInputDimension()];

        int N = dim;

        double sum = 0.0;
        for (int n = 1; n <= N; ++n) {
            sum += unconstrainedValues[n - 1] / Math.sqrt(n * (n + 1));
        }
        sumToZeroValues[1 - 1] = sum;

        for (int n = 1; n <= N; ++n) {
            sum = 0.0;
            for (int i = n + 1; i <= N; ++i) {
                sum += unconstrainedValues[i - 1] / Math.sqrt(i * (i + 1));
            }
            sumToZeroValues[n + 1 - 1] = sum - n * unconstrainedValues[n - 1] / Math.sqrt(n * (n + 1));
        }

        double[] z = new double[dim + 1];
        double sumW = 0.0;
        for (int n = 1; n <= N; ++n) {
            int i = N - n + 1;
            double w = unconstrainedValues[i - 1] * inv_sqrt(i * (i + 1));
            sumW += w;
            z[i - 1] += sumW;
            z[i + 1 - 1] -= w * i;
        }

        return sumToZeroValues;
    }




    private static double inv_sqrt(double x) {
        return 1 / Math.sqrt(x);
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

        if (Double.isNaN(logDeterminant)) {
            double[][] tmp = computeJacobianMatrixInverse(null);
            double[][] jacobian = new double[tmp.length + 1][];
            System.arraycopy(tmp, 0, jacobian, 0, tmp.length);
            jacobian[tmp.length] = new double[tmp.length + 1];
            jacobian[tmp.length][tmp.length] = 1.0;

            Matrix J = new Matrix(jacobian);

            try {
                logDeterminant = J.logDeterminant();
            } catch (IllegalDimension e) {
                throw new RuntimeException("Not yet implemented");
            }
        }

        return logDeterminant;
    }

    @Override
    public double[] getGradientLogJacobianInverse(double[] valuesOnReals) {
        return new double[outputDimension];
    }

    // ************************************************************************* //
    // Computation of the (transposed) Jacobian matrix
    // ************************************************************************* //

    public double[][] computeJacobianMatrixInverse(double[] valuesOnReals) {

        if (constantJacobian == null) {
            double[][] jacobian = new double[outputDimension][inputDimension];

            for (int i = 0; i < outputDimension; ++i) {
                double factor = 1.0 / Math.sqrt((i + 1) * (i + 2));
                for (int j = 0; j <= i; ++j) {
                    jacobian[i][j] = factor;
                }
                jacobian[i][i + 1] = -(i + 1) * factor;
            }
            constantJacobian = jacobian;
        }

        return constantJacobian;
    }

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        private static final String DIMENSION = "dimension";
        private static final String UNIT_SIMPLEX_VIA_SOFTMAX_TRANSFORM = "SumToZeroTransform";

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            int dim = xo.getIntegerAttribute(DIMENSION);
            return new SumToZeroTransform(dim);
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
            return SumToZeroTransform.class;
        }

        @Override
        public String getParserName() {
            return UNIT_SIMPLEX_VIA_SOFTMAX_TRANSFORM;
        }
    };
}

