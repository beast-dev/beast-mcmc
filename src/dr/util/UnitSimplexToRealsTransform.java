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

import dr.xml.*;

/**
 * @author Marc Suchard
 * @author Xinghua Tao
 *
 * @notes <a href="https://mc-stan.org/docs/2_19/reference-manual/simplex-transform-section.html"/>
 */

public class UnitSimplexToRealsTransform extends Transform.MultivariateTransform {

    private final double[] z;
    private final double[] logOffsets;

    public UnitSimplexToRealsTransform(int dim) {
        super(dim);

        this.z = new double[dim - 1];
        this.logOffsets = new double[dim];

        for (int i = 0; i < dim - 1; ++i) {
            logOffsets[i] = -Math.log(dim - i - 1);
        }
    }

    @Override
    protected double[] transform(double[] valuesOnSimplex) {

        computeZFromValuesOnSimplex(valuesOnSimplex);

        double[] valuesOnReals = new double[dim];
        for (int i = 0; i < dim - 1; ++i) {
            valuesOnReals[i] = logit(z[i]) - logOffsets[i];
        }

        return valuesOnReals;
    }

    @Override
    protected double[] inverse(double[] valuesOnReals) {

        computeZFromValuesOnReals(valuesOnReals);

        double stickRemainder = 1.0;
        double[] valuesOnSimplex = new double[dim];

        valuesOnSimplex[0] = z[0];
        for (int i = 1; i < dim - 1; ++i) {
            stickRemainder -= valuesOnSimplex[i - 1];
            valuesOnSimplex[i] = stickRemainder * z[i];
        }
        stickRemainder -= valuesOnSimplex[dim - 2];
        valuesOnSimplex[dim - 1] = stickRemainder;

        return valuesOnSimplex;
    }
    
    private void computeZFromValuesOnSimplex(double[] valuesOnSimplex) {
        double stickRemainder = 1.0;
        z[0] = valuesOnSimplex[0];
        for (int i = 1; i < dim - 1; ++i) {
            stickRemainder -= valuesOnSimplex[i - 1];
            z[i] = valuesOnSimplex[i] / stickRemainder;
        }
    }

    private void computeZFromValuesOnReals(double[] valuesOnReals) {
        for (int i = 0; i < dim - 1; ++i) {
            z[i] = logitInverse(valuesOnReals[i] + logOffsets[i]);
        }
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

        computeZFromValuesOnSimplex(valuesOnSimplex);

        double logDet = 0.0;

        double stickRemainder = 1.0;
        for (int i = 0; i < dim - 1; ++i) {
            logDet += Math.log(z[i]) + Math.log(1.0 - z[i]) + Math.log(stickRemainder);
            stickRemainder -= valuesOnSimplex[i];
        }

        return -logDet;
    }

    @Override
    public double[] getGradientLogJacobianInverse(double[] valuesOnReals) {

        double[] valuesOSimplex = inverse(valuesOnReals); // also computes z[]
        double[][] jacobian = computeJacobianMatrixInverse(valuesOnReals);

        double[] gradient = new double[dim];
        for (int i = 0; i < dim - 1; ++i) {
            gradient[i] = 1 - 2 * z[i];
        }

        double stickRemainder = 1.0;
        double[] acc = new double[dim];
        for (int i = 0; i < dim - 2; ++i) {
            stickRemainder -= valuesOSimplex[i];

            for (int j = 0; j <= i; ++j) {
                acc[j] += jacobian[j][i];
                gradient[j] -= acc[j] / stickRemainder;
            }
        }

        return gradient;
    }

    // ************************************************************************* //
    // Computation of the (transposed) Jacobian matrix
    // ************************************************************************* //

    public double[][] computeJacobianMatrixInverse(double[] valuesOnReals) {

        double[] valuesOnSimplex = inverse(valuesOnReals); // also computes z[]

        double[][] jacobian = new double[dim][dim];
        double[] accumulativeDifferential = new double[dim];

        double stickRemainder = 1.0;
        for (int i = 0; i < dim - 1; ++i) {

            for (int j = 0; j < i; ++j) {
                jacobian[j][i] = -z[i] * accumulativeDifferential[j];
            }
            jacobian[i][i] = stickRemainder * z[i] * (1.0 - z[i]);

            for (int j = 0; j <= i; ++j) {
                accumulativeDifferential[j] += jacobian[j][i];
            }

            stickRemainder -= valuesOnSimplex[i];
        }

        for (int i = 0; i < dim - 1; ++i) {
            jacobian[i][dim - 1] = -accumulativeDifferential[i];
        }

        return jacobian;
    }

    private static double logit(double x) {
        return Math.log(x / (1.0 - x));
    }

    private static double logitInverse(double y) {
        double expY = Math.exp(y);
        return expY / (1.0 + expY);
    }

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        private static final String DIMENSION = "dimension";
        private static final String UNIT_SIMPLEX_TRANSFORM = "UnitSimplexTransform";

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            int dim = xo.getIntegerAttribute(DIMENSION);
            return new UnitSimplexToRealsTransform(dim);
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
            return UnitSimplexToRealsTransform.class;
        }

        @Override
        public String getParserName() {
            return UNIT_SIMPLEX_TRANSFORM;
        }
    };
}

