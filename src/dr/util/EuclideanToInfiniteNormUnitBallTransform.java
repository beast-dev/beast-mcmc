/*
 * EuclideanToInfiniteNormBallTransform.java
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

package dr.util;

import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Paul Bastide
 */

public class EuclideanToInfiniteNormUnitBallTransform extends Transform.MultivariateTransform {
    /**
     * Transform a vector $x from the euclidean norm unit ball (||x||_2 <= 1.0)
     * to a vector of the infinite norm unit ball (max(|x_i|) <= 1.0).
     * Uses a transformation similar to the one in Lewandowski, Kurowicka, and Joe (2009).
     */

    public EuclideanToInfiniteNormUnitBallTransform(int dim) {
        super(dim);
    }

    // values = vector of euclidean unit ball
    @Override
    protected double[] transform(double[] values) {
        assert isInEuclideanUnitBall(values) : "Initial vector is not in the Euclidean unit ball.";

        double[] transformedValues = new double[values.length];

        double acc = 1.0;
        double temp;
        for (int i = 0; i < dim; i++) {
            temp = values[i] / acc;
            transformedValues[i] = temp;
            acc *= Math.sqrt(1 - temp * temp);
        }

        return transformedValues;
    }

    // values = vector of invinite unit ball
    @Override
    protected double[] inverse(double[] values) {
        assert isInInfiniteUnitBall(values) : "Initial vector is not in the Euclidean unit ball.";

        double[] transformedValues = new double[values.length];

        double acc = 1.0;
        double temp;
        for (int i = 0; i < dim; i++) {
            temp = values[i];
            transformedValues[i] = temp * acc;
            acc *= Math.sqrt(1 - temp * temp);
        }

        return transformedValues;
    }

    private boolean isInEuclideanUnitBall(double[] x) {
        return (squaredNorm(x) <= 1.0);
    }

    private boolean isInInfiniteUnitBall(double[] x) {
        for (int k = 0; k < dim; k++) {
            if (!(x[k] <= 1.0 && x[k] >= -1.0)) return false;
        }
        return true;
    }

    public static double squaredNorm(double[] x) {
        return squaredNorm(x, 0, x.length);
    }

    public static double squaredNorm(double[] x, int offset, int length) {
        double norm = 0.0;
        for (int i = 0; i < length; i++) {
            norm += x[offset + i] * x[offset + i];
        }
        return norm;
    }

    public static double projection(final double[] x) {
        return projection(x, 0, x.length);
    }

    public static double projection(final double[] x, final int offset, final int length) {
        return Math.sqrt(1.0 - squaredNorm(x, offset, length));
    }

    @Override
    public double[] inverse(double[] values, int from, int to, double sum) {
        throw new RuntimeException("Not relevant.");
    }

    public String getTransformName() {
        return "EuclideanToInfiniteNormUnitBallTransform";
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
    protected double getLogJacobian(double[] values) {
        double[] transformedValues = transform(values);

        double logJacobian = 0;
        for (int i = 0; i < dim - 1; i++) {
            logJacobian += (dim - i - 1) * Math.log(1.0 - Math.pow(transformedValues[i], 2));
        }
        return -0.5 * logJacobian;
    }

    @Override
    protected double[] getGradientLogJacobianInverse(double[] values) {
        double[] gradientLogJacobian = new double[values.length];
        for (int i = 0; i < dim - 1; i++) { // Sizes of conditioning sets
            gradientLogJacobian[i] = -(dim - i - 1) * values[i] / (1.0 - Math.pow(values[i], 2));
        }
        return gradientLogJacobian;
    }

    // ************************************************************************* //
    // Computation of the jacobian matrix
    // ************************************************************************* //

    // Returns the *transpose* of the Jacobian matrix: jacobian[j][i] = d x_i / d y_j
    @Override
    public double[][] computeJacobianMatrixInverse(double[] values) {
        double[][] jacobian = new double[dim][dim];

        for (int j = 0; j < dim; j++) {
            double acc = 1.0;
            double temp;
            // before j
            for (int i = 0; i < j; i++) {
                acc *= Math.sqrt(1 - Math.pow(values[i], 2));
            }
            // j
            jacobian[j][j] = acc;
            temp = values[j];
            acc *= -temp / Math.sqrt(1 - Math.pow(temp, 2));
            // after j
            for (int i = j + 1; i < dim; i++) {
                temp = values[i];
                jacobian[j][i] = temp * acc;
                acc *= Math.sqrt(1 - Math.pow(temp, 2));
            }
        }

        return jacobian;
    }

    //************************************************************************
    // Parser
    //************************************************************************

    public static final String NAME = "sphericalTransform";
    public static final String DIMENSION = "dim";
//    private static final String IS_MATRIX = "isMatrix";

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            int dim = xo.getIntegerAttribute(DIMENSION);

            // Fisher Z  (Infinte norm unit ball to unconstrained)
            List<Transform> transforms = new ArrayList<Transform>();
            for (int i = 0; i < dim * (dim + 1); i++) {
                transforms.add(Transform.FISHER_Z);
            }
            Transform.Array fisherZTransforms = new Transform.Array(transforms, null);

            // Spherical (Euclidean to Infinite norm unit ball)
            List<MultivariateTransform> transformsMul = new ArrayList<MultivariateTransform>();
            for (int i = 0; i < dim + 1; i++) {
                transformsMul.add(new EuclideanToInfiniteNormUnitBallTransform(dim));
            }
            Transform.MultivariateTransform sphericalTransform = new Transform.MultivariateArray(transformsMul);

            // Compose
            return new Transform.ComposeMultivariable(fisherZTransforms, sphericalTransform);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A spherical transform using Fisher Z and LKJ.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newIntegerRule(DIMENSION, false),
//                AttributeRule.newBooleanRule(IS_MATRIX, true)
        };

        public Class getReturnType() {
            return Transform.ComposeMultivariable.class;
        }
    };

}
