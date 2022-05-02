/*
 * EuclideanBallToRTransform.java
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

package dr.util;

import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

import static dr.util.EuclideanToInfiniteNormUnitBallTransform.squaredNorm;

/**
 * @author Marc A. Suchard
 * @author Paul Bastide
 */

public class EuclideanBallToRTransform extends Transform.MultivariateTransform {
    /**
     * Transform a vector $x from the euclidean norm unit ball (||x||_2 <= 1.0)
     * to an unconstrained vector of R^n.
     */

    public EuclideanBallToRTransform(int dim) {
        super(dim);
    }

    // values = vector of euclidean unit ball
    @Override
    protected double[] transform(double[] values) {
        assert isInEuclideanUnitBall(values) : "Initial vector is not in the Euclidean unit ball.";

        double[] transformedValues = new double[values.length];
        double factor = Math.pow(1.0 - squaredNorm(values), -0.5);
        for (int i = 0; i < values.length; i++) {
            transformedValues[i] = values[i] * factor;
        }

        return transformedValues;
    }

    // values = vector of R^n
    @Override
    protected double[] inverse(double[] transformedValues) {
        double[] values = new double[transformedValues.length];
        double factor = Math.pow(1.0 + squaredNorm(transformedValues), -0.5);
        for (int i = 0; i < transformedValues.length; i++) {
            values[i] = transformedValues[i] * factor;
        }

        return values;
    }

    private boolean isInEuclideanUnitBall(double[] x) {
        return (squaredNorm(x) <= 1.0);
    }

    private boolean isInStrictEuclideanUnitBall(double[] x) {
        return (squaredNorm(x) <= 1.0);
    }

    @Override
    public boolean isInInteriorDomain(double[] values) {
        return isInStrictEuclideanUnitBall(values);
    }

    @Override
    public double[] inverse(double[] values, int from, int to, double sum) {
        throw new RuntimeException("Not relevant.");
    }

    public String getTransformName() {
        return "EuclideanBallToRTransform";
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
        return -0.5 * (dim + 2) * Math.log(1.0 - squaredNorm(values));
    }

    @Override
    protected double[] getGradientLogJacobianInverse(double[] transformedValues) {
        throw new RuntimeException("Not yet implemented");
//        double[] gradientLogJacobian = new double[transformedValues.length];
//        double factor = 1.0 + squaredNorm(transformedValues);
//        for (int i = 0; i < dim - 2; i++) {
//            gradientLogJacobian[i] = 2 * (dim - 1) * transformedValues[i] / factor;
//        }
//        return gradientLogJacobian;
    }

    // ************************************************************************* //
    // Computation of the jacobian matrix
    // ************************************************************************* //

    // Returns the *transpose* of the Jacobian matrix: jacobian[j][i] = d x_i / d y_j
    @Override
    public double[][] computeJacobianMatrixInverse(double[] transformedValues) {
        double[][] jacobian = new double[dim][dim];
        double cst = 1.0 + squaredNorm(transformedValues);
        double factor = Math.pow(cst, -1.5);

        for (int j = 0; j < dim; j++) {
            jacobian[j][j] = (cst - transformedValues[j] * transformedValues[j]) * factor;
            for (int i = j + 1; i < dim; i++) {
                jacobian[j][i] = -transformedValues[j] * transformedValues[i] * factor;
                jacobian[i][j] = jacobian[j][i];
            }
        }

        return jacobian;
    }

    //************************************************************************
    // Parser
    //************************************************************************

    public static final String NAME = "sphericalTransform2";
    public static final String DIMENSION = "dim";

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            int dim = xo.getIntegerAttribute(DIMENSION);

            // Spherical (Euclidean to unconstrained)
            List<MultivariableTransform> transformsMul = new ArrayList<MultivariableTransform>();
            for (int i = 0; i < dim + 1; i++) {
                transformsMul.add(new EuclideanBallToRTransform(dim));
            }

            return new Transform.MultivariateArray(transformsMul);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A spherical transform.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newIntegerRule(DIMENSION, false),
        };

        public Class getReturnType() {
            return Transform.MultivariateArray.class;
        }
    };

}

