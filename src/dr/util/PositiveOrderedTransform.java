/*
 * PositiveOrdered.java
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

public class PositiveOrderedTransform extends Transform.MultivariateTransform {

    public PositiveOrderedTransform(int dim) {
        super(dim);
    }

    @Override
    // x (positive ordered) -> y (unconstrained)
    protected double[] transform(double[] values) {
        double[] result = new double[dim];
        result[0] = Math.log(values[0]);
        for (int i = 1; i < dim; i++) {
            result[i] = Math.log(values[i] - values[i - 1]);
        }
        return result;
    }

    public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value, int from, int to) {
        throw new RuntimeException("not implemented yet");
    }

    @Override
    public double[][] updateHessianLogDensity(double[][] hessian, double[][] transformationHessian, double[] gradient, double[] value, int from, int to) {
        throw new RuntimeException("not implemented yet");
    }

    @Override
    public double[] updateGradientInverseUnWeightedLogDensity(double[] gradient, double[] value, int from, int to) {
        throw new RuntimeException("not implemented yet");
    }

    @Override
    public double[] updateGradientUnWeightedLogDensity(double[] gradient, double[] value, int from, int to) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected double[] inverse(double[] values) {
        double[] result = new double[dim];
        result[0] = Math.exp(values[0]);
        for (int i = 1; i < dim; i++) {
            result[i] = result[i - 1] + Math.exp(values[i]);
        }
        return result;
    }

    @Override
    public double[] inverse(double[] values, int from, int to, double sum) {
        throw new RuntimeException("Not relevant.");
    }

    public String getTransformName() {
        return "PositiveOrdered";
    }

    @Override
    public double[] gradient(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    public double[] gradientInverse(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    protected double getLogJacobian(double[] values) {
        double result = Math.log(values[0]);
        for (int i = 1; i < dim; i++) {
            result += Math.log(values[i] - values[i - 1]);
        }
        return -result;
    }

    @Override
    public double[] getGradientLogJacobianInverse(double[] values) {
        double[] result = new double[dim];
        for (int i = 0; i < dim; i++) {
            result[i] = 1.0;
        }
        return result;
    }

    @Override
    // jacobian[j][i] = d x_i / d y_j
    public double[][] computeJacobianMatrixInverse(double[] values) {
        double[][] jacobian = new double[dim][dim];
        for (int i = 0; i < dim; i++) {
            for (int j = i; j < dim; j++) {
                jacobian[j][i] = Math.exp(values[j]);
            }
        }
        return jacobian;
    }

    //************************************************************************
    // Parser
    //************************************************************************

    public static final String NAME = "positiveOrderedTransform";
    public static final String DIMENSION = "dim";
//    private static final String IS_MATRIX = "isMatrix";

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            int dim = xo.getIntegerAttribute(DIMENSION);

            return new PositiveOrderedTransform(dim);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A compound matrix parametrized by its eigen values and vectors.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newIntegerRule(DIMENSION, false)
        };

        public Class getReturnType() {
            return Transform.MultivariateTransform.class;
        }
    };
}