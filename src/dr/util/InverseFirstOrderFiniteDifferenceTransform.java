/*
 * InverseFirstOrderFiniteDifferenceTransform.java
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
import dr.util.Transform;

import static dr.util.Transform.Util.parseTransform;

public class InverseFirstOrderFiniteDifferenceTransform extends Transform.MultivariateTransform implements Reportable {

    @Override
    public String getReport() {
        return null;
    }

    private final Transform.UnivariableTransform incrementTransform;
    private final FirstOrderFiniteDifferenceTransform firstOrderFiniteDifferenceTransform;

    public InverseFirstOrderFiniteDifferenceTransform(int dim, Transform.UnivariableTransform incrementTransform) {
        super(dim);
        this.incrementTransform = incrementTransform;
        firstOrderFiniteDifferenceTransform = new FirstOrderFiniteDifferenceTransform(dim, incrementTransform);
    }

    @Override
    protected double[] transform(double[] values) {
        return firstOrderFiniteDifferenceTransform.inverse(values);
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
        return firstOrderFiniteDifferenceTransform.transform(values);
    }

    @Override
    public double[] inverse(double[] values, int from, int to, double sum) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public boolean isInInteriorDomain(double[] values) {
        return incrementTransform.isInInteriorDomain(values, 0, dim);
    }

    public String getTransformName() {
        return "inverseFirstOrderFiniteDifference";
    }

    @Override
    public double[] gradient(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    public double[] gradientInverse(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented.");
    }

//    @Override
//    protected double[] updateGradientLogDensity(double[] gradient, double[] values) {
//
//        double[] updated = new double[dim];
//        double[] transformedValues = transform(values);
//        double[] components = new double[dim];
//        for (int i = 0; i < dim; i++) {
//            components[i] = gradient[i] * incrementTransform.derivativeOfTransformWrtValue(transformedValues[i]);
//        }
//        double[] gradLogJacobian = getGradientLogJacobianInverse(values);
//        updated[dim - 1] = components[dim - 1]- gradLogJacobian[dim - 1];
//        for (int i = dim - 2; i > -1; i--) {
//            updated[i] = components[i] - components[i + 1] - gradLogJacobian[i];
//        }
//        return updated;
//
//    }

    @Override
    public double getLogJacobian(double[] values) {
        double logJacobian = 0.0;
        double s = 0.0;
        // transform is lower triangular
        for (int i = 0; i < values.length; i++) {
            s += values[i];
            logJacobian += Math.log(incrementTransform.gradientInverse(s));
        }
        // Why is this inverted?
        return -logJacobian;
    }

    @Override
    public double[] getGradientLogJacobianInverse(double[] values) {

        double[] gradLogJacobian = firstOrderFiniteDifferenceTransform.getGradientLogJacobianInverse(values);
        double[] gradient = new double[dim];

        for (int i = 0; i < dim - 1; i++) {
            gradient[i] = (gradLogJacobian[i] - gradLogJacobian[i+1]) * incrementTransform.derivativeOfTransformWrtValue(values[i]);
        }
        gradient[dim - 1] = gradLogJacobian[dim - 1] * incrementTransform.derivativeOfTransformWrtValue(values[dim - 1]);

        return gradient;
    }

    @Override
    // jacobian[j][i] = d x_i / d y_j
    public double[][] computeJacobianMatrixInverse(double[] values) {
        double[][] jacobian = new double[dim][dim];
        // Is this transposed wrong?
        for (int i = 0; i < dim - 1; i++) {
            jacobian[i][i] = incrementTransform.derivativeOfTransformWrtValue(values[i]);
            jacobian[i][i+1] = -incrementTransform.derivativeOfTransformWrtValue(values[i]);
        }
        jacobian[dim - 1][dim - 1] = incrementTransform.derivativeOfTransformWrtValue(values[dim - 1]);
        return jacobian;
    }

    //************************************************************************
    // Parser
    //************************************************************************

    public static final String NAME = "inverseFirstOrderFiniteDifferenceTransform";
    public static final String INCREMENT_TRANSFORM = "type";
    public static final String DIMENSION = "dim";

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

//            Transform transform = parseTransform(xo);
//            if (transform == null) {
//                transform = new Transform.NoTransform();
//            }
            double upper = xo.getAttribute("upper", 1.0);
            double lower = xo.getAttribute("lower", 0.0);

            String ttype = (String) xo.getAttribute(INCREMENT_TRANSFORM, "none");
            UnivariableTransform incrementTransform;
            if ( ttype.equalsIgnoreCase("none") ) {
                incrementTransform = new Transform.NoTransform();
            } else if ( ttype.equalsIgnoreCase("log") ) {
                incrementTransform = new Transform.LogTransform();
            } else if ( ttype.equalsIgnoreCase("logit") ) {
                incrementTransform = new Transform.ScaledLogitTransform(lower, upper);
            } else {
                throw new RuntimeException("Invalid option for "+ INCREMENT_TRANSFORM);
            }

            int dim = xo.getIntegerAttribute(DIMENSION);

            return new InverseFirstOrderFiniteDifferenceTransform(dim, incrementTransform);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "For vector X, creates new vector Y with Y[i] = incrementTransform^-1(sum_j=1^i X[j]).";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newIntegerRule(DIMENSION, false),
                AttributeRule.newStringRule(INCREMENT_TRANSFORM,true,"The transformation on the increments.")
//                new ElementRule(Transform.MultivariableTransformWithParameter.class, true),
        };

        public Class getReturnType() {
            return Transform.MultivariateTransform.class;
        }
    };
}