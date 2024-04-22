/*
 * FirstOrderFiniteDifferenceTransform.java
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

public class FirstOrderFiniteDifferenceTransform extends Transform.MultivariateTransform implements Reportable {

    @Override
    public String getReport() {
        return null;
    }

    private final Transform.UnivariableTransform incrementTransform;

    public FirstOrderFiniteDifferenceTransform(int dim, Transform.UnivariableTransform incrementTransform) {
        super(dim);
        this.incrementTransform = incrementTransform;
    }

    @Override
    protected double[] transform(double[] values) {
        double[] incrementTransformedValues = incrementTransform.transform(values, 0, dim);
        double[] increments = new double[values.length];
        increments[0] = incrementTransformedValues[0];
        for (int i = 1; i < values.length; i++) {
            increments[i] = incrementTransformedValues[i] - incrementTransformedValues[i - 1];
        }
        return increments;
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
        double[] fx = new double[values.length];
        fx[0] = incrementTransform.inverse(values[0]);
        double s = values[0];
        for (int i = 1; i < values.length; i++) {
            s += values[i];
            fx[i] = incrementTransform.inverse(s);
        }
        return fx;
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
        return "firstOrderFiniteDifference";
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
    // This gets called on the untransformed values to update the gradient to the transformed space
    // values: the untransformed values
    // gradient: the gradient wrt the (untransformed) values
    // see for example updateMomentum(), line 789 of HamiltonianMonteCarloOperator.java
    protected double[] updateGradientLogDensity(double[] gradient, double[] values) {
        double[] updated = new double[dim];
        double[] transformedValues = transform(values);

        double[] jacobianInverseFirstRow = new double[dim];
        double s = 0.0;
        for (int i = 0; i < dim; i++) {
            s += transformedValues[i];
            jacobianInverseFirstRow[i] = incrementTransform.gradientInverse(s);
        }

        updated[dim - 1] = gradient[dim - 1] * jacobianInverseFirstRow[dim - 1];
        for (int i = dim - 2; i > -1; i--) {
            updated[i] = gradient[i] * jacobianInverseFirstRow[i] + updated[i+1];
        }

        double[] gradLogJacobian = getGradientLogJacobianInverse(transformedValues);
        for (int i = dim - 1; i > -1; i--) {
            updated[i] += gradLogJacobian[i];
        }

        return updated;

    }

    @Override
    public double getLogJacobian(double[] values) {
        double logJacobian = 0.0;
        // Transform is lower triangular
        for (int i = 0; i < values.length; i++) {
            logJacobian += Math.log(incrementTransform.derivativeOfTransformWrtValue(values[i]));
        }
        return logJacobian;
    }

//    @Override
//    // called on the transformed values
//    public double[] getGradientLogJacobianInverse(double[] values) {
//        // jacobianDiagonal == diagonal of Jacobian of inverse transform
//        double[] jacobianDiagonal = new double[values.length];
//        double s = 0.0;
//        double[] cumSum = new double[dim];
//        for (int i = 0; i < values.length; i++) {
//            s += values[i];
//            cumSum[i] = s;
//            jacobianDiagonal[i] = incrementTransform.gradientInverse(s);
//        }
//        double[] gradient = new double[values.length];
//        double tmp = 0.0;
//        for (int i = values.length - 1; i > -1; i--) {
//            tmp += (1.0 / jacobianDiagonal[i]) * incrementTransform.secondDerivativeOfInverseTransformWrtValue(cumSum[i]);
//            gradient[i] = tmp;
//        }
//
//        return gradient;
//    }

    @Override
    // called on the transformed values
    public double[] getGradientLogJacobianInverse(double[] values) {
        // jacobianDiagonal == diagonal of Jacobian of inverse transform
        double[] logJacobianDiagonal = new double[values.length];
        double s = 0.0;
        double[] cumSum = new double[dim];
        for (int i = 0; i < values.length; i++) {
            s += values[i];
            cumSum[i] = s;
            logJacobianDiagonal[i] = incrementTransform.logGradientInverse(s);
        }
        double[] gradient = new double[values.length];
        double tmp = 0.0;
        for (int i = values.length - 1; i > -1; i--) {
            // Avoid Inf/Inf errors
            tmp += Math.exp(-logJacobianDiagonal[i] + incrementTransform.logSecondDerivativeOfInverseTransformWrtValue(cumSum[i]));
            gradient[i] = tmp;
        }

        return gradient;
    }
    @Override
    // Gets called on *transformed* values, as per Transform.updateGradientInverseUnWeightedLogDensity()
    // jacobian[j][i] = d x_i / d y_j
    public double[][] computeJacobianMatrixInverse(double[] values) {
        double[][] jacobian = new double[dim][dim];
        // Is this transposed wrong?
        double s = 0.0;
        for (int i = 0; i < dim; i++) {
            s += values[i];
            for (int j = 0; j < i + 1; j++) {
                jacobian[j][i] = incrementTransform.gradientInverse(s);
            }
        }
        return jacobian;
    }

    //************************************************************************
    // Parser
    //************************************************************************

    public static final String NAME = "firstOrderFiniteDifferenceTransform";
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

            return new FirstOrderFiniteDifferenceTransform(dim, incrementTransform);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "For vector X, creates new vector Y with Y[0] = incrementTransform(X[0]), Y[i] = incrementTransform(X[i]) - incrementTransform(X[i - i]).";
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