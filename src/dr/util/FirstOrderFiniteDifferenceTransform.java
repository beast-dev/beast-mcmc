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

    //TODO: seems like we should be able to replace this with use of Transform.java functions?
    public enum IncrementTransform {
        LOG("log") {
            public double transform(double x, double upper, double lower) {
                return Math.log(x);
            }

            public double inverse(double x, double upper, double lower) {
                return Math.exp(x);
            }

//            public double derivativeOfTransform(double x, double upper, double lower) {
//                return 1.0 / x;
//            }

            public double grad(double x, double upper, double lower) {
                return x;
            }

            public double updateGradientLogDensity(double gradient, double value, double upper, double lower) {
                return gradient * value + 1.0;
            }

            public double derivativeOfTransform(double x, double upper, double lower) {
                return 1.0 / x;
            }

            public double secondDerivativeOfTransform(double x, double upper, double lower) {
                return 1.0 / (x * x);
            }

            public double derivativeOfInverseTransform(double x, double upper, double lower) {
                return Math.exp(x);
            }

            public double secondDerivativeOfInverseTransform(double x, double upper, double lower) {
                return Math.exp(x);
            }

            public boolean isInteriorDomain(double x, double upper, double lower) {
                return x >= 0;
            }

        },
        LOGIT("logit") {
            public double transform(double x, double upper, double lower) {
                throw new RuntimeException("not yet implemented");
            }

            public double inverse(double x, double upper, double lower) {
                throw new RuntimeException("not yet implemented");
            }

            public double grad(double x, double upper, double lower) {
                throw new RuntimeException("not yet implemented");
            }

            public double updateGradientLogDensity(double gradient, double value, double upper, double lower) {
                throw new RuntimeException("not yet implemented");
            }

            public double derivativeOfTransform(double x, double upper, double lower) {
                double scaledX = (x - lower) / (upper - lower);
                return (upper - lower) * scaledX * (1-scaledX);
            }

            public double secondDerivativeOfTransform(double x, double upper, double lower) {
                throw new RuntimeException("Not yet implemented");
            }

            public double derivativeOfInverseTransform(double y, double upper, double lower) {
                double inverse = 1 / (1 + Math.exp(-y));
                return (upper-lower) * inverse * (1 - inverse);
            }

            public double secondDerivativeOfInverseTransform(double x, double upper, double lower) {
                throw new RuntimeException("Not yet implemented");
            }

            public boolean isInteriorDomain(double x, double upper, double lower) {
                return x >= lower && x <= upper;
            }

            // Helper functions
            private double getScaledLogit(double x, double upper, double lower) {
                double u = (x - lower / (upper - lower));
                return Math.log(u / (1 - u));
            }

            private double getScaledSigmoid(double y, double upper, double lower) {
                double inverse = 1 / (1 + Math.exp(-y));
                return lower + (upper-lower) * inverse;
            }

        },
        NONE("none") {
            public double transform(double x, double upper, double lower) {
                return x;
            }

            public double inverse(double x, double upper, double lower) {
                return x;
            }

            public double grad(double x, double upper, double lower) {
                throw new RuntimeException("not yet implemented");
            }

            public double updateGradientLogDensity(double gradient, double value, double upper, double lower) {
                throw new RuntimeException("not yet implemented");
            }

            public double derivativeOfTransform(double x, double upper, double lower) {
                return 1.0;
            }

            public double secondDerivativeOfTransform(double x, double upper, double lower) {
                return 0.0;
            }

            public double derivativeOfInverseTransform(double x, double upper, double lower) {
                return 1.0;
            }

            public double secondDerivativeOfInverseTransform(double x, double upper, double lower) {
                return 0.0;
            }

            public boolean isInteriorDomain(double x, double upper, double lower) {
                return true;
            }

        };

        IncrementTransform(String transformType) {
            this.transformType = transformType;
        }

        public String getTransformType() {return transformType;}

        private String transformType;
        public abstract double transform(double x, double upper, double lower);
        public abstract double inverse(double x, double upper, double lower);
        public abstract double grad(double x, double upper, double lower);
        public abstract double updateGradientLogDensity(double gradient, double value, double upper, double lower);
        public abstract double derivativeOfTransform(double x, double upper, double lower);
        public abstract double secondDerivativeOfTransform(double x, double upper, double lower);
        public abstract double derivativeOfInverseTransform(double x, double upper, double lower);
        public abstract double secondDerivativeOfInverseTransform(double x, double upper, double lower);
        public abstract boolean isInteriorDomain(double x, double upper, double lower);

        public static FirstOrderFiniteDifferenceTransform.IncrementTransform factory(String match) {
            for (FirstOrderFiniteDifferenceTransform.IncrementTransform transform : FirstOrderFiniteDifferenceTransform.IncrementTransform.values()) {
                if (match.equalsIgnoreCase(transform.getTransformType())) {
                    return transform;
                }
            }
            return null;
        }
    }

    private final IncrementTransform incrementTransform;
    private final double upper;
    private final double lower;

    public FirstOrderFiniteDifferenceTransform(int dim, IncrementTransform incrementTransform, double upper, double lower) {
        super(dim);
        this.incrementTransform = incrementTransform;
        this.upper = upper;
        this.lower = lower;
    }

    @Override
    protected double[] transform(double[] values) {
        double[] incrementTransformedValues = values.clone();
        for (int i = 0; i < values.length; i++) {
            incrementTransformedValues[i] = incrementTransform.transform(values[i], upper, lower);
        }
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
        fx[0] = incrementTransform.inverse(values[0], upper, lower);
        double s = values[0];
        for (int i = 1; i < values.length; i++) {
            s += values[i];
            fx[i] = incrementTransform.inverse(s, upper, lower);
        }
        return fx;
    }

    @Override
    public double[] inverse(double[] values, int from, int to, double sum) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public boolean isInInteriorDomain(double[] values) {
        for (int i = 0; i < values.length; i++) {
            if ( !incrementTransform.isInteriorDomain(values[i], upper, lower)) {
                return false;
            }
        }
        return true;
    }

    public String getTransformName() {
        return "FiniteDifference";
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
//        updated[dim - 1] = gradient[dim - 1] * incrementTransform.derivativeOfInverseTransform(transformedValues[dim - 1], upper, lower);
//        for (int i = dim - 2; i > -1; i--) {
//            updated[i] = gradient[i] * incrementTransform.derivativeOfInverseTransform(transformedValues[i], upper, lower) + updated[i + 1];
//        }
//        double[] gradLogJacobian = getGradientLogJacobianInverse(values);
//        for (int i = dim - 1; i > -1; i--) {
//            updated[i] += gradLogJacobian[i];
//        }
//
//        return updated;
//
//    }

    @Override
    public double getLogJacobian(double[] values) {
        double logJacobian = 0.0;
        // Transform is lower triangular
        for (int i = 0; i < values.length; i++) {
            logJacobian += Math.log(incrementTransform.derivativeOfTransform(values[i], upper, lower));
        }
        return logJacobian;
    }

    @Override
    public double[] getGradientLogJacobianInverse(double[] values) {

//        // If x are the original values and y=f(x) are the transformed values,
//        // this is d/dy det(log(Jacobian(y->x))
//        // Call the incrementTransform g
//        // y[i] = g(x[i]) - g(x[i-1])
//        // x[i] = g^-1(sum_k=1^i y[i])
//        // The transform is triangular, so we only need the diagonals
        // jacobianDiagonal == diagonal of Jacobian of inverse transform
        double[] jacobianDiagonal = new double[values.length];
        double s = 0.0;
        double[] cumSum = new double[dim];
        for (int i = 0; i < values.length; i++) {
            s += values[i];
            cumSum[i] = s;
            jacobianDiagonal[i] = incrementTransform.derivativeOfInverseTransform(s, upper, lower);
        }
        double[] gradient = new double[values.length];
        double tmp = 0.0;
        for (int i = values.length - 1; i > -1; i--) {
            tmp += (1.0 / jacobianDiagonal[i]) * incrementTransform.secondDerivativeOfInverseTransform(cumSum[i], upper, lower);
            gradient[i] = tmp;
//            gradient[i] = 1.0 / jacobianDiagonal[i] * incrementTransform.secondDerivativeOfInverseTransform(s, upper, lower);
        }

//        double[] numGrad = getNumericalGradientLogJacobianInverse(values);
//        System.err.println("\n");
//        System.err.println("Values:              " + new dr.math.matrixAlgebra.Vector(values));
//        System.err.println("Transformed values:  " + new dr.math.matrixAlgebra.Vector(transform(values)));
//        System.err.println("Jacobian:            " + getLogJacobian(values));
//        System.err.println("Jacobian's Diagonal: " + new dr.math.matrixAlgebra.Vector(jacobianDiagonal));
//        System.err.println("Jacobian of inverse: " + getLogJacobianInverse(values));
//        System.err.println("Numerical gradient:  " + new dr.math.matrixAlgebra.Vector(numGrad));
//        System.err.println("Analytical gradient: " + new dr.math.matrixAlgebra.Vector(gradient));
//        System.err.println("\n");

        return gradient;
    }

    public double getLogJacobianInverse(double[] values) {
        double logJacobian = 0.0;
        double s = 0.0;
        // Inverse transform is lower triangular
        for (int i = 0; i < values.length; i++) {
            s += values[i];
            logJacobian += Math.log(incrementTransform.derivativeOfInverseTransform(s, upper, lower));
        }
        return logJacobian;
    }

    @Override
    // jacobian[j][i] = d x_i / d y_j
    public double[][] computeJacobianMatrixInverse(double[] values) {
        double[][] jacobian = new double[dim][dim];
        // Is this transposed wrong?
        double s = 0.0;
        for (int i = 0; i < dim; i++) {
            s += values[i];
            for (int j = 0; j < i + 1; j++) {
                jacobian[j][i] = incrementTransform.derivativeOfInverseTransform(s, upper, lower);
            }
        }
//        System.err.println(new dr.math.matrixAlgebra.Matrix(jacobian));
//        double[][] numJacob = computeNumericalJacobianInverse(values);
//        System.err.println(new dr.math.matrixAlgebra.Matrix(numJacob));
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
            double upper = xo.getAttribute("upper", Double.POSITIVE_INFINITY);
            double lower = xo.getAttribute("lower", Double.NEGATIVE_INFINITY);

            String ttype = (String) xo.getAttribute(INCREMENT_TRANSFORM, "none");
            FirstOrderFiniteDifferenceTransform.IncrementTransform transform = FirstOrderFiniteDifferenceTransform.IncrementTransform.factory(ttype);
            if (transform == FirstOrderFiniteDifferenceTransform.IncrementTransform.factory("logit") && !(Double.isFinite(lower) && Double.isFinite(upper))) {
                throw new RuntimeException("Logit transform on increments requires finite upper and lower bounds.");
            }

            int dim = xo.getIntegerAttribute(DIMENSION);

            return new FirstOrderFiniteDifferenceTransform(dim, transform, upper, lower);

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