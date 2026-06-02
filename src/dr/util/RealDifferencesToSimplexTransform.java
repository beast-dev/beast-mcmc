/*
 * RealDifferencesToSimplexTransform.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
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

import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.xml.*;

public class RealDifferencesToSimplexTransform extends Transform.MultivariateTransform {

    private Parameter weights;

    public RealDifferencesToSimplexTransform(int dim, Parameter weights) {
        super(dim);
        this.weights = weights;
    }

    public RealDifferencesToSimplexTransform(int dim) {
        super(dim);
        weights = new Parameter.Default(dim, (double) 1 /dim);
    }

    @Override
    public double[] inverse(double[] values, int from, int to, double sum) {
        throw new RuntimeException("not implemented");
    }

    @Override
    public double[] gradient(double[] values, int from, int to) {
        throw new RuntimeException("not implemented");
    }

    @Override
    public double[] gradientInverse(double[] values, int from, int to) {
        throw new RuntimeException("not implemented");
    }

    @Override
    public String getTransformName() {
        return "realDifferencesToSimplex";
    }

    @Override
    protected double[] transform(double[] values) {

        // This is a transformation of an n-dimensional vector to another n-dimensional vector but what comes _out_ is a
        // simplex of one greater dimension. The weights parameter has dimension n+1.

        double[] out = new double[values.length + 1];

        double denominator = 0;

        for(int i=0; i<=dim; i++){
            double innerSum = 0;
            for(int j=0; j<=i; j++) {
                if(j<dim){
                    innerSum += values[j];
                } else {
                    innerSum += 1;

                }
            }
            denominator += innerSum*weights.getParameterValue(i);
        }

        for(int i=0; i<dim; i++){
            if(i==0){
                out[i] = values[i]/denominator;
            } else if(i<dim) {
                out[i] = (out[i-1]*denominator + values[i])/denominator;
            }
        }

        double subtotal = 0;

        for(int i=0; i<dim; i++){
            subtotal += out[i]*weights.getParameterValue(i);
        }

        out[dim] = (1-subtotal)/weights.getParameterValue(dim);

        return(out);

    }

    @Override
    protected double[] inverse(double[] values) {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    protected double getLogJacobian(double[] values) {
        Matrix partialsMatrix = new Matrix(dim, dim);

        double sqrtDenominator = 0;
        for(int i=0; i<=dim; i++){
            double innerSum = 0;
            for(int j=0; j<=i; j++) {
                if(j<dim){
                    innerSum += values[j];
                } else {
                    innerSum += 1;
                }
            }
            sqrtDenominator += innerSum*weights.getParameterValue(i);
        }

        double denominator = sqrtDenominator * sqrtDenominator;

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                double weightSum = 0;
                for (int d = j; d <= dim; d++) {
                    weightSum += weights.getParameterValue(d);
                }
                double valueSum = 0;
                for (int a = 0; a < i + 1; a++) {
                    valueSum += values[a];
                }
                if (j <= i) {
                    partialsMatrix.set(i, j, (sqrtDenominator - weightSum * valueSum) / denominator);
                } else {
                    partialsMatrix.set(i, j, (-weightSum * valueSum) / denominator);
                }
            }
        }

        double logJacobian = 0;
        try {
            logJacobian = Math.log(Math.abs(partialsMatrix.determinant()));
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }
        return logJacobian;
    }

    @Override
    protected double[] getGradientLogJacobianInverse(double[] values) {
        throw new RuntimeException("not implemented");
    }

    @Override
    public double[][] computeJacobianMatrixInverse(double[] values) {
        throw new RuntimeException("not implemented");
    }

    @Override
    protected boolean isInInteriorDomain(double[] values) {
        throw new RuntimeException("not implemented");
    }

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public static final String WEIGHTS = "weights";
        public static final String DIMENSION = "dimension";
        public static final String REAL_DIFFERENCES_TO_SIMPLEX_TRANSFORM = "realDifferencesToSimplexTransform";
        @Override
        public Class getReturnType() {
            return RealDifferencesToSimplexTransform.class;
        }

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            if(xo.hasChildNamed(WEIGHTS)){
                Parameter weights = (Parameter) xo.getElementFirstChild(WEIGHTS);
                return new RealDifferencesToSimplexTransform(weights.getDimension()-1, weights);
            } else if(xo.hasAttribute(DIMENSION)) {
                int dimension = xo.getIntegerAttribute(DIMENSION);
                return new RealDifferencesToSimplexTransform(dimension);
            } else {
                throw new XMLParseException("RealDifferencesToSimplex must have either a dimension attribute or" +
                        "a set of weights");
            }
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new OrRule(
                    AttributeRule.newDoubleRule(DIMENSION),
                    new ElementRule(WEIGHTS, Parameter.class))
            };
        }

        @Override
        public String getParserDescription() {
            return "Transform from a (n-1)-dimensional vector of positive real numbers to a n-dimensional weighted " +
                    "simplex";
        }

        @Override
        public String getParserName() {
            return REAL_DIFFERENCES_TO_SIMPLEX_TRANSFORM;
        }


    };

}
