/*
 * FreeRateSimplexTransform.java
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

import dr.evomodel.substmodel.InfinitesimalRatesLogger;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.xml.*;
import org.apache.commons.math.stat.descriptive.moment.Mean;

import java.math.BigDecimal;
import java.math.RoundingMode;

public class FreeRateSimplexTransform extends Transform.MultivariateTransform {

    private Parameter weights;

    public FreeRateSimplexTransform(int dim, Parameter weights) {
        super(dim);
        this.weights = weights;
        this.outputDimension = dim + 1;
    }

    public FreeRateSimplexTransform(int dim) {
        super(dim);
        weights = new Parameter.Default(dim, (double) 1 /dim);
        this.outputDimension = dim + 1;
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
        return "freeRateSimplex";
    }


    @Override
    protected double[] transform(double[] values) {

        // This is a transformation of an n-dimensional vector to another n-dimensional vector but what comes _out_ is a
        // simplex of one greater dimension. The weights parameter has dimension n+1.

        BigDecimal[] bigOut = new BigDecimal[values.length + 1];

        BigDecimal denominator = BigDecimal.valueOf(0);

        for(int i=0; i<dim; i++){
            denominator = denominator.add(BigDecimal.valueOf(values[i]).multiply(BigDecimal.valueOf(weights.getParameterValue(i))));
        }

        denominator = denominator.add(BigDecimal.valueOf(weights.getParameterValue(dim)));

        for(int i=0; i<dim; i++){
            bigOut[i] = BigDecimal.valueOf(values[i]).divide(denominator, 12, RoundingMode.HALF_UP);
        }

        bigOut[dim] = BigDecimal.valueOf(1).divide(denominator, 12, RoundingMode.HALF_UP);
        double[] out = new double[dim+1];

        for(int i=0; i<=dim; i++){
            out[i] = bigOut[i].doubleValue();
        }

        // just testing the reverse works!

        BigDecimal[] bigIn = new BigDecimal[values.length];

        BigDecimal reverseDenominatorTemp = BigDecimal.valueOf(1);

        for(int i=0; i<dim; i++){
            reverseDenominatorTemp = reverseDenominatorTemp.subtract(bigOut[i].multiply(BigDecimal.valueOf(weights.getParameterValue(i))));
        }

        for(int i=0; i<dim; i++){
            bigIn[i] = bigOut[i].multiply(BigDecimal.valueOf(weights.getParameterValue(dim))).divide(reverseDenominatorTemp, 12, RoundingMode.HALF_UP);

        }

        for(int i=0; i<dim; i++){
            BigDecimal tmp = BigDecimal.valueOf(values[i]);
        }

        return(out);

    }

    @Override
    protected double[] inverse(double[] values) {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    protected double getLogJacobian(double[] values) {
        Matrix partialsMatrix = new Matrix(dim, dim);

        double logSqrtDenominator;

        // reminder: dim is one less than the dimension of the simplex

        double tempSum = 0;

        for(int i=0; i<dim; i++){
            tempSum += weights.getParameterValue(i)*values[i];
        }
        logSqrtDenominator =  Math.log1p(-tempSum);

        assert logSqrtDenominator > -Double.NEGATIVE_INFINITY;

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                if (j == i) {
                    partialsMatrix.set(i, j, ((1-tempSum) + values[i]*weights.getParameterValue(i)));
                } else {
                    partialsMatrix.set(i, j, (values[i]*weights.getParameterValue(j)));
                }
            }
        }

        double logJacobian = 0;
        try {

            logJacobian =  dim*Math.log(weights.getParameterValue(dim)) + Math.log(Math.abs(partialsMatrix.determinant())) - dim*2*logSqrtDenominator;
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
        public static final String FREERATE_SIMPLEX_TRANSFORM = "freeRateSimplexTransform";
        @Override
        public Class getReturnType() {
            return FreeRateSimplexTransform.class;
        }

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            if(xo.hasChildNamed(WEIGHTS)){
                Parameter weights = (Parameter) xo.getElementFirstChild(WEIGHTS);
                return new FreeRateSimplexTransform(weights.getDimension()-1, weights);
            } else if(xo.hasAttribute(DIMENSION)) {
                int dimension = xo.getIntegerAttribute(DIMENSION);
                return new FreeRateSimplexTransform(dimension);
            } else {
                throw new XMLParseException("FreeRateSimplexTransform must have either a dimension attribute or" +
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
            return FREERATE_SIMPLEX_TRANSFORM;
        }


    };

}
