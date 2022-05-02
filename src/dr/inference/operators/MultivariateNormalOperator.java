/*
 * MultivariateNormalOperator.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.operators;

import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.*;
import dr.xml.*;


/**
 * @author Marc Suchard
 */
public class MultivariateNormalOperator extends AbstractAdaptableOperator {

    public static final String MVN_OPERATOR = "mvnOperator";
    public static final String SCALE_FACTOR = "scaleFactor";
    public static final String VARIANCE_MATRIX = "varMatrix";
    public static final String FORM_XTX = "formXtXInverse";

    private double scaleFactor;
    private final Parameter parameter;
    private final int dim;

    private double[][] cholesky;

    public MultivariateNormalOperator(Parameter parameter, double scaleFactor, double[][] inMatrix, double weight,
                                      AdaptationMode mode, boolean isVarianceMatrix) {

        super(mode);
        this.scaleFactor = scaleFactor;
        this.parameter = parameter;
        setWeight(weight);
        dim = parameter.getDimension();

        SingularValueDecomposition svd = new SingularValueDecomposition(new DenseDoubleMatrix2D(inMatrix));
        if (inMatrix[0].length != svd.rank()) {
            throw new RuntimeException("Variance matrix in mvnOperator is not of full rank");
        }

        final double[][] matrix;
        if (isVarianceMatrix) {
            matrix = inMatrix;
        } else {
            matrix = formXtXInverse(inMatrix);
        }

        try {
            cholesky = (new CholeskyDecomposition(matrix)).getL();
        } catch (IllegalDimension illegalDimension) {
            throw new RuntimeException("Unable to decompose matrix in mvnOperator");
        }
    }

    public MultivariateNormalOperator(Parameter parameter, double scaleFactor,
                                      MatrixParameter varMatrix, double weight, AdaptationMode mode, boolean isVariance) {
        this(parameter, scaleFactor, varMatrix.getParameterAsMatrix(), weight, mode, isVariance);
    }

    private double[][] formXtXInverse(double[][] X) {
        int N = X.length;
        int P = X[0].length;

        double[][] matrix = new double[P][P];
        for (int i = 0; i < P; i++) {
            for (int j = i; j < P; j++) {
                double total = 0.0;
                for (int k = 0; k < N; k++) {
                    total += X[k][i] * X[k][j];
                }
                matrix[j][i] = matrix[i][j] = total;
            }
        }

        // Take inverse
        double[][] inverse = new SymmetricMatrix(matrix).inverse().toComponents();

        // Force symmetric
        for (int i = 0; i < inverse.length; ++i) {
            for (int j = i; j < inverse[i].length; ++j) {
                inverse[j][i] = inverse[i][j] = (inverse[j][i] + inverse[i][j]) / 2;
            }
        }

        return inverse;
    }

    public double doOperation() {

        double[] x = parameter.getParameterValues();
        double[] epsilon = new double[dim];
        //double[] y = new double[dim];
        for (int i = 0; i < dim; i++)
            epsilon[i] = scaleFactor * MathUtils.nextGaussian();

        for (int i = 0; i < dim; i++) {
            for (int j = i; j < dim; j++) {
                x[i] += cholesky[j][i] * epsilon[j];
                // caution: decomposition returns lower triangular
            }
            parameter.setParameterValueQuietly(i, x[i]);
        }
        parameter.fireParameterChangedEvent();

        return 0;
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return parameter.getParameterName();
    }

    @Override
    protected double getAdaptableParameterValue() {
        return Math.log(scaleFactor);
    }

    public void setAdaptableParameterValue(double value) {
        scaleFactor = Math.exp(value);
    }

    public double getRawParameter() {
        return scaleFactor;
    }

    public double getScaleFactor() {
        return scaleFactor;
    }

    public String getAdaptableParameterName() {
        return "scaleFactor";
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MVN_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            AdaptationMode mode = AdaptationMode.parseMode(xo);

            double weight = xo.getDoubleAttribute(WEIGHT);
            double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);

            if (scaleFactor <= 0.0) {
                throw new XMLParseException("scaleFactor must be greater than 0.0");
            }

            Parameter parameter = (Parameter) xo.getChild(Parameter.class);

            boolean formXtXInverse = xo.getAttribute(FORM_XTX, false);

            XMLObject cxo = xo.getChild(VARIANCE_MATRIX);
            MatrixParameter varMatrix = (MatrixParameter) cxo.getChild(MatrixParameter.class);

            // Make sure varMatrix is square and dim(varMatrix) = dim(parameter)

            if (!formXtXInverse) {
                if (varMatrix.getColumnDimension() != varMatrix.getRowDimension())
                    throw new XMLParseException("The variance matrix is not square");
            }

            if (varMatrix.getColumnDimension() != parameter.getDimension())
                throw new XMLParseException("The parameter and variance matrix have differing dimensions");

            return new MultivariateNormalOperator(parameter, scaleFactor, varMatrix, weight, mode, !formXtXInverse);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a multivariate normal random walk operator on a given parameter.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(SCALE_FACTOR),
                AttributeRule.newDoubleRule(WEIGHT),
                AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
                AttributeRule.newBooleanRule(FORM_XTX, true),
                new ElementRule(Parameter.class),
                new ElementRule(VARIANCE_MATRIX,
                        new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)})

        };

    };
}
