/*
 * PrecisionMatrixGibbsOperator.java
 *
 * Copyright (C) 2002-2007 Alexei Drummond and Andrew Rambaut
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

import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.xml.*;


/**
 * @author Marc Suchard
 */
public class MultivariateNormalOperator extends AbstractCoercableOperator {

    public static final String MVN_OPERATOR = "mvnOperator";
    public static final String SCALE_FACTOR = "scaleFactor";
    public static final String VARIANCE_MATRIX = "varMatrix";

    private double scaleFactor;
    private final Parameter parameter;
    private final int dim;

    private double[][] cholesky;

    public MultivariateNormalOperator(Parameter parameter, double scaleFactor,
                                      MatrixParameter varMatrix, double weight, CoercionMode mode) {
        super(mode);
        this.scaleFactor = scaleFactor;
        this.parameter = parameter;
        setWeight(weight);
        dim = parameter.getDimension();
        cholesky = new double[dim][dim];
        for (int i = 0; i < dim; i++) {
            for (int j = i; j < dim; j++)
                cholesky[i][j] = cholesky[j][i] = varMatrix.getParameterValue(i, j);
        }

        try {
            cholesky = (new CholeskyDecomposition(cholesky)).getL();
        } catch (IllegalDimension illegalDimension) {
        }

    }

    public double doOperation() throws OperatorFailedException {

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
            parameter.setParameterValue(i, x[i]);
//            System.out.println(i+" : "+x[i]);
        }
//                    System.exit(-1);
        return 0;
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return parameter.getParameterName();
    }

    public double getCoercableParameter() {
        return Math.log(scaleFactor);
    }

    public void setCoercableParameter(double value) {
        scaleFactor = Math.exp(value);
    }

    public double getRawParameter() {
        return scaleFactor;
    }

    public double getScaleFactor() {
        return scaleFactor;
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    public double getMinimumAcceptanceLevel() {
        return 0.1;
    }

    public double getMaximumAcceptanceLevel() {
        return 0.4;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return 0.20;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return 0.30;
    }

    public final String getPerformanceSuggestion() {

        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();
        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
        double sf = OperatorUtils.optimizeWindowSize(scaleFactor, prob, targetProb);
        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return MVN_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            CoercionMode mode = CoercionMode.parseMode(xo);

            double weight = xo.getDoubleAttribute(WEIGHT);
            double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);

            if (scaleFactor <= 0.0) {
                throw new XMLParseException("scaleFactor must be greater than 0.0");
            }

            Parameter parameter = (Parameter) xo.getChild(Parameter.class);


            XMLObject cxo = (XMLObject) xo.getChild(VARIANCE_MATRIX);
            MatrixParameter varMatrix = (MatrixParameter) cxo.getChild(MatrixParameter.class);

            // Make sure varMatrix is square and dim(varMatrix) = dim(parameter)

            if (varMatrix.getColumnDimension() != varMatrix.getRowDimension())
                throw new XMLParseException("The variance matrix is not square");

            if (varMatrix.getColumnDimension() != parameter.getDimension())
                throw new XMLParseException("The parameter and variance matrix have differing dimensions");

            return new MultivariateNormalOperator(parameter, scaleFactor, varMatrix, weight, mode);
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
                new ElementRule(Parameter.class),
                new ElementRule(VARIANCE_MATRIX,
                        new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)})

        };

    };
}
