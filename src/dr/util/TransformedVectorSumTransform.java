/*
 * TransformedVectorSumTransform.java
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
import dr.inference.model.TransformedMultivariateParameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class TransformedVectorSumTransform extends Transform.MultivariateTransform {

    private static final String NAME = "transformedVectorSumTransform";
    private static final String PARSER_NAME2 = "vectorScanTransformedParameter";
    private final Transform incrementTransform;

    public TransformedVectorSumTransform(int dim, Transform incrementTransform) {
        super(dim);
        this.incrementTransform = incrementTransform;
    }

    @Override
    public double[] inverse(double[] values, int from, int to, double sum) {
        throw new RuntimeException("Not yet implemented.");
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
    public String getTransformName() {
        return NAME;
    }

    @Override
    protected double[] updateGradientLogDensity(double[] transformedGradient, double[] untransformedValues) {

        final int dim = untransformedValues.length;

        double[] transformedValues = transform(untransformedValues); // TODO This seem unnecessary; maybe change interface to pass these values?
        double[] untransformedGradient = new double[dim];

        untransformedGradient[dim - 1] = transformedGradient[dim - 1] *
                incrementTransform.gradient(transformedValues[dim - 1]);
        for (int i = dim - 2; i >= 0; --i) {
            untransformedGradient[i] = transformedGradient[i] *
                    incrementTransform.gradient(transformedValues[i]) + untransformedGradient[i + 1];
        }

        return untransformedGradient;
    }

    @Override
    protected double[] transform(double[] values) {
        double[] fx = new double[values.length];
        fx[0] = values[0];
        for (int i = 1; i < values.length; i++) {
            fx[i] = fx[i-1] + values[i];
        }
        return incrementTransform.inverse(fx, 0, values.length);
    }

    @Override
    protected double[] inverse(double[] values) {
        values = incrementTransform.transform(values,0, values.length);
        double[] increments = new double[values.length];
        increments[0] = values[0];
        for (int i = 1; i < values.length; i++) {
            increments[i] = values[i] - values[i - 1];
        }
        return increments;
    }

    @Override
    protected double getLogJacobian(double[] values) {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    protected double[] getGradientLogJacobianInverse(double[] values) {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    public double[][] computeJacobianMatrixInverse(double[] values) {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    protected boolean isInInteriorDomain(double[] values) {
        throw new RuntimeException("Not yet implemented.");
    }

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        public final static String INCREMENT_TRANSFORM = "incrementTransformType";

        @Override
        public TransformedMultivariateParameter parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter param = (Parameter) xo.getChild(Parameter.class);

            int dim = param.getDimension();

            double upper = Double.POSITIVE_INFINITY;
            double lower = Double.NEGATIVE_INFINITY;
            if( xo.hasAttribute("upper") && xo.hasAttribute("lower")) {
                upper = xo.getDoubleAttribute("upper");
                lower = xo.getDoubleAttribute("lower");
            }

            Transform incrementTransform;
            String transformType = (String) xo.getAttribute(INCREMENT_TRANSFORM);
            if (transformType.equalsIgnoreCase("log")) {
                incrementTransform = Transform.LOG;
            } else if (transformType.equalsIgnoreCase("logit")) {
                incrementTransform = new Transform.ScaledLogitTransform(upper, lower);
            } else if (transformType.equalsIgnoreCase("none")) {
                incrementTransform = new Transform.NoTransform();
            } else {
                throw new RuntimeException("Invalid transform type");
            }

            TransformedVectorSumTransform transform = new TransformedVectorSumTransform(dim, incrementTransform);

            return new TransformedMultivariateParameter(param, transform);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[0];
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return TransformedMultivariateParameter.class;
        }

        @Override
        public String getParserName() {
            return NAME;
        }

        @Override
        public String[] getParserNames() { return new String[]{NAME, PARSER_NAME2}; }
    };
}