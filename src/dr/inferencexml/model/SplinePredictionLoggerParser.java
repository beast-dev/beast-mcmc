/*
 * MatrixMatrixProductParser.java
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

package dr.inferencexml.model;

import dr.inference.model.Parameter;
import dr.math.distributions.SplineBasisMatrix;
import dr.math.distributions.SplinePredictionLogger;
import dr.xml.*;

import java.util.List;

/**
 * @author Marc Suchard
 * @author Pratyusa Datta
 */

public class SplinePredictionLoggerParser extends AbstractXMLObjectParser {

    private final static String SPLINE_PREDICTION = "bSplinePrediction";
    private final static String PREDICTION_POINTS = "predictionPoints";
    private final static String MIN = "min";
    private final static String MAX = "max";
    private final static String LENGTH = "length";
    private final static String FORMAT = "format";
    private final static String INTERCEPT = "intercept";

    private final XMLSyntaxRule[] rules = {
            new ElementRule(SplineBasisMatrix.class),
            new ElementRule(Parameter.class, 2, 2),
            AttributeRule.newStringRule(FORMAT, true),
            new XORRule(
                    new ElementRule(PREDICTION_POINTS, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
                    new AndRule(new XMLSyntaxRule[]{
                            new ElementRule(MIN, new XMLSyntaxRule[]{
                                    new XORRule(
                                            new ElementRule(Parameter.class),
                                            new ElementRule(Double.class)
                                    )}
                            ),
                            new ElementRule(MAX, new XMLSyntaxRule[]{
                                    new XORRule(
                                            new ElementRule(Parameter.class),
                                            new ElementRule(Double.class)
                                    )}
                            ),
                            new ElementRule(LENGTH, new XMLSyntaxRule[]{
                                    new ElementRule(Integer.class)
                            }),
                    })
            ),
            new ElementRule(INTERCEPT, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class),
            }, true),
    };

    private double getMin(Parameter p) {
        double min = Double.POSITIVE_INFINITY;
        for (int i = 0; i < p.getDimension(); ++i) {
            double x = p.getParameterValue(i);
            if (x < min) {
                min = x;
            }
        }
        return min;
    }

    private double getMax(Parameter p) {
        double max = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < p.getDimension(); ++i) {
            double x = p.getParameterValue(i);
            if (x > max) {
                max = x;
            }
        }
        return max;
    }

    @Override
    public Object parseXMLObject (XMLObject xo) throws XMLParseException {

        SplineBasisMatrix splineBasis;
        Parameter coefficients;

        List<Parameter> parameterList = xo.getAllChildren(Parameter.class);
        if (parameterList.get(0) instanceof SplineBasisMatrix) {
            splineBasis = (SplineBasisMatrix) parameterList.get(0);
            coefficients = parameterList.get(1);
        } else {
            coefficients = parameterList.get(0);
            splineBasis = (SplineBasisMatrix) parameterList.get(1);
        }

        Parameter predictionPoints;
        if (xo.hasChildNamed(PREDICTION_POINTS)) {
            predictionPoints = (Parameter) xo.getElementFirstChild(PREDICTION_POINTS);
        } else {

            XMLObject minXo = xo.getChild(MIN);
            double min = minXo.getChild(0) instanceof Parameter ?
                    getMin((Parameter) minXo.getChild(Parameter.class)) :
                    minXo.getDoubleChild(0);

            XMLObject maxXo = xo.getChild(MAX);
            double max = minXo.getChild(0) instanceof Parameter ?
                    getMax((Parameter) maxXo.getChild(Parameter.class)) :
                    maxXo.getDoubleChild(0);

            int length = xo.getChild(LENGTH).getIntegerChild(0);

            double[] values = new  double[length];
            for (int i = 0; i < length - 1; ++i) {
                values[i] = min + (max - min) / (length + 1) * i;
            }
            values[length - 1] = max;

            predictionPoints = new Parameter.Default(values);
        }

        Parameter intercept = null;
        if (xo.hasChildNamed(INTERCEPT)) {
            intercept = (Parameter) xo.getElementFirstChild(INTERCEPT);
        }

        String name;
        if (xo.hasId()) {
            name = xo.getId();
        } else {
            name = splineBasis.getId() + ".prediction";
        }

        String format = null;
        if (xo.hasAttribute(FORMAT)) {
            format = xo.getStringAttribute(FORMAT);
        }

        return new SplinePredictionLogger(name, splineBasis, coefficients, predictionPoints, intercept, format);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules () {
        return rules;
    }

    @Override
    public String getParserDescription () {
        return "Gets Latent Factor Model to return data with residuals computed";
    }

    @Override
    public Class getReturnType () {
        return SplinePredictionLogger.class;
    }

    @Override
    public String getParserName () {
        return SPLINE_PREDICTION;
    }
}
