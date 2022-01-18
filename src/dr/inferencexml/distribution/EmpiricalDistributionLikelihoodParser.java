/*
 * EmpiricalDistributionLikelihoodParser.java
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

package dr.inferencexml.distribution;

import dr.inference.distribution.EmpiricalDistributionData;
import dr.inference.distribution.EmpiricalDistributionLikelihood;
import dr.inference.distribution.SplineInterpolatedLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Marc A. Suchard
 */
public class EmpiricalDistributionLikelihoodParser extends AbstractXMLObjectParser {

    public static final String FILE_NAME = "fileName";
    public static final String DATA = "data";
    public static final String FROM = "from";
    public static final String TO = "to";
    private static final String SPLINE_INTERPOLATION = "splineInterpolation";
    private static final String DEGREE = "degree";
    private static final String INVERSE = "inverse";
    private static final String READ_BY_COLUMN = "readByColumn";
    public static final String OFFSET="offset";
    public static final String LOWER = "lower";
    public static final String UPPER = "upper";

    private static final String FILE_INFORMATION = "fileInformation";
    private static final String RAW_VALUES = "grid";
    private static final String LIKELIHOOD = "logLikelihood";
    private static final String VALUES = "value";
    private static final String DENSITY_IN_LOG_SPACE = "densityInLogSpace";

    public String getParserName() {
        return EmpiricalDistributionLikelihood.EMPIRICAL_DISTRIBUTION_LIKELIHOOD;
    }

    private Parameter getRawValues(String name, XMLObject xo) throws XMLParseException {
        return getRawValues(xo.getChild(name));
    }

    private Parameter getRawValues(XMLObject cxo) throws XMLParseException {
        Parameter parameter;
        if (cxo.getChild(0) instanceof Parameter) {
            parameter = (Parameter) cxo.getChild(Parameter.class);
        } else {
            parameter = new Parameter.Default(cxo.getDoubleArrayChild(0));
        }
        return parameter;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean splineInterpolation = xo.getAttribute(SPLINE_INTERPOLATION, false);
        int degree = xo.getAttribute(DEGREE, 3); // Default is cubic-spline
        boolean inverse = xo.getAttribute(INVERSE, false);

        EmpiricalDistributionLikelihood likelihood;
        if (xo.hasChildNamed(FILE_INFORMATION)) {

            String fileName = xo.getStringAttribute(FILE_NAME);
            boolean byColumn = xo.getAttribute(READ_BY_COLUMN, true);

            if (splineInterpolation) {
                if (degree < 1)
                    throw new XMLParseException("Spline degree must be greater than zero!");
                likelihood = new SplineInterpolatedLikelihood(fileName, degree, inverse, byColumn);
            } else
                //likelihood = new EmpiricalDistributionLikelihood(fileName,inverse,byColumn);
                throw new XMLParseException("Only spline-interpolated empirical distributions are currently support");

        } else {

            List<EmpiricalDistributionData> dataList = new ArrayList<>();

            for (int i = 0; i < xo.getChildCount(); ++i) {
                XMLObject cxo = (XMLObject) xo.getChild(i);
                if (cxo.getName().equals(RAW_VALUES)) {
                    Parameter thetaParameter = getRawValues(VALUES, cxo);
                    Parameter likelihoodParameter = getRawValues(LIKELIHOOD, cxo);

                    if (likelihoodParameter.getDimension() != thetaParameter.getDimension()) {
                        throw new XMLParseException("Unequal grid lengths");
                    }

                    boolean densityInLogSpace = cxo.getAttribute(DENSITY_IN_LOG_SPACE, true);

                    dataList.add(new EmpiricalDistributionData(
                            thetaParameter.getParameterValues(), likelihoodParameter.getParameterValues(),
                            densityInLogSpace));
                }
            }

            likelihood = new SplineInterpolatedLikelihood(dataList, degree, inverse);
        }

        XMLObject cxo1 = xo.getChild(DATA);
        final int from = cxo1.getAttribute(FROM, -1);
        int to = cxo1.getAttribute(TO, -1);
        if (from >= 0 || to >= 0) {
            if (to < 0) {
                to = Integer.MAX_VALUE;
            }
            if (!(from >= 0 && to >= 0 && from < to)) {
                throw new XMLParseException("ill formed from-to");
            }
            likelihood.setRange(from, to);
        }

        for (int j = 0; j < cxo1.getChildCount(); j++) {
            if (cxo1.getChild(j) instanceof Statistic) {

                likelihood.addData((Statistic) cxo1.getChild(j));
            } else {
                throw new XMLParseException("illegal element in " + cxo1.getName() + " element");
            }
        }

        if (likelihood.getDistributionDimension() != 1 &&
                (likelihood.getDistributionDimension() != likelihood.getDimension())) {
            throw new XMLParseException("Data dimension != distribution dimension");
        }

        double offset = cxo1.getAttribute(OFFSET,0); 
        likelihood.setOffset(offset);

        if (cxo1.hasAttribute(LOWER) || cxo1.hasAttribute(UPPER)) {
            likelihood.setBounds(
                    cxo1.getAttribute(LOWER, Double.NEGATIVE_INFINITY),
                    cxo1.getAttribute(UPPER, Double.POSITIVE_INFINITY)
            );
        }

        return likelihood;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {

            AttributeRule.newBooleanRule(SPLINE_INTERPOLATION, true),
            AttributeRule.newIntegerRule(DEGREE, true),
            AttributeRule.newBooleanRule(INVERSE, true),
            new XORRule(
                    new ElementRule(FILE_INFORMATION, new XMLSyntaxRule[]{
                            AttributeRule.newBooleanRule(READ_BY_COLUMN, true),
                            AttributeRule.newStringRule(FILE_NAME),
                    }),
                    new ElementRule(RAW_VALUES, new XMLSyntaxRule[]{
                            AttributeRule.newBooleanRule(DENSITY_IN_LOG_SPACE, true),
                            new ElementRule(LIKELIHOOD,
                                    new XMLSyntaxRule[]{
                                            new XORRule(
                                                    new ElementRule(Parameter.class),
                                                    new ElementRule(Double.class) // TODO Fix for array
                                            )}),
                            new ElementRule(VALUES,
                                    new XMLSyntaxRule[]{
                                            new XORRule(
                                                    new ElementRule(Parameter.class),
                                                    new ElementRule(Double.class) // TODO Fix for array
                                            )}
                            ),

                    }, 1, Integer.MAX_VALUE)),
            new ElementRule(DATA, new XMLSyntaxRule[]{
                    AttributeRule.newIntegerRule(FROM, true),
                    AttributeRule.newIntegerRule(TO, true),
                    AttributeRule.newDoubleRule(OFFSET, true),
                    AttributeRule.newDoubleRule(LOWER, true),
                    AttributeRule.newDoubleRule(UPPER, true),
                    new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
            })
    };

    public String getParserDescription() {
        return "Calculates the likelihood of some data given some empirically-generated distribution.";
    }

    public Class getReturnType() {
        return Likelihood.class;
    }
}