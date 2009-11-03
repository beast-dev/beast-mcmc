/*
 * DistributionModelParser.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inferencexml;

import dr.inference.distribution.ExponentialDistributionModel;
import dr.inference.distribution.GammaDistributionModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 */
public abstract class DistributionModelParser extends AbstractXMLObjectParser {

    public static final String OFFSET = "offset";
    public static final String MEAN = "mean";
    public static final String SHAPE = "shape";
    public static final String SCALE = "scale";

    /**
     * @param parameters an array of the parsed parameters, in order of the getParameterNames() array.
     * @param offset     the parsed offset
     * @return a distribution model constructed from provided parameters and offset
     */
    abstract ParametricDistributionModel parseDistributionModel(Parameter[] parameters, double offset);

    /**
     * @return a list of xml element names for parameters of this distribution.
     */
    abstract String[] getParameterNames();

    abstract boolean allowOffset();

    public final Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double offset = xo.getAttribute(OFFSET, 0.0);

        String[] names = getParameterNames();
        Parameter[] parameters = new Parameter[names.length];
        for (int i = 0; i < names.length; i++) {
            parameters[i] = getParameter(xo, names[i]);
        }

        return parseDistributionModel(parameters, offset);
    }

    private Parameter getParameter(XMLObject xo, String parameterName) throws XMLParseException {
        final XMLObject cxo = xo.getChild(parameterName);
        return cxo.getChild(0) instanceof Parameter ?
                (Parameter) cxo.getChild(Parameter.class) : new Parameter.Default(cxo.getDoubleChild(0));
    }

    public XMLSyntaxRule[] getSyntaxRules() {

        String[] names = getParameterNames();

        XMLSyntaxRule[] rules = new XMLSyntaxRule[names.length + (allowOffset() ? 1 : 0)];
        for (int i = 0; i < names.length; i++) {
            rules[i] = new XORRule(
                    new ElementRule(names[i], Double.class),
                    new ElementRule(names[i], Parameter.class)
            );
        }
        if (allowOffset()) {
            rules[rules.length - 1] = AttributeRule.newDoubleRule(OFFSET, true);
        }

        return rules;
    }

    public static final XMLObjectParser EXPONENTIAL_DISTRIBUTION_PARSER = new DistributionModelParser() {

        public String getParserName() {
            return ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL;
        }

        ParametricDistributionModel parseDistributionModel(Parameter[] parameters, double offset) {
            return new ExponentialDistributionModel(parameters[0], offset);
        }

        public String[] getParameterNames() {
            return new String[]{MEAN};
        }

        public String getParserDescription() {
            return "A model of an exponential distribution.";
        }

        public boolean allowOffset() {
            return true;
        }

        public Class getReturnType() {
            return ExponentialDistributionModel.class;
        }
    };

    public static final XMLObjectParser GAMMA_DISTRIBUTION_PARSER = new DistributionModelParser() {
        public String getParserName() {
            return GammaDistributionModel.GAMMA_DISTRIBUTION_MODEL;
        }

        ParametricDistributionModel parseDistributionModel(Parameter[] parameters, double offset) {
            return new GammaDistributionModel(parameters[0], parameters[1]);
        }

        public String[] getParameterNames() {
            return new String[]{SHAPE, SCALE};
        }

        public String getParserDescription() {
            return "A model of a gamma distribution.";
        }

        public boolean allowOffset() {
            return false;
        }

        public Class getReturnType() {
            return GammaDistributionModel.class;
        }
    };

     public static final XMLObjectParser ONE_P_GAMMA_DISTRIBUTION_MODEL = new DistributionModelParser() {
        public String getParserName() {
            return GammaDistributionModel.ONE_P_GAMMA_DISTRIBUTION_MODEL;
        }

        ParametricDistributionModel parseDistributionModel(Parameter[] parameters, double offset) {
            return new GammaDistributionModel(parameters[0], null);
        }

        public String[] getParameterNames() {
            return new String[]{SHAPE};
        }

        public String getParserDescription() {
            return "A model of a one parameter gamma distribution.";
        }

        public boolean allowOffset() {
            return false;
        }

        public Class getReturnType() {
            return GammaDistributionModel.class;
        }
    };
}
