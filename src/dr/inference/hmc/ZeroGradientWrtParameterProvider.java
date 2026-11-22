/*
 * ZeroGradientWrtParameterProvider.java
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

package dr.inference.hmc;

import dr.inference.model.ConstantLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Supplies a zero gradient for a parameter block while preserving the block's
 * dimension in a compound gradient.
 *
 * @author Filippo Monti
 */
public class ZeroGradientWrtParameterProvider implements GradientWrtParameterProvider {

    private static final String ZERO_GRADIENT = "zeroGradient";

    private final Parameter parameter;
    private final Likelihood likelihood;

    public ZeroGradientWrtParameterProvider(Parameter parameter) {
        this.parameter = parameter;
        this.likelihood = new ConstantLikelihood(0.0);
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return parameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        return new double[parameter.getDimension()];
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public String getParserName() {
            return ZERO_GRADIENT;
        }

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            Parameter parameter = (Parameter) xo.getChild(Parameter.class);
            return new ZeroGradientWrtParameterProvider(parameter);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        @Override
        public String getParserDescription() {
            return "Creates a zero gradient provider for the supplied parameter.";
        }

        @Override
        public Class getReturnType() {
            return ZeroGradientWrtParameterProvider.class;
        }

        private final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newStringRule("id", true),
                new ElementRule(Parameter.class),
        };
    };
}
