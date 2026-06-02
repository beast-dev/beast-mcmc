/*
 * GammaDistributionModelParser.java
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

package dr.inferencexml.distribution;

import dr.inference.distribution.GammaDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

public class GammaDistributionModelParser extends AbstractXMLObjectParser {

    public static final String MEAN = "mean";
    public static final String SHAPE = "shape";
    public static final String SCALE = "scale";
    public static final String RATE = "rate";
    public static final String OFFSET = "offset";

    public String getParserName() {
        return GammaDistributionModel.GAMMA_DISTRIBUTION_MODEL;
    }

    private Parameter getParameterOrValue(String name, XMLObject xo) throws XMLParseException {
        Parameter parameter;
        if (xo.hasChildNamed(name)) {
            XMLObject cxo = xo.getChild(name);

            parameter = (Parameter)cxo.getChild(Parameter.class);
            if (parameter == null) {
                if (cxo.getChildCount() < 1) {
                    throw new XMLParseException("Distribution parameter, " + name + ", is missing a value or parameter element");
                }
                try {
                    double value = cxo.getDoubleChild(0);
                    parameter = new Parameter.Default(value);
                } catch (XMLParseException xpe) {
                    throw new XMLParseException("Distribution parameter, " + name + ", has bad value: " + xpe.getMessage());
                }
            }

            return parameter;
        } else {
            return null;
        }

    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double offset = xo.getAttribute(OFFSET, 0.0);

        Parameter shapeParameter = getParameterOrValue(SHAPE, xo);

        Parameter parameter2;
        GammaDistributionModel.GammaParameterizationType parameterization;

        if (xo.hasChildNamed(SCALE)) {
            parameter2 = getParameterOrValue(SCALE, xo);
            parameterization = GammaDistributionModel.GammaParameterizationType.ShapeScale;
        } else if (xo.hasChildNamed(RATE)) {
            parameter2 = getParameterOrValue(RATE, xo);
            parameterization = GammaDistributionModel.GammaParameterizationType.ShapeRate;
        } else if (xo.hasChildNamed(MEAN)) {
            parameter2 = getParameterOrValue(MEAN, xo);
            parameterization = GammaDistributionModel.GammaParameterizationType.ShapeMean;
        } else {
            parameter2 = null;
            parameterization = GammaDistributionModel.GammaParameterizationType.OneParameter;
        }

        return new GammaDistributionModel(parameterization, shapeParameter, parameter2, offset);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(SHAPE,  new XMLSyntaxRule[]{new ElementRule(Parameter.class, true)}, "Shape parameter"),
            new ElementRule(SCALE,  new XMLSyntaxRule[]{new ElementRule(Parameter.class, true)}, "Scale parameter", true),
            new ElementRule(RATE,  new XMLSyntaxRule[]{new ElementRule(Parameter.class, true)}, "Rate parameter", true),
            new ElementRule(MEAN,  new XMLSyntaxRule[]{new ElementRule(Parameter.class, true)}, "Mean parameter", true),
            AttributeRule.newDoubleRule(OFFSET, true)
    };

    public String getParserDescription() {
        return "The gamma probability distribution.";
    }

    public Class getReturnType() {
        return GammaDistributionModel.class;
    }
}
