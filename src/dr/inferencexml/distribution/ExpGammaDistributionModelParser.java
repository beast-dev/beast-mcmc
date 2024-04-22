/*
 * ExpGammaDistributionModelParser.java
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

import dr.inference.distribution.ExpGammaDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

public class ExpGammaDistributionModelParser extends AbstractXMLObjectParser {

    public static final String SHAPE = "shape";
    public static final String SCALE = "scale";

    public String getParserName() {
        return ExpGammaDistributionModel.EXP_GAMMA_DISTRIBUTION_MODEL;
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

        Parameter shapeParameter = getParameterOrValue(SHAPE, xo);
        Parameter scaleParameter = getParameterOrValue(SCALE, xo);

        return new ExpGammaDistributionModel(shapeParameter, scaleParameter);
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
    };

    public String getParserDescription() {
        return "The probability distribution for the natural logarithm of a gamma-distributed random variable.";
    }

    public Class getReturnType() {
        return ExpGammaDistributionModel.class;
    }
}
