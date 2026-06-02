/*
 * UniformDistributionModelParser.java
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

import dr.inference.distribution.UniformDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Reads a normal distribution model from a DOM Document element.
 */
public class UniformDistributionModelParser extends AbstractXMLObjectParser {

    public static final String UNIFORM_DISTRIBUTION_MODEL = "uniformDistributionModel";
    public static final String LOWER = "lower";
    public static final String UPPER = "upper";

    public String getParserName() {
        return UNIFORM_DISTRIBUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter upperParam;
        Parameter lowerParam;

        XMLObject cxo = xo.getChild(LOWER);
        if (cxo.getChild(0) instanceof Parameter) {
            lowerParam = (Parameter) cxo.getChild(Parameter.class);
        } else {
            lowerParam = new Parameter.Default(cxo.getDoubleChild(0));
        }

        cxo = xo.getChild(UPPER);
        if (cxo.getChild(0) instanceof Parameter) {
            upperParam = (Parameter) cxo.getChild(Parameter.class);
        } else {
            upperParam = new Parameter.Default(cxo.getDoubleChild(0));
        }

        return new UniformDistributionModel(lowerParam, upperParam);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(LOWER,
                    new XMLSyntaxRule[]{
                            new XORRule(
                                    new ElementRule(Parameter.class),
                                    new ElementRule(Double.class)
                            )}
            ),
            new ElementRule(UPPER,
                    new XMLSyntaxRule[]{
                            new XORRule(
                                    new ElementRule(Parameter.class),
                                    new ElementRule(Double.class)
                            )}
            )
    };

    public String getParserDescription() {
        return "Describes a uniform distribution with a given lower and upper bounds ";
    }

    public Class getReturnType() {
        return UniformDistributionModel.class;
    }
}
