/*
 * GumbelDistributionModelParser.java
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

import dr.inference.distribution.GumbelDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Reads a (minimum-oriented) Gumbel type I distribution model from a DOM Document element.
 */
public class GumbelDistributionModelParser extends AbstractXMLObjectParser {

    public static final String GUMBEL_DISTRIBUTION_MODEL = "gumbelDistributionModel";
    public static final String LOCATION = "location";
    public static final String SCALE = "scale";

    public String getParserName() {
        return GUMBEL_DISTRIBUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter locationParam;
        Parameter scaleParam;

        XMLObject cxo = xo.getChild(LOCATION);
        if (cxo.getChild(0) instanceof Parameter) {
            locationParam = (Parameter) cxo.getChild(Parameter.class);
        } else {
            locationParam = new Parameter.Default(cxo.getDoubleChild(0));
        }

        cxo = xo.getChild(SCALE);
        if (cxo.getChild(0) instanceof Parameter) {
            scaleParam = (Parameter) cxo.getChild(Parameter.class);
        } else {
            scaleParam = new Parameter.Default(cxo.getDoubleChild(0));
        }

        return new GumbelDistributionModel(locationParam, scaleParam);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(LOCATION,
                    new XMLSyntaxRule[]{
                            new XORRule(
                                    new ElementRule(Parameter.class),
                                    new ElementRule(Double.class)
                            )}
            ),
            new ElementRule(SCALE,
                    new XMLSyntaxRule[]{
                            new XORRule(
                                    new ElementRule(Parameter.class),
                                    new ElementRule(Double.class)
                            )}
            )
    };

    public String getParserDescription() {
        return "Describes a (minimum-oriented) Gumbel type I distribution with a given location and scale " +
                "that can be used in a distributionLikelihood element. If R ~ Exponential(mean = theta), " +
                "then log(R) ~ this distribution with location = log(theta), scale = 1.";
    }

    public Class getReturnType() {
        return GumbelDistributionModel.class;
    }

}
