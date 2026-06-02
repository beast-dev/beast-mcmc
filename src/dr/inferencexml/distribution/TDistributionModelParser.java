/*
 * TDistributionModelParser.java
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

import dr.inference.distribution.TDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Reads a normal distribution model from a DOM Document element.
 */
public class TDistributionModelParser extends AbstractXMLObjectParser {

    public static final String T_DISTRIBUTION_MODEL = "tDistributionModel";
    public static final String LOCATION = "location";
    public static final String SCALE = "scale";
    public static final String DF = "df";

    public String getParserName() {
        return T_DISTRIBUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter locationParam;
        Parameter scaleParam;
        Parameter dfParam;

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


        cxo = xo.getChild(DF);
        if (cxo.getChild(0) instanceof Parameter) {
            dfParam = (Parameter) cxo.getChild(Parameter.class);
        } else {
            dfParam = new Parameter.Default(cxo.getDoubleChild(0));
        }
        return new TDistributionModel(locationParam, scaleParam, dfParam);
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
            ),
            new ElementRule(DF,
                    new XMLSyntaxRule[]{
                            new XORRule(
                                    new ElementRule(Parameter.class),
                                    new ElementRule(Double.class)
                            )}
            )
    };

    public String getParserDescription() {
        return "Describes a normal distribution with a given mean and standard deviation " +
                "that can be used in a distributionLikelihood element";
    }

    public Class getReturnType() {
        return TDistributionModel.class;
    }

}
