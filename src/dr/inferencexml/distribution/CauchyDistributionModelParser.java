/*
 * NormalDistributionModelParser.java
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

import dr.inference.distribution.CauchyDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

public class CauchyDistributionModelParser extends AbstractXMLObjectParser {

    public static final String CAUCHY_DISTRIBUTION_MODEL = "cauchyDistributionModel";
    public static final String MEDIAN = "median";
    public static final String SCALE = "scale";

    public String getParserName() {
        return CAUCHY_DISTRIBUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final Parameter medianParam;
        final Parameter scaleParam;

        XMLObject cxo = xo.getChild(MEDIAN);
        if (cxo.getChild(0) instanceof Parameter) {
            medianParam = (Parameter) cxo.getChild(Parameter.class);
        } else {
            medianParam = new Parameter.Default(cxo.getDoubleChild(0));
        }

        cxo = xo.getChild(SCALE);
        if (cxo.getChild(0) instanceof Parameter) {
            scaleParam = (Parameter) cxo.getChild(Parameter.class);
        } else {
            scaleParam = new Parameter.Default(cxo.getDoubleChild(0));
        }

        return new CauchyDistributionModel(medianParam, scaleParam);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MEDIAN,
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
    };

    public String getParserDescription() {
        return "Describes a Cauchy distribution with a given median and scale " +
                "that can be used in a distributionLikelihood element";
    }

    public Class getReturnType() {
        return CauchyDistributionModel.class;
    }

}
