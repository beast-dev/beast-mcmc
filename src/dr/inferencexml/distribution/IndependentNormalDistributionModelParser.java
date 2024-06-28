/*
 * IndependentNormalDistributionModelParser.java
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

import dr.inference.distribution.IndependentNormalDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import static dr.inference.distribution.IndependentNormalDistributionModel.INDEPENDENT_NORMAL_DISTRIBUTION_MODEL;

public class IndependentNormalDistributionModelParser extends AbstractXMLObjectParser {
    public static String MEAN = "mean";
    public static String VARIANCE = "variance";
    public static String PRECISION = "precision";
    public static String DATA = "data";
    public static String ID = "id";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String id = xo.getStringAttribute(ID);
        Parameter mean = (Parameter) xo.getChild(MEAN).getChild(Parameter.class);
        Parameter precision = null;
        if (xo.getChild(PRECISION) != null) {
            precision = (Parameter) xo.getChild(PRECISION).getChild(Parameter.class);
        }
        Parameter variance = null;
        if (xo.getChild(VARIANCE) != null) {
            variance = (Parameter) xo.getChild(VARIANCE).getChild(Parameter.class);
        }
        Parameter data = (Parameter) xo.getChild(DATA).getChild(Parameter.class);

        return new IndependentNormalDistributionModel(id, mean, variance, precision, data);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(MEAN, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new XORRule(
                    new ElementRule(VARIANCE, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }),
                    new ElementRule(PRECISION, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    })
            ),
            new ElementRule(DATA, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            })
    };

    @Override
    public String getParserDescription() {
        return "A series of independent normal distribution models";
    }

    @Override
    public Class getReturnType() {
        return IndependentNormalDistributionModel.class;
    }

    @Override
    public String getParserName() {
        return INDEPENDENT_NORMAL_DISTRIBUTION_MODEL;
    }
}
