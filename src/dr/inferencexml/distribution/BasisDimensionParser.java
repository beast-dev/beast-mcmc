/*
 * GaussianProcessFieldParser.java
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

import dr.inference.distribution.RandomField;
import dr.inference.model.DesignMatrix;
import dr.math.distributions.gp.*;
import dr.xml.*;
import org.w3c.dom.NamedNodeMap;

import java.util.HashMap;
import java.util.Map;

public class BasisDimensionParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "basis";
    public static final String WEIGHTFUNCTION = "weightFunction";
    public static final String TYPE = "type";
    public static final String NORMALIZED = "normalized";
    public static final String UNITARY_VARIANCE = "unitaryVariance";

    public String getParserName() { return PARSER_NAME; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        GaussianProcessKernel kernel = (GaussianProcessKernel) xo.getChild(GaussianProcessKernel.class);

        DesignMatrix design = (DesignMatrix) xo.getChild(DesignMatrix.class);
        if (design != null) {
            WeightFunction weightFunction = null;
            if (xo.hasChildNamed(WEIGHTFUNCTION)) {
                XMLObject wfElement = xo.getChild(WEIGHTFUNCTION);

                String type;
                try {
                    type = (String) wfElement.getAttribute(TYPE);
                } catch (XMLParseException e) {
                    throw new RuntimeException(e);
                }

                Map<String, Double> params = new HashMap<>();
                NamedNodeMap attributeList = wfElement.getAttributes();
                for (int i = 0; i < attributeList.getLength(); i++) {
                    String attrName = attributeList.item(i).getNodeName();
                    if (!attrName.equalsIgnoreCase(TYPE)) {
                        params.put(attrName, Double.parseDouble(attributeList.item(i).getNodeValue()));
                    }
                }
                weightFunction = WeightFunction.WeightFunctionFactory.create(type, params);
            }
            if (xo.getAttribute(NORMALIZED, false)) {
                return new NormalizedBasisDimension(kernel, design, design, weightFunction);
            } else if (xo.getAttribute(UNITARY_VARIANCE, false)) {
                return new UnitaryVarianceBasisDimension(kernel, design, design, weightFunction);
            } else {
                return new BasisDimension(kernel, design, design, weightFunction);
            }

        } else {
            RandomField.WeightProvider weights = (RandomField.WeightProvider)
                    xo.getChild(RandomField.WeightProvider.class);
            if (weights == null) {
                throw new XMLParseException("Either a design matrix or weights must be specified for a basis dimension");
            }

            return new BasisDimension(kernel, weights);
        }
    }


    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(NORMALIZED, true),
            AttributeRule.newBooleanRule(UNITARY_VARIANCE, true),
            new XORRule(
                new ElementRule(DesignMatrix.class),
                new ElementRule(RandomField.WeightProvider.class)),
            new ElementRule(GaussianProcessKernel.class),
            new ElementRule(WeightFunction.class, true),
    };

    public String getParserDescription() { // TODO update
        return "Creates a basis for a Guassian process " +
                "using a specified kernel and design matrix.";
    }

    public Class getReturnType() { return RandomField.class; }
}
