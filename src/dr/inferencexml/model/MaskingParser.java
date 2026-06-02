/*
 * MaskingParser.java
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

package dr.inferencexml.model;

import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Marc A. Suchard
 */
public class MaskingParser extends AbstractXMLObjectParser {

    private static final String MASKING = "masking";
    private static final String DIMENSION = "dimension";
    private static final String NAME = "name";
    private static final String COMPLEMENT = "complement";

    public class MaskingParameter extends Parameter.Default {

        private final Parameter parameter;

        public MaskingParameter(int dimension, Parameter parameter) {
            super(dimension, 0.0);
            this.parameter = parameter;
        }
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        List<String> names = new ArrayList<>();
        for (XMLObject cxo : xo.getAllChildren(DIMENSION)) {
            names.add(cxo.getStringAttribute(NAME));
        }

        Map<String,Integer> dimMap = new HashMap<>();
        for (int i = 0; i < parameter.getDimension(); ++i) {
            String n = parameter.getDimensionName(i);
            dimMap.put(n, i);
        }

        Parameter result = new MaskingParameter(parameter.getDimension(), parameter);

        for (String name : names) {
            if (dimMap.containsKey(name)) {
                int index = dimMap.get(name);
                result.setParameterValue(index, 1);
            } else {
                throw new XMLParseException("Dimension '" + name + "' not found");

            }
        }

        if (xo.getAttribute(COMPLEMENT, false)) {
            for (int i = 0; i < parameter.getDimension(); ++i) {
                result.setParameterValue(i,1 - result.getParameterValue(i));
            }
        }

        return result;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Parameter.class),
            AttributeRule.newBooleanRule(COMPLEMENT, true),
            new ElementRule(DIMENSION, new XMLSyntaxRule[]{
                    AttributeRule.newStringRule(NAME),
            }, 1, Integer.MAX_VALUE),
    };

    public String getParserDescription() {
        return "A binary masking as a parameter.";
    }

    public Class getReturnType() {
        return MaskingParameter.class;
    }

    public String getParserName() {
        return MASKING;
    }
}
