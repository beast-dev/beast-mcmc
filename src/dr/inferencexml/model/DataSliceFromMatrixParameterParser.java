/*
 * SliceDataFromMatrixParameterParser.java
 *
 * Copyright © 2002-2025 the BEAST Development Team
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

import dr.inference.model.DataSliceFromMatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.util.Attribute;
import dr.xml.*;

import static dr.inference.model.DataSliceFromMatrixParameter.SliceDirection;

public class DataSliceFromMatrixParameterParser extends AbstractXMLObjectParser {

    private static final String SLICE_DATA = "sliceData";
    private static final String DIRECTION = "by";

    @Override
    public String getParserName() {
        return SLICE_DATA;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        MatrixParameterInterface parameter = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
        String directionString = xo.getStringAttribute(DIRECTION);
        SliceDirection direction = parseDirection(directionString);

        int dim = direction.getCount(parameter);
        Parameter[] attributes = new Parameter[dim];
        for (int i = 0; i < dim; ++i) {
            attributes[i] = new DataSliceFromMatrixParameter(parameter, i, direction);
        }

        return attributes;
    }

    private SliceDirection parseDirection(String name) throws XMLParseException {
        for (SliceDirection d : SliceDirection.values()) {
            if (name.equalsIgnoreCase(d.getName())) {
                return d;
            }
        }
        throw new XMLParseException("Unknown slice-direction '" + name + "'");
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                AttributeRule.newStringRule(DIRECTION),
                new ElementRule(MatrixParameterInterface.class),
        };
    }

    @Override
    public String getParserDescription() {
        return SLICE_DATA;
    }

    @Override
    public Class getReturnType() { return Parameter[].class; }
}
