/*
 * PowerTransformParser.java
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

package dr.util;

import dr.xml.*;


/**
 * @author Filippo Monti
 */

public class AffineTransformParser extends AbstractXMLObjectParser {
    public static final String AFFINE_TRANSFORM = "affineTransform";
    private static final String DIMENSION = "dim";

    private static final String LOCATION = "location";
    private static final String SCALE = "scale";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double location = xo.hasAttribute(LOCATION) ? xo.getDoubleAttribute(LOCATION) : 0.0;
        double scale = xo.hasAttribute(SCALE) ? xo.getDoubleAttribute(SCALE) : 1.0;
        if (scale <= 0.0) {
            throw new XMLParseException("The '" + SCALE + "' attribute of the " +
                    TransformParsers.TRANSFORM + " xml element must be greater than 0.");
        }

        Transform.ParsedTransform parsedTransform = new Transform.ParsedTransform();
        if (xo.hasAttribute(DIMENSION)) {
            int dim = xo.getIntegerAttribute(DIMENSION);
            if (dim < 1) {
                throw new XMLParseException("The '" + DIMENSION + "' attribute of the " +
                        TransformParsers.TRANSFORM + " xml element must be greater than 0.");
            }
            parsedTransform.transform = new Transform.Array(new Transform.AffineTransform(location, scale), dim, null);
        } else {
            parsedTransform.transform = new Transform.AffineTransform(location, scale);
        }
        return parsedTransform;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(LOCATION, true),
                AttributeRule.newDoubleRule(SCALE, true),
                AttributeRule.newIntegerRule(DIMENSION, true),
        };
    }

    @Override
    public String getParserDescription() {
        return "Computes: (parameter - location) / scale.";
    }

    @Override
    public Class getReturnType() {
        return Transform.ParsedTransform.class;
    }

    @Override
    public String getParserName() {
        return AFFINE_TRANSFORM;
    }
}
