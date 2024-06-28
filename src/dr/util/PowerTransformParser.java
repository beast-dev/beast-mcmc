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
 * @author Gabriel Hassler
 * @author Marc Suchard
 */

public class PowerTransformParser extends AbstractXMLObjectParser {
    public static final String POWER_TRANSFORM = "powerTransform";
    private static final String POWER = "power";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Transform.ParsedTransform parsedTransform = (Transform.ParsedTransform) xo.getChild(Transform.ParsedTransform.class);
        if (!(parsedTransform.transform instanceof Transform.PowerTransform)) {
            throw new XMLParseException("The '" + TransformParsers.TYPE + "' attribute of the " +
                    TransformParsers.TRANSFORM + " xml element must be '" + Transform.Type.POWER.getName() + "'.");
        }

        double power = xo.getDoubleAttribute(POWER);

        parsedTransform.transform = new Transform.PowerTransform(power);
        return parsedTransform;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(POWER),
                new ElementRule(Transform.ParsedTransform.class)
        };
    }

    @Override
    public String getParserDescription() {
        return "Raises a parameter to a given power.";
    }

    @Override
    public Class getReturnType() {
        return Transform.ParsedTransform.class;
    }

    @Override
    public String getParserName() {
        return POWER_TRANSFORM;
    }
}
