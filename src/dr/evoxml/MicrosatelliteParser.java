/*
 * MicrosatelliteParser.java
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

package dr.evoxml;

import dr.evolution.datatype.Microsatellite;
import dr.xml.*;


/**
 * @author Chieh-Hsi Wu
 *
 * Microsatellite data type parser
 *
 */
public class MicrosatelliteParser extends AbstractXMLObjectParser {

    public static final String MICROSAT = "microsatellite";
    public static final String MIN = "min";
    public static final String MAX = "max";
    public static final String UNIT_LENGTH = "unitLength";

    public String getParserName() {
        return MICROSAT;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        int min = xo.getIntegerAttribute(MIN);
        int max = xo.getIntegerAttribute(MAX);
        int unitLength = xo.hasAttribute(UNIT_LENGTH) ? xo.getIntegerAttribute(UNIT_LENGTH) : 1;
        String name = xo.getId();

        return new Microsatellite(name, min, max, unitLength);
    }

    public String getParserDescription() {
        return "This element represents a microsatellite data type.";
    }

    public String getExample() {
        return "<microsatellite min=\"0\" max=\"20\" unitLength=\"2\"/>";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new AndRule(new XMLSyntaxRule[]{
                    AttributeRule.newStringRule(XMLObject.ID),
                    AttributeRule.newIntegerRule(MIN),
                    AttributeRule.newIntegerRule(MAX),
                    AttributeRule.newIntegerRule(UNIT_LENGTH, true)})
        };
    }


    public Class getReturnType() {
        return Microsatellite.class;
    }
}