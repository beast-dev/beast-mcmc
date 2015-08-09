/*
 * PropertyParser.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.xml;

import dr.util.Property;

public class PropertyParser extends AbstractXMLObjectParser {

    public String getParserName() {
        return "property";
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Object object = xo.getChild(0);
        String name = xo.getStringAttribute("name");

        Property property;

        if (xo.hasAttribute("index")) {
            int index = xo.getIntegerAttribute("index");
            property = new Property(object, name, index);
        } else if (xo.hasAttribute("label")) {
            String label = xo.getStringAttribute("label");
            property = new Property(object, name, label);
        } else {
            property = new Property(object, name);
        }

        if (property.getGetter() == null)
            throw new XMLParseException("unknown property, " + name + ", for object, " + object + ", in property element");

        return property;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns an object representing the named property of the given child object.";
    }

    public Class getReturnType() {
        return Object.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule("name", "name of the property", "length"),
            new ElementRule(Object.class)
    };

}
