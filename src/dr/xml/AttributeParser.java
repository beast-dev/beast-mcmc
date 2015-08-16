/*
 * AttributeParser.java
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

import dr.util.Attribute;

public class AttributeParser extends AbstractXMLObjectParser {

	public final static String ATTRIBUTE = "attr";
	public final static String NAME = "name";
	public final static String VALUE = "value";

	public String getParserName() { return ATTRIBUTE; }
		
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final String name = xo.getStringAttribute(NAME);
        if( xo.hasAttribute(VALUE) ) {
            return new Attribute.Default<Object>(name, xo.getAttribute(VALUE));
        }
        final Object value = xo.getChild(0);

        return new Attribute.Default<Object>(name, value);
	}
	
	//************************************************************************
	// AbstractXMLObjectParser implementation
	//************************************************************************

	public String getParserDescription() {
		return "This element represents a name/value pair.";
	}
	
	public Class getReturnType() { return Attribute.class; }

	public XMLSyntaxRule[] getSyntaxRules() { return rules; }

	private final XMLSyntaxRule[] rules = {
		new StringAttributeRule("name", "The name to give to this attribute"),
		new ElementRule(Object.class )
	};
}
