/*
 * AttributesParser.java
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

public class AttributesParser extends AbstractXMLObjectParser {

	public final static String ATTRIBUTES = "attributes";
	public final static String NAMES = "names";
	public final static String VALUES = "values";

	public String getParserName() { return ATTRIBUTES; }
		
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		
		String[] names = ((XMLObject)xo.getChild(NAMES)).getStringArrayChild(0);
		String[] values =((XMLObject)xo.getChild(VALUES)).getStringArrayChild(0);
		
		if (names.length != values.length) {
			throw new XMLParseException("The number of names and values must match.");
		}
		
		Attribute[] attributes = new Attribute[names.length];
		for (int i =0; i < attributes.length; i++) {
			attributes[i] = new Attribute.Default(names[i], values[i]);
		}
		
		return attributes;
	}

	//************************************************************************
	// AbstractXMLObjectParser implementation
	//************************************************************************

	public String getParserDescription() {
		return "This element represents an array of name/value pairs.";
	}
	
	public Class getReturnType() { return Attribute[].class; }

	public XMLSyntaxRule[] getSyntaxRules() { return rules; }

	private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
		AttributeRule.newStringArrayRule("names"),
		AttributeRule.newStringArrayRule("values" )
	};
	
}
