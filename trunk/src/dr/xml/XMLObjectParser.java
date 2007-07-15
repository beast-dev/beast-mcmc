/*
 * XMLObjectParser.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

public interface XMLObjectParser {

	/**
	 * @param store contains all named objects that have already been parsed.
	 */
	Object parseXMLObject(XMLObject xo, String id, ObjectStore store) throws XMLParseException;

	String getParserName();

    String[] getParserNames();

	String getParserDescription();

	boolean hasExample();

	String getExample();


	/**
	 * @return a description of this parser as HTML.
	 */
	String toHTML(XMLDocumentationHandler handler);

	String toWiki(XMLDocumentationHandler handler);
	
	Class getReturnType();

	XMLSyntaxRule[] getSyntaxRules();
}
