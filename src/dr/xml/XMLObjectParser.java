/*
 * XMLObjectParser.java
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

import java.util.Map;

public interface XMLObjectParser {

    /**
     *
     * @return  (Java) class of parsed element (i.e. class of object retuned in  parseXMLObject)
     */
    Class getReturnType();

    /**
     * @param store contains all named objects that have already been parsed.
     */
    Object parseXMLObject(XMLObject xo, String id, Map<String, XMLObject> store, boolean strictXML) throws XMLParseException;

    /**
     *
     * @return Parser name, which is identical to name of xml element parsed by it.
     */
    String getParserName();

    /**
     *
     * @return A list of parser name synonyms (including name returned by getParserName)
     */
    String[] getParserNames();

    /**
     *
     * @return  Human readable description of xml element parsed by parser.
     */
    String getParserDescription();

    /**
     *
     * @return true if an example is available.
     */
    boolean hasExample();

    /**
     *
     * @return element example
     */
    String getExample();


    /**
     * @return a description of this parser as HTML.
     */
    String toHTML(XMLDocumentationHandler handler);

    String toWiki(XMLDocumentationHandler handler);

    String toMarkdown(XMLDocumentationHandler handler);

    /**
     * @return an array of syntax rules required by this element.
     * Order is not important.
     */
    XMLSyntaxRule[] getSyntaxRules();

    /**
     * True if a top level element with that name is allowed.
     * @param elementName
     * @return
     */
    boolean isAllowed(String elementName);

    boolean hasSyntaxRules();
}
