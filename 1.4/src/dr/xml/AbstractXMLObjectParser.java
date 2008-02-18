/*
 * AbstractXMLObjectParser.java
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

import org.w3c.dom.NamedNodeMap;

import java.io.PrintWriter;
import java.io.StringWriter;

public abstract class AbstractXMLObjectParser implements XMLObjectParser {

    public final Object parseXMLObject(XMLObject xo, String id, ObjectStore store) throws XMLParseException {

        this.store = store;

        if (hasSyntaxRules()) {
            XMLSyntaxRule[] rules = getSyntaxRules();
            for (int i =0; i < rules.length; i++) {
                if (!rules[i].isSatisfied(xo)) {
                    // System.err.println(this);
                    throw new XMLParseException("The '<" + getParserName() +
                            ">' element with id, '" + id +
                            "', is incorrectly constructed.\nThe following was expected:\n" +
                            rules[i].ruleString(xo));
                }
            }

            // Look for undeclared attributes and issue a warning
            final NamedNodeMap attributes = xo.getAttributes();
            for(int k = 0; k < attributes.getLength(); ++k) {
                String name = attributes.item(k).getNodeName();
                if( name.equals(XMLObject.ID) ) continue;

                for (int i = 0; i < rules.length; i++) {
                    if( rules[i].containsAttribute(name) ) {
                        name = null;
                        break;
                    }
                }
                if( name != null ) {
                    System.err.println("WARNING: unhandled attribute (typo?) " + name + " in " + xo);
                }
            }
        }

        try {
            return parseXMLObject(xo);
        } catch(XMLParseException xpe) {
            throw new XMLParseException("Error parsing '<" + getParserName() +
                    ">' element with id, '" + id + "':\n" +
                    xpe.getMessage());
        }
    }

    public String[] getParserNames() { return new String[] { getParserName() }; }

    public final void throwUnrecognizedElement(XMLObject xo) throws XMLParseException {
        throw new XMLParseException("Unrecognized element '<" + xo.getName() + ">' in element '<" + getParserName()+">'");
    }

    public abstract Object parseXMLObject(XMLObject xo) throws XMLParseException;
    /**
     * @return an array of syntax rules required by this element.
     * Order is not important.
     */
    public abstract XMLSyntaxRule[] getSyntaxRules();

    public abstract String getParserDescription();

    public abstract Class getReturnType();

    public final boolean hasExample() {
        return getExample() != null;
    }

    public String getExample() { return null; }

    public final ObjectStore getStore() {
        return store;
    }

    /**
     * @return a description of this parser as a string.
     */
    public final String toHTML(XMLDocumentationHandler handler) {
        StringBuffer buffer = new StringBuffer();
        buffer.append("<div id=\"").append(getParserName()).append("\" class=\"element\">\n");
        buffer.append("  <div class=\"elementheader\">\n");
        buffer.append("    <span class=\"elementname\">").append(getParserName()).append("</span> element\n");
        buffer.append("    <div class=\"description\">\n");
        buffer.append("      ").append(getParserDescription()).append("\n");
        buffer.append("    </div>\n");
        buffer.append("  </div>\n");
        if (hasSyntaxRules()) {
            XMLSyntaxRule[] rules = getSyntaxRules();
            buffer.append("  <div class=\"rules\">\n");
            for (int i = 0; i < rules.length; i++) {
                buffer.append(rules[i].htmlRuleString(handler));
            }
            buffer.append("  </div>\n");
        }
        buffer.append("<div class=\"example\">");
        StringWriter sw = new StringWriter();
        PrintWriter pw = new PrintWriter(sw);
        handler.outputExampleXML(pw, this);
        pw.flush(); pw.close();
        buffer.append(sw.toString());
        buffer.append("</div>\n");
        buffer.append("</div>\n");
        return buffer.toString();
    }

    /**
     * @return a description of this parser as a string.
     */
    public String toString() {
        StringBuffer buffer = new StringBuffer();
        buffer.append("\nELEMENT ").append(getParserName()).append("\n");

        if (hasSyntaxRules()) {
            XMLSyntaxRule[] rules = getSyntaxRules();
            for (int i = 0; i < rules.length; i++) {
                buffer.append("  ").append(rules[i].ruleString()).append("\n");
            }
        }
        return buffer.toString();
    }

    //************************************************************************
    // private methods
    //************************************************************************

    public final boolean hasSyntaxRules() {
        XMLSyntaxRule[] rules = getSyntaxRules();
        return (rules != null && rules.length > 0);
    }

    private ObjectStore store = null;
}
