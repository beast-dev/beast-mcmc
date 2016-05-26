/*
 * SimpleNodeParser.java
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

package dr.evoxml;

import dr.evolution.tree.SimpleNode;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.util.Attributable;
import dr.util.Attribute;
import dr.xml.*;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 * @version $Id: SimpleNodeParser.java,v 1.2 2005/05/24 20:25:59 rambaut Exp $
 */
public class SimpleNodeParser extends AbstractXMLObjectParser {

    public final static String NODE = "node";
    public final static String HEIGHT = "height";
    public final static String RATE = "rate";

    public String getParserName() { return NODE; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        SimpleNode node = new SimpleNode();

        Taxon taxon = null;

        if (xo.hasAttribute(HEIGHT)) {
            node.setHeight(xo.getDoubleAttribute(HEIGHT));
        }

        if (xo.hasAttribute(RATE)) {
            node.setRate(xo.getDoubleAttribute(RATE));
        }

        for (int i = 0; i < xo.getChildCount(); i++) {
            Object child = xo.getChild(i);
            if (child instanceof dr.evolution.tree.SimpleNode) {
                node.addChild((dr.evolution.tree.SimpleNode)child);
            } else if (child instanceof Taxon) {
                taxon = (Taxon)child;
            } else if (child instanceof Date) {
                node.setAttribute("date", child);
            } else if (child instanceof Attribute) {
                Attribute attr = (Attribute)child;
                String name = attr.getAttributeName();
                Object value = attr.getAttributeValue();
                node.setAttribute(name, value);
            } else if (child instanceof Attribute[]) {
                Attribute[] attrs = (Attribute[])child;
                for (int j =0; j < attrs.length; j++) {
                    String name = attrs[j].getAttributeName();
                    Object value = attrs[j].getAttributeValue();
                    node.setAttribute(name, value);
                }
            } else if (child instanceof XMLObject) {
                XMLObject xoc = (XMLObject)child;
                if (xoc.getName().equals(Attributable.ATTRIBUTE)) {
                    node.setAttribute(xoc.getStringAttribute(Attributable.NAME), xoc.getAttribute(Attributable.VALUE));
                } else {
                    throw new XMLParseException("Unrecognized element" + xoc.getName() + " found in node element!");
                }
            } else {
                throw new XMLParseException("Unrecognized element found in node element!");
            }
        }

        if (taxon != null) {
            node.setTaxon(taxon);
        }

        return node;
    }

    public String getParserDescription() {
        return "This element represents a node in a tree.";
    }

    public Class getReturnType() { return SimpleNode.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
        AttributeRule.newDoubleRule(HEIGHT, true, "the age of the node"),
        AttributeRule.newDoubleRule(RATE, true, "the relative rate of evolution at this node - default is 1.0"),
        new XORRule(
            new ElementRule(Taxon.class, "The taxon of this leaf node"),
            new ElementRule(SimpleNode.class, "The children of this internal node", 2, 2))

    };
}
