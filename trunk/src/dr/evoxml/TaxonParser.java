/*
 * TaxonParser.java
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

package dr.evoxml;

import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.util.Attribute;
import dr.util.Identifiable;
import dr.xml.*;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 * @version $Id: TaxonParser.java,v 1.2 2005/05/24 20:25:59 rambaut Exp $
 */
public class TaxonParser extends AbstractXMLObjectParser {

    public final static String TAXON = "taxon";

    public String getParserName() { return TAXON; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Taxon taxon = new Taxon(xo.getStringAttribute(Identifiable.ID));

        for (int i = 0; i < xo.getChildCount(); i++) {
            Object child = xo.getChild(i);

            if (child instanceof Date) {
                taxon.setDate((Date)child);
            } else if (child instanceof Attribute) {
                Attribute attr = (Attribute)child;
                String name = attr.getAttributeName();
                Object value = attr.getAttributeValue();
                taxon.setAttribute(name, value);
            } else if (child instanceof Attribute[]) {
                Attribute[] attrs = (Attribute[])child;
                for (int j =0; j < attrs.length; j++) {
                    String name = attrs[j].getAttributeName();
                    Object value = attrs[j].getAttributeValue();
                    taxon.setAttribute(name, value);
                }
            } else {
                throw new XMLParseException("Unrecognized element found in taxon element!");
            }
        }

        return taxon;
    }

    public String getParserDescription() { return ""; }
    public Class getReturnType() { return Taxon.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
        new StringAttributeRule(Identifiable.ID, "A unique identifier for this taxon"),
        new ElementRule(dr.evolution.util.Date.class, true)
    };

}
