/*
 * WanderingTaxonLoggerParser.java
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

package dr.evomodelxml;

import dr.xml.*;
import dr.evomodel.tree.WanderingTaxonLogger;
import dr.evolution.util.Taxon;

/**
 * @author Marc A. Suchard
 */
public class WanderingTaxonLoggerParser extends AbstractXMLObjectParser {

    public static final String NAME = "name";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(NAME, xo.getId());

        WanderingTaxonLogger.Relative relative = WanderingTaxonLogger.Relative.SISTER;
        if (xo.getAttribute(WanderingTaxonLogger.RELATIVE,"sister").equalsIgnoreCase("parent")) {
            relative = WanderingTaxonLogger.Relative.PARENT;
        }

        Taxon taxon = (Taxon) xo.getChild(Taxon.class);

        return new WanderingTaxonLogger(name,taxon,relative);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                AttributeRule.newStringRule(NAME,true),
                AttributeRule.newStringRule(WanderingTaxonLogger.RELATIVE,true),
                new ElementRule(Taxon.class),
        };
    }

    public String getParserDescription() {
        return null;
    }

    public Class getReturnType() {
        return WanderingTaxonLogger.class;
    }

    public String getParserName() {
        return WanderingTaxonLogger.WANDERER;
    }
}
