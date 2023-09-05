/*
 * SitePatternsParser.java
 *
 * Copyright (c) 2002-2021 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.alignment.TaxaFilteredSitePatterns;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Marc A. Suchard
 */
public class TaxaFilteredSitePatternsParser extends AbstractXMLObjectParser {

    private static final String PATTERNS = "taxaFilteredPatterns";
    private static final String OPERATION = "operation";

    public String getParserName() {
        return PATTERNS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        SitePatterns original = (SitePatterns) xo.getChild(SitePatterns.class);
        Taxon taxon = (Taxon) xo.getChild(Taxon.class);
        String operation = xo.getStringAttribute(OPERATION);

        List<Taxon> include = null;
        List<Taxon> exclude = null;

        if (operation.equalsIgnoreCase("include")) {

            include = new ArrayList<>();
            include.add(taxon);

        } else if (operation.equalsIgnoreCase("exclude")) {

            exclude = new ArrayList<>();
            exclude.add(taxon);

        } else {
            throw new XMLParseException("Unknown operation");
        }

        return new TaxaFilteredSitePatterns(original, include, exclude);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(OPERATION),
            new ElementRule(SitePatterns.class),
            new ElementRule(Taxon.class),
    };

    public String getParserDescription() {
        return "A filtered site pattern list";
    }

    public Class getReturnType() {
        return TaxaFilteredSitePatterns.class;
    }
}
