/*
 * ReportSitePatternMapParser.java
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

package dr.evoxml;

import dr.evolution.alignment.SitePatterns;
import dr.util.FileHelpers;
import dr.xml.*;

import java.io.PrintWriter;

import static dr.inferencexml.loggers.LoggerParser.getLogFile;

/**
 * @author Marc A. Suchard
 */
public class ReportSitePatternMapParser extends AbstractXMLObjectParser {

    private static final String SITE_PATTERN_MAP = "reportSitePatternMap";
    private static final String FILE_NAME = FileHelpers.FILE_NAME;

    public String getParserName() {
        return SITE_PATTERN_MAP;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final PrintWriter pw = getLogFile(xo, getParserName());
        SitePatterns sitePatterns = (SitePatterns) xo.getChild(SitePatterns.class);

        int totalSiteCount = sitePatterns.getSiteCount();
        int totalPatternCount = sitePatterns.getPatternCount();

        pw.println("Site count    = " + totalSiteCount);
        pw.println("Pattern count = " + totalPatternCount);
        pw.println("Map (site -> pattern):");

        for (int s = 0; s < totalSiteCount; ++s) {
            int p = sitePatterns.getPatternIndex(s);
            pw.println((s + 1) + " -> " + (p + 1));
        }

        pw.close();
        return null;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(FILE_NAME, true),
            new ElementRule(SitePatterns.class),
    };

    public String getParserDescription() {
        return "A weighted list of the unique site patterns (unique columns) in an alignment.";
    }

    public Class getReturnType() {
        return Object.class;
    }

}
