/*
 * MaskedPatternsParser.java
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

import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.Patterns;
import dr.evolution.alignment.SiteList;
import dr.evolution.alignment.SitePatterns;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 *
 * @version $Id$
 */
public class MaskedPatternsParser extends AbstractXMLObjectParser {

    public static final String MASKED_PATTERNS = "maskedPatterns";
    public static final String MASK = "mask";
    public static final String NEGATIVE = "negative";
    public static final String INVERSE = "inverse";

    public String getParserName() { return MASKED_PATTERNS; }

    /**
     * Parses a patterns element and returns a patterns object.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        SiteList siteList = (SiteList)xo.getChild(SiteList.class);

        boolean inverseMask = xo.getBooleanAttribute(INVERSE, false) || xo.getBooleanAttribute(NEGATIVE, false);
        String maskString = (String)xo.getElementFirstChild(MASK);

        boolean[] mask = new boolean[siteList.getSiteCount()];
        int k = 0;
        for (char c : maskString.toCharArray()) {
            if (Character.isDigit(c)) {
                if (k >= mask.length) {
                    break;
                }
                mask[k] = (c == '0' ? inverseMask : !inverseMask);
                k++;
            }
        }

        if (k != mask.length) {
            throw new XMLParseException("The mask needs to be the same length as the alignment (spaces are ignored)");
        }

        throw new UnsupportedOperationException("This has not been implemented");

        // @todo - work out how to do this. Create a MaskedSitePatterns which can then be put into a SitePatterns
        // for compression etc. Or MaskedAlignment may make more sense here.
//        SitePatterns patterns = new MaskedSitePatterns(siteList, mask);
//
//        if (patterns == null) {
//            throw new XMLParseException("The mask needs include at least one pattern");
//        }
//
//        if (xo.hasAttribute(XMLParser.ID)) {
//            final Logger logger = Logger.getLogger("dr.evoxml");
//            logger.info("Site patterns '" + xo.getId() + "' created by masking alignment with id '" + siteList.getId() + "'");
//            logger.info("  pattern count = " + patterns.getPatternCount());
//        }
//
//        return patterns;
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(INVERSE, true),
            AttributeRule.newBooleanRule(NEGATIVE, true),
            new ElementRule(SiteList.class),
            new ElementRule(MASK, String.class, "A parameter of 1s and 0s that represent included and excluded sites")
    };

    public String getParserDescription() {
        return "A weighted list of the unique site patterns (unique columns) in an alignment.";
    }

    public Class getReturnType() { return PatternList.class; }

}
