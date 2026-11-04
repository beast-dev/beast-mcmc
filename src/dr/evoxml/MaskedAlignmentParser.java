/*
 * MaskedPatternsParser.java
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

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.MaskedAlignment;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.SiteList;
import dr.xml.*;

/**
 * @author Andrew Rambaut
 *
 */
public class MaskedAlignmentParser extends AbstractXMLObjectParser {

    public static final String MASKED_ALIGNMENT = "maskedAlignment";
    public static final String MASK = "mask";
    public static final String NEGATIVE = "negative";
    public static final String INVERSE = "inverse";

    public String getParserName() { return MASKED_ALIGNMENT; }

    /**
     * Parses an alignment element and returns a masked alignment object.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Alignment alignment = (Alignment)xo.getChild(Alignment.class);

        boolean inverseMask = xo.getBooleanAttribute(INVERSE, false) || xo.getBooleanAttribute(NEGATIVE, false);
        String maskString = (String)xo.getElementFirstChild(MASK);

        boolean[] mask = new boolean[alignment.getSiteCount()];
        int k = 0;
        int onCount = 0;
        for (char c : maskString.toCharArray()) {
            if (Character.isDigit(c)) {
                if (k >= mask.length) {
                    break;
                }
                mask[k] = (c == '0' ? inverseMask : !inverseMask);
                if (mask[k]) {
                    onCount += 1;
                }
                k++;
            }
        }

        if (k != mask.length) {
            throw new XMLParseException("The mask needs to be the same length as the alignment (spaces are ignored)");
        }

        if (onCount == 0) {
            throw new XMLParseException("The mask needs to have at least one site unmasked");
        }

        return new MaskedAlignment(alignment, mask);
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(INVERSE, true),
            AttributeRule.newBooleanRule(NEGATIVE, true),
            new ElementRule(SiteList.class),
            new ElementRule(MASK, String.class, "A parameter of 1s and 0s that represent included and excluded sites")
    };

    public String getParserDescription() {
        return "Applies a mask to an alignment so that masked sites are ambiguous.";
    }

    public Class getReturnType() { return MaskedAlignment.class; }

}
