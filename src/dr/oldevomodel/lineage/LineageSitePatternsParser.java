/*
 * LineageSitePatternsParser.java
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

package dr.oldevomodel.lineage;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.PatternList;
import dr.evolution.util.TaxonList;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexander Alekseyenko
 * @author Jack O'Brien
 */
public class LineageSitePatternsParser extends AbstractXMLObjectParser {

    public static final String PATTERNS = "lineageSNPS";
    public static final String FROM = "from";
    public static final String TO = "to";
    public static final String EVERY = "every";
    public static final String TAXON_LIST = "taxonList";
    public static final String STRIP = "strip";
    public static final String UNIQUE = "unique";


    public String getParserName() {
        return PATTERNS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Alignment alignment = (Alignment) xo.getChild(Alignment.class);
        TaxonList taxa = null;

        int from = 0;
        int to = -1;
        int every = xo.getAttribute(EVERY, 1);

        boolean strip = xo.getAttribute(STRIP, true);

        boolean unique = xo.getAttribute(UNIQUE, true);

        if (xo.hasAttribute(FROM)) {
            from = xo.getIntegerAttribute(FROM) - 1;
            if (from < 0)
                throw new XMLParseException("illegal 'from' attribute in patterns element");
        }

        if (xo.hasAttribute(TO)) {
            to = xo.getIntegerAttribute(TO) - 1;
            if (to < 0 || to < from)
                throw new XMLParseException("illegal 'to' attribute in patterns element");
        }

        if (every <= 0) throw new XMLParseException("illegal 'every' attribute in patterns element");

        if (xo.hasChildNamed(TAXON_LIST)) {
            taxa = (TaxonList) xo.getElementFirstChild(TAXON_LIST);
        }

        if (from > alignment.getSiteCount())
            throw new XMLParseException("illegal 'from' attribute in patterns element");

        if (to > alignment.getSiteCount())
            throw new XMLParseException("illegal 'to' attribute in patterns element");

        LineageSitePatterns patterns = new LineageSitePatterns(alignment, taxa, from, to, every, strip, unique);

        int f = from + 1;
        int t = to + 1; // fixed a *display* error by adding + 1 for consistency with f = from + 1
        if (to == -1) t = alignment.getSiteCount();

        if (xo.hasAttribute(XMLParser.ID)) {
            final Logger logger = Logger.getLogger("dr.evoxml");
            logger.info("Site patterns '" + xo.getId() + "' created from positions " +
                    Integer.toString(f) + "-" + Integer.toString(t) +
                    " of alignment '" + alignment.getId() + "'");

            if (every > 1) {
                logger.info("  only using every " + every + " site");
            }
            logger.info("  pattern count = " + patterns.getPatternCount());
        }

        return patterns;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(FROM, true, "The site position to start at, default is 1 (the first position)"),
            AttributeRule.newIntegerRule(TO, true, "The site position to finish at, must be greater than <b>" + FROM + "</b>, default is length of given alignment"),
            AttributeRule.newIntegerRule(EVERY, true, "Determines how many sites are selected. A value of 3 will select every third site starting from <b>" + FROM + "</b>, default is 1 (every site)"),
            new ElementRule(TAXON_LIST,
                    new XMLSyntaxRule[]{new ElementRule(TaxonList.class)}, true),
            new ElementRule(Alignment.class),
            AttributeRule.newBooleanRule(STRIP, true, "Strip out completely ambiguous sites"),
            AttributeRule.newBooleanRule(UNIQUE, true, "Return a weight list of unique patterns"),
    };

    public String getParserDescription() {
        return "A weighted list of the unique site patterns (unique columns) in an alignment.";
    }

    public Class getReturnType() {
        return PatternList.class;
    }

}
