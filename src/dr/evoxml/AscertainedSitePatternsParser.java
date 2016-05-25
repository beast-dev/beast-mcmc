/*
 * AscertainedSitePatternsParser.java
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

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.AscertainedSitePatterns;
import dr.evolution.alignment.PatternList;
import dr.evolution.util.TaxonList;
import dr.xml.*;
import dr.util.Citable;

import java.util.logging.Logger;

/**
 * Package: AscertainedSitePatternsParser
 * Description:
 * <p/>
 * <p/>
 * Created by
 * Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * Date: Mar 10, 2008
 * Time: 3:06:52 PM
 */
public class AscertainedSitePatternsParser extends AbstractXMLObjectParser {

    public static final String APATTERNS = "ascertainedPatterns";
    public static final String FROM = "from";
    public static final String TO = "to";
    public static final String EVERY = "every";
    public static final String TAXON_LIST = "taxonList";
    public static final String INCLUDE = "includePatterns";
    public static final String EXCLUDE = "excludePatterns";

    public String getParserName() {
        return APATTERNS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Alignment alignment = (Alignment) xo.getChild(Alignment.class);
        XMLObject xoc;
        TaxonList taxa = null;

        int from = -1;
        int to = -1;
        int every = xo.getAttribute(EVERY, 1);
        if (every <= 0) throw new XMLParseException("illegal 'every' attribute in patterns element");

        int startInclude = -1;
        int stopInclude = -1;
        int startExclude = -1;
        int stopExclude = -1;

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


        if (xo.hasChildNamed(TAXON_LIST)) {
            taxa = (TaxonList) xo.getElementFirstChild(TAXON_LIST);
        }

        if (from > alignment.getSiteCount())
            throw new XMLParseException("illegal 'from' attribute in patterns element");

        if (to > alignment.getSiteCount())
            throw new XMLParseException("illegal 'to' attribute in patterns element");

        if (from < 0) from = 0;
        if (to < 0) to = alignment.getSiteCount() - 1;

//        if (xo.hasAttribute(XMLParser.ID)) {
            Logger.getLogger("dr.evoxml").info("Creating ascertained site patterns '" + xo.getId() + "' from positions " +
                    Integer.toString(from + 1) + "-" + Integer.toString(to + 1) +
                    " of alignment '" + alignment.getId() + "'");
            if (every > 1) {
                Logger.getLogger("dr.evoxml").info("  only using every " + every + " site");
            }
//        }

        if (xo.hasChildNamed(INCLUDE)) {
            xoc = xo.getChild(INCLUDE);
            if (xoc.hasAttribute(FROM) && xoc.hasAttribute(TO)) {
                startInclude = xoc.getIntegerAttribute(FROM) - 1;
                stopInclude = xoc.getIntegerAttribute(TO);
            } else {
                throw new XMLParseException("both from and to attributes are required for includePatterns");
            }

            if (startInclude < 0 || stopInclude < startInclude) {
                throw new XMLParseException("invalid 'from' and 'to' attributes in includePatterns");
            }
            Logger.getLogger("dr.evoxml").info("\tAscertainment: Patterns in columns " + (startInclude + 1) + " to " + (stopInclude) + " are only possible. ");
        }

        if (xo.hasChildNamed(EXCLUDE)) {
            xoc = xo.getChild(EXCLUDE);
            if (xoc.hasAttribute(FROM) && xoc.hasAttribute(TO)) {
                startExclude = xoc.getIntegerAttribute(FROM) - 1;
                stopExclude = xoc.getIntegerAttribute(TO);
            } else {
                throw new XMLParseException("both from and to attributes are required for excludePatterns");
            }

            if (startExclude < 0 || stopExclude < startExclude) {
                throw new XMLParseException("invalid 'from' and 'to' attributes in includePatterns");
            }
            Logger.getLogger("dr.evoxml").info("\tAscertainment: Patterns in columns " + (startExclude + 1) + " to " + (stopExclude) + " are not possible. ");
        }

        AscertainedSitePatterns patterns = new AscertainedSitePatterns(alignment, taxa,
                from, to, every,
                startInclude, stopInclude,
                startExclude, stopExclude);

        Logger.getLogger("dr.evoxml").info("\tThere are " + patterns.getPatternCount() + " patterns in total.");
        
        Logger.getLogger("dr.evoxml").info("\tPlease cite:\n" + Citable.Utils.getCitationString(patterns));

        return patterns;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newIntegerRule(FROM, true, "The site position to start at, default is 1 (the first position)"),
            AttributeRule.newIntegerRule(TO, true, "The site position to finish at, must be greater than <b>" + FROM + "</b>, default is length of given alignment"),
            AttributeRule.newIntegerRule(EVERY, true, "Determines how many sites are selected. A value of 3 will select every third site starting from <b>" + FROM + "</b>, default is 1 (every site)"),
            new ElementRule(TAXON_LIST,
                    new XMLSyntaxRule[]{new ElementRule(TaxonList.class)}, true),
            new ElementRule(Alignment.class),
            new ContentRule("<includePatterns from=\"Z\" to=\"X\"/>"),
            new ContentRule("<excludePatterns from=\"Z\" to=\"X\"/>")
    };

    public String getParserDescription() {
        return "A weighted list of the unique site patterns (unique columns) in an alignment.";
    }

    public Class getReturnType() {
        return PatternList.class;
    }

}
