/*
 * PatternSubSetParser.java
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
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: PatternsParser.java,v 1.2 2005/05/24 20:25:59 rambaut Exp $
 */
public class PatternSubSetParser extends AbstractXMLObjectParser {

    public static final String PATTERNS_SUB_SET = "patternSubSet";
    public static final String SUB_SET = "subSet";
    public static final String SUB_SET_COUNT = "subSetCount";

    public String getParserName() {
        return PATTERNS_SUB_SET;
    }

    /**
     * Parses a patterns element and returns a patterns object.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        SiteList patterns = (SiteList) xo.getChild(SiteList.class);

        int subSet = 0;
        int subSetCount = 0;

        if (xo.hasAttribute(SUB_SET)) {
            subSet = xo.getIntegerAttribute(SUB_SET) - 1;
            if (subSet < 0)
                throw new XMLParseException("illegal 'subSet' attribute in patterns element");
        }

        if (xo.hasAttribute(SUB_SET_COUNT)) {
            subSetCount = xo.getIntegerAttribute(SUB_SET_COUNT);
            if (subSetCount < 0)
                throw new XMLParseException("illegal 'subSetCount' attribute in patterns element");
        }

        Patterns subPatterns = new Patterns(patterns, 0, 0, 1, subSet, subSetCount);
        
        if (xo.hasAttribute(XMLParser.ID)) {
            final Logger logger = Logger.getLogger("dr.evoxml");
            logger.info("Pattern subset '" + xo.getId() + "' created from '" + patterns.getId() +"' ("+(subSet+1)+"/"+subSetCount+")");
            logger.info("  pattern count = " + subPatterns.getPatternCount());
        }
        
        return subPatterns;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(SUB_SET, true, "Which subset of patterns to use (out of subSetCount)"),
            AttributeRule.newIntegerRule(SUB_SET_COUNT, true, "The number of subsets"),

            new ElementRule(SiteList.class)
    };

    public String getParserDescription() {
        return "A weighted list of the unique site patterns (unique columns) in an alignment.";
    }

    public Class getReturnType() {
        return PatternList.class;
    }

}
