/*
 * MergePatternsParser.java
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

import dr.evolution.alignment.*;
import dr.xml.*;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 * @version $Id: MergePatternsParser.java,v 1.1 2005/07/08 11:27:53 rambaut Exp $
 */
public class MergePatternsParser extends AbstractXMLObjectParser {

    public static final String MERGE_PATTERNS = "mergePatterns";

    public String getParserName() { return MERGE_PATTERNS; }

    /**
     * Parses a patterns element and returns a patterns object.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

	    PatternList patternList = (PatternList)xo.getChild(0);
	    Patterns patterns = new Patterns(patternList);
	    for (int i = 1; i < xo.getChildCount(); i++) {
		    patterns.addPatterns((PatternList)xo.getChild(i));
	    }

        return patterns;
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
        new ElementRule(PatternList.class, 1, Integer.MAX_VALUE)
    };

    public String getParserDescription() {
        return "A weighted list of the unique site patterns (unique columns) in an alignment.";
    }

    public Class getReturnType() { return PatternList.class; }

}
