/*
 * ConstantPatternsParser.java
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

import dr.evolution.alignment.*;
import dr.evolution.datatype.DataType;
import dr.evolution.util.TaxonList;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * Creates a set of patterns for constant sites with weights as provided
 * Can be merged with other pattern lists to pad out polymorphic sites
 *
 * @author Andrew Rambaut
 * @version $Id: $
 */
public class ConstantPatternsParser extends AbstractXMLObjectParser {

    public static final String CONSTANT_PATTERNS = "constantPatterns";
    private static final String COUNTS = "counts";

    public String getParserName() {
        return CONSTANT_PATTERNS;
    }

    /**
     * Parses a patterns element and returns a patterns object.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        PatternList source = (PatternList)xo.getChild(PatternList.class);
        Parameter constantPatternCounts = (Parameter) xo.getElementFirstChild(COUNTS);

        Patterns patterns = new Patterns(source.getDataType(), source);

        for (int i = 0; i < source.getDataType().getStateCount(); i++) {
            int[] pattern = new int[patterns.getPatternLength()];
            for (int j = 0; j < pattern.length; j++) {
                pattern[j] = i;
            }
            patterns.addPattern(pattern, constantPatternCounts.getParameterValue(i));
        }

        return patterns;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(PatternList.class),
            new ElementRule(COUNTS,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    })
    };

    public String getParserDescription() {
        return "Creates a set of patterns for constant sites with weights as provided.";
    }

    public Class getReturnType() {
        return PatternList.class;
    }

}
