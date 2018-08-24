/*
 * MsatPatternStatisticParser.java
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

import dr.xml.*;
import dr.evolution.alignment.Patterns;
import dr.evolution.alignment.MicrosatellitePatternStatistic;

/**
 * @author Chieh-Hsi Wu
 * Parser for computing the statistics of msat pattern
 */
public class MicrosatellitePatternStatisticParser extends AbstractXMLObjectParser {
    public static final String MSAT_PATTERN_STATISTIC_PARSER = "msatPatternStatistic";
    public static final String MODE = "mode";

    public String getParserName(){
        return MSAT_PATTERN_STATISTIC_PARSER;
    }


    public Object parseXMLObject(XMLObject xo) throws XMLParseException {


        Patterns pats = (Patterns)xo.getChild(Patterns.class);
        if(xo.hasAttribute(MODE)){
            return new MicrosatellitePatternStatistic(pats, xo.getStringAttribute(MODE));
        }


        return new MicrosatellitePatternStatistic(pats);

    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(Patterns.class),
                AttributeRule.newStringRule(MODE, true)
        };
    }



    public String getParserDescription(){
        return "Returns MsatPatternStatistic object";
    }

    public Class getReturnType(){
        return MicrosatellitePatternStatistic.class;
    }
}
