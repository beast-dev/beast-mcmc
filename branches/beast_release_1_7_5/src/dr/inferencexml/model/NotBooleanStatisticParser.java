/*
 * NotBooleanStatisticParser.java
 *
 * Copyright (c) 2002-2012 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.model;

import dr.inference.model.BooleanStatistic;
import dr.inference.model.NotBooleanStatistic;
import dr.xml.*;

/**
 */
public class NotBooleanStatisticParser extends AbstractXMLObjectParser {

    public static String NEGATE_STATISTIC = "notBooleanStatistic";

    public String getParserName() {
        return NEGATE_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        BooleanStatistic booleanStatistic = (BooleanStatistic) xo.getChild(BooleanStatistic.class);
        return new NotBooleanStatistic(booleanStatistic);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a statistic that is the element-wise negation of the child boolean statistic.";
    }

    public Class getReturnType() {
        return BooleanStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(BooleanStatistic.class),
    };

}
